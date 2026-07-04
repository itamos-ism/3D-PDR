#!/usr/bin/env python3
"""
3D-PDR RAM estimator.

Reads the coolant list and HEALPix level from a params.dat file, fetches
NLEV and NTEMP from each LAMDA file, and reports a per-component RAM
breakdown for a 3D model of user-specified resolution.

Components accounted for (based on initialization.F90):
  * Coolant level-population arrays   : line(N,N)+pop+solution+relchange+pophist  [per cell]
  * Ray-tracing path arrays           : epoint/length + projected + AV columns    [per cell]
  * Chemical abundances               : abundance(nspec)                           [per cell]
  * Other per-cell scalars            : Tgas, brackets, positions, flags, etc.
  * Static LAMDA data                 : A,B,freq matrices + collision tables       [once]
  * Temporary evalpop (peak)          : nrays*maxpoints*NLEV per thread            [during cooling]

Usage:
  python ram_estimate.py --res 128 [options]

Options:
  --params FILE    params.dat file  (default: params.dat)
  --res N          Grid resolution: N³ cells  (required for per-cell estimates)
  --healpix L      HEALPix level override  (default: read from params.dat)
  --lmax N         Cap max energy levels per coolant (simulates -lmax=N; 0 = no cap)
  --raytheia 0|1|2 Ray mode matching config.mk RAYTHEIA= value:
                              0=legacy (epoint per cell, maxpoints=1000)
                              1=RAYTHEIA_MO (global path arrays, recommended)
                              2=RAYTHEIA without MO (length per cell, maxpoints=3N-2)
                              (default: 0)
  --threads T      OpenMP thread count for peak evalpop estimate  (default: 1)
  --nspec S        Chemical species count  (default: auto-detect, else 33)
  --chemdir DIR    chemfiles directory  (default: chemfiles/ next to params.dat)
"""

import argparse
import math
import os
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------
DB = 8          # bytes per real(kind=dp)
DI = 4          # bytes per integer(kind=i4b)

# 3DPDR.F90 (right after "ITERATION = 0", before the main loop) allocates the
# per-cell cooling/heating arrays like this:
#     do i=1,coo
#       allocate(pdr(p)%cooling(i))     ! re-allocates a POINTER, size 1,2,...,coo
#     enddo
#     do i=1,12
#       allocate(pdr(p)%heating(i))     ! re-allocates a POINTER, size 1,2,...,12
#     enddo
# cooling/heating are POINTER (not ALLOCATABLE), so each re-ALLOCATE silently
# orphans the previous target instead of erroring or freeing it: every
# intermediate size (1..coo-1 and 1..11) is a permanent, unreachable leak for
# the life of the run, on top of the final correctly-sized array actually
# used. This should be a single `allocate(pdr(p)%cooling(coo))` / `(12)` per
# cell instead. Measured empirically (glibc malloc, gcc -O2, this host) by
# replaying the exact allocation pattern and reading VmRSS before/after:
#   buggy (19 malloc calls/cell, as in the code today) : ~1104 bytes/cell
#   fixed (2 malloc calls/cell, correctly sized)        : ~176 bytes/cell
# i.e. the leak alone wastes ~928 bytes/cell -- for a 256^3 grid that is
# ~14.5 GiB of real RAM for nothing. This is allocator/host-dependent, so
# treat it as an order-of-magnitude figure, not an exact one.
COOLHEAT_BYTES_PER_CELL_BUGGY = 1104
COOLHEAT_BYTES_PER_CELL_FIXED = 176

def detect_coolheat_bug_fixed(srcdir):
    """Inspect 3DPDR.F90 to see whether the cooling/heating re-ALLOCATE-on-
    POINTER leak (see COOLHEAT_BYTES_PER_CELL_* comment above) is still
    present. Returns (fixed: bool, note: str). Defaults to "buggy" (the
    conservative, higher-RAM assumption) if the source can't be inspected."""
    path = os.path.join(srcdir or '', '3DPDR.F90')
    if not os.path.isfile(path):
        return False, f"3DPDR.F90 not found under {srcdir}; assuming unfixed (conservative)"
    with open(path) as fh:
        src = fh.read()
    if 'allocate(pdr(p)%cooling(i))' in src or 'allocate(pdr(p)%heating(i))' in src:
        return False, f"leak pattern found in {path}"
    if 'allocate(pdr(p)%cooling(coo))' in src and 'allocate(pdr(p)%heating(12))' in src:
        return True, f"fixed allocation found in {path}"
    return False, f"could not recognise the allocation pattern in {path}; assuming unfixed (conservative)"

# config.mk NETWORK value -> chemfiles/{species,rates}_<suffix> filename suffix
# (mirrors the #ifdef chain in input_parameters.F90; note case sensitivity
#  matches the #elif branches there, e.g. 'myname' is lowercase)
NETWORK_SUFFIX = {
    'REDUCED': 'reduced.d', 'MEDIUM': 'medium.d', 'MEDIUMS': 'mediumS.d',
    'SCUTUM': 'scutum.d', 'FULL': 'full.d', 'MYNETWORK': 'mynetwork.d',
    'OURANIA': 'OURANIA.d', 'myname': 'myname.d', 'NEWNETWORK': 'NEWNETWORK.d',
    'NETWITHGAN': 'NETWITHGAN.d', 'network3': 'network3.d',
    'NEWNETWORKgroup': 'NEWNETWORKgroup.d', 'mycustomnet': 'mycustomnet.d',
}

# ---------------------------------------------------------------------------
# config.mk parser
# ---------------------------------------------------------------------------
def parse_config_mk(path):
    """Return dict of KEY -> VALUE (both str) from a config.mk file, ignoring
    comment lines (#) and blank lines. Returns {} if the file is missing."""
    flags = {}
    if not path or not os.path.isfile(path):
        return flags
    with open(path) as fh:
        for raw in fh:
            line = raw.split('#', 1)[0].strip()
            if not line or '=' not in line:
                continue
            key, val = line.split('=', 1)
            flags[key.strip()] = val.strip()
    return flags


def flag_on(config, key, default=True):
    """True if config[key] == '1'; if key absent, return default (so the
    estimate degrades gracefully to previous behaviour without a config.mk)."""
    if key not in config:
        return default
    return config[key] == '1'


def count_datfile_lines(path):
    """Count non-blank lines (mirrors Fortran list-directed read-to-EOF)."""
    with open(path) as fh:
        return sum(1 for l in fh if l.strip())

# ---------------------------------------------------------------------------
# LAMDA header reader
# ---------------------------------------------------------------------------
def read_lamda_header(path):
    """Return (nlev, ntemp, species_name) from a LAMDA file."""
    data_lines = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if line and not line.startswith('!'):
                data_lines.append(line)
            if len(data_lines) >= 3:
                break
    if len(data_lines) < 3:
        raise ValueError(f"Too few data lines in {path}")
    name = data_lines[0]
    # data_lines[1] = molecular weight
    parts = data_lines[2].split()
    nlev  = int(float(parts[0]))
    ntemp = int(float(parts[1])) if len(parts) > 1 else 1
    return nlev, ntemp, name

# ---------------------------------------------------------------------------
# params.dat parser
# ---------------------------------------------------------------------------
def parse_params(params_file):
    """
    Return dict with keys:
      healpix_level  (int)
      coolant_files  (list of bare filenames, e.g. ['12co.dat', ...])

    params.dat structure:
        =====         ← separator
        Section title
        =====         ← separator
        data line
        data line
        =====         ← separator (starts next section)
        ...
    State machine: track whether the previous non-blank content was a
    separator, so we can identify title lines vs data lines.
    """
    result = {'healpix_level': None, 'coolant_files': []}

    # Build a flat list of (kind, text): kind is 'sep' or 'data'.
    tokens = []
    with open(params_file) as fh:
        for raw in fh:
            s = raw.strip()
            if not s:
                continue
            if s.startswith('=') or s.startswith('-'):
                tokens.append(('sep', s))
            else:
                val = s.split('!')[0].strip()
                if val:
                    tokens.append(('data', val))

    # State machine: sep → title → sep → data ... → sep (next section)
    STATE_WANT_TITLE = 'want_title'
    STATE_WANT_SEP   = 'want_sep'
    STATE_IN_DATA    = 'in_data'

    state   = STATE_WANT_TITLE
    section = ''

    for kind, text in tokens:
        if kind == 'sep':
            if state == STATE_WANT_TITLE:
                pass                          # consecutive seps, stay
            elif state == STATE_WANT_SEP:
                state = STATE_IN_DATA         # sep after title → data region opens
            elif state == STATE_IN_DATA:
                state = STATE_WANT_TITLE      # sep closes data region
                section = ''
        else:
            if state == STATE_WANT_TITLE:
                section = text.lower()
                state   = STATE_WANT_SEP
            elif state == STATE_WANT_SEP:
                # title without closing sep → treat as data immediately
                state = STATE_IN_DATA
                _consume(result, section, text)
            elif state == STATE_IN_DATA:
                _consume(result, section, text)

    return result


def _consume(result, section, text):
    """Store a data token into result based on the current section title."""
    if 'ray' in section and 'trac' in section:
        if result['healpix_level'] is None:
            try:
                result['healpix_level'] = int(float(text))
            except ValueError:
                pass
    elif 'coolant' in section:
        if text.lower().endswith('.dat'):
            result['coolant_files'].append(text)

# ---------------------------------------------------------------------------
# auto-detect nspec from species file
# ---------------------------------------------------------------------------
def detect_nspec(chemdir):
    """Count lines in species_*.dat files under chemdir. Returns None if not found."""
    if not os.path.isdir(chemdir):
        return None
    for fname in os.listdir(chemdir):
        if fname.startswith('species_') and fname.endswith('.dat'):
            path = os.path.join(chemdir, fname)
            try:
                with open(path) as fh:
                    count = sum(1 for l in fh if l.strip() and not l.startswith('!'))
                if count > 0:
                    return count, fname
            except OSError:
                pass
    return None

# ---------------------------------------------------------------------------
# per-cell Fortran struct/descriptor-overhead probe
# ---------------------------------------------------------------------------
# pdr_node (one per grid cell) and pdr_excit (one per coolant per grid cell)
# hold several POINTER/ALLOCATABLE array components (epoint, AV, abundance,
# pop, line, pophist, ...). Each such component costs an array *descriptor*
# in the struct itself (base address, dtype, per-dimension stride/bound)
# in addition to the heap memory backing the array once allocated -- and the
# per-cell components below account only for the backing memory, not the
# descriptors. This probe compiles a throwaway program against the project's
# own modules.F90 (with the same preprocessor flags used to build the real
# binary) and asks the compiler for STORAGE_SIZE, so the overhead is measured
# rather than guessed and stays correct if the derived types change.
STRUCT_PROBE_SRC = """
program sizetest
  use maincode_module
  implicit none
  type(pdr_node), allocatable :: t(:)
  type(pdr_excit) :: e
  allocate(t(1))
  print '(A,I0)', "PDR_NODE_BYTES=", storage_size(t(1))/8
  print '(A,I0)', "PDR_EXCIT_BYTES=", storage_size(e)/8
end program
"""

def probe_struct_overhead(srcdir, config):
    """Compile modules.F90 from srcdir with gfortran to measure the real
    per-cell (pdr_node) and per-coolant-per-cell (pdr_excit) struct overhead.
    Returns (pdr_node_bytes, pdr_excit_bytes, note) or (None, None, reason)."""
    needed = ('definitions.F90', 'healpix_types.F90', 'modules.F90')
    if not srcdir or not all(os.path.isfile(os.path.join(srcdir, f)) for f in needed):
        return None, None, f"source files not found under {srcdir}"

    raytheia = int(config.get('RAYTHEIA', 0) or 0)
    defines = []
    if raytheia >= 1:
        defines.append('-DRAYTHEIA')
    if raytheia == 1:
        defines.append('-DRAYTHEIA_MO')
    for key in ('THERMALBALANCE', 'NGACCEL', 'NGRELAX', 'FREEZECOOL',
                'ILLINOIS', 'REBRACKET', 'RESTART'):
        if flag_on(config, key, default=True):
            defines.append(f'-D{key}')

    with tempfile.TemporaryDirectory(prefix='ram_estimate_probe_') as tmp:
        try:
            for f in needed:
                with open(os.path.join(tmp, f)) as _:
                    pass
        except FileNotFoundError:
            pass
        for f in needed:
            with open(os.path.join(srcdir, f)) as src, \
                 open(os.path.join(tmp, f), 'w') as dst:
                dst.write(src.read())
        probe_path = os.path.join(tmp, 'sizetest.F90')
        with open(probe_path, 'w') as fh:
            fh.write(STRUCT_PROBE_SRC)

        compile_cmd = (['gfortran', '-cpp'] + defines +
                       ['-c', 'definitions.F90', 'healpix_types.F90', 'modules.F90'])
        try:
            r = subprocess.run(compile_cmd, cwd=tmp, capture_output=True, text=True, timeout=60)
            if r.returncode != 0:
                return None, None, f"gfortran compile failed: {r.stderr.strip()[-300:]}"
            r = subprocess.run(['gfortran', '-cpp'] + defines +
                               ['sizetest.F90', 'definitions.o', 'healpix_types.o',
                                'modules.o', '-o', 'sizetest'],
                               cwd=tmp, capture_output=True, text=True, timeout=60)
            if r.returncode != 0:
                return None, None, f"gfortran link failed: {r.stderr.strip()[-300:]}"
            r = subprocess.run(['./sizetest'], cwd=tmp, capture_output=True, text=True, timeout=10)
        except (OSError, subprocess.SubprocessError) as ex:
            return None, None, f"gfortran unavailable: {ex}"

        node_bytes = excit_bytes = None
        for line in r.stdout.splitlines():
            if line.startswith('PDR_NODE_BYTES='):
                node_bytes = int(line.split('=', 1)[1])
            elif line.startswith('PDR_EXCIT_BYTES='):
                excit_bytes = int(line.split('=', 1)[1])
        if node_bytes is None or excit_bytes is None:
            return None, None, f"unexpected probe output: {r.stdout!r}"
        flag_str = ' '.join(defines) or '(none)'
        return node_bytes, excit_bytes, f"measured via gfortran ({flag_str})"

# ---------------------------------------------------------------------------
# formatting helpers
# ---------------------------------------------------------------------------
def fmt_bytes(n):
    """Human-readable byte count."""
    for unit in ('B', 'kB', 'MB', 'GB', 'TB'):
        if abs(n) < 1024.0:
            return f"{n:.2f} {unit}"
        n /= 1024.0
    return f"{n:.2f} PB"

def bar(frac, width=30):
    filled = int(round(frac * width))
    return '[' + '#' * filled + '.' * (width - filled) + ']'

# ---------------------------------------------------------------------------
# RAM computation
# ---------------------------------------------------------------------------
def compute_ram(coolants, res, healpix_level, lmax, raytheia, nthreads, nspec,
                struct_overhead=(None, None), chemanalysis_nreac=None,
                coolheat_bug_fixed=False):
    """
    coolants : list of (nlev, ntemp, name, path)
    res      : grid resolution per side (N → N³ cells)
    struct_overhead : (pdr_node_bytes, pdr_excit_bytes) measured via
                       probe_struct_overhead(), or (None, None) to fall back
                       to a rough fixed guess for pdr_node and 0 for pdr_excit.
    chemanalysis_nreac : if not None, adds the CHEMANALYSIS temp_rate(nreac,
                          pdr_ptot) global array (see allocations.F90).
    Returns dict of component → bytes
    """
    ncells  = res ** 3
    nside   = 2 ** healpix_level
    nrays   = 12 * nside ** 2

    # Effective levels per coolant after lmax cap
    eff = []
    for nlev, ntemp, name, path in coolants:
        n = min(nlev, lmax) if lmax > 0 else nlev
        eff.append((n, ntemp, name, path, nlev))

    # ------------------------------------------------------------------ #
    #  maxpoints                                                           #
    # ------------------------------------------------------------------ #
    # config.mk RAYTHEIA=1 → -DRAYTHEIA -DRAYTHEIA_MO  (global arrays)
    # config.mk RAYTHEIA=2 → -DRAYTHEIA                 (per-cell length)
    # config.mk RAYTHEIA=0 →  neither                   (epoint per cell)
    if raytheia == 0:
        maxpoints = 1000           # hardcoded in readdensity.F90 line 89
        mp_note   = "hardcoded (RAYTHEIA=0 in config.mk)"
    elif raytheia == 1:            # RAYTHEIA_MO: computed from octree traversal
        maxpoints = int(math.ceil(math.sqrt(3) * res))
        mp_note   = f"~{maxpoints} estimated (RAYTHEIA_MO, octree; actual printed at run time)"
    else:                          # raytheia == 2: 3N-2 formula, readdensity.F90 line 87
        maxpoints = 3 * res - 2
        mp_note   = f"3×{res}−2 = {maxpoints} (RAYTHEIA=2 in config.mk)"

    # ------------------------------------------------------------------ #
    #  COMPONENT 1: coolant level-population arrays (per cell, dominant)  #
    # ------------------------------------------------------------------ #
    # initialization.F90 lines 33-46:
    #   pop(N), line(N,N), solution(N), relativechange(N)   [always]
    #   pophist(N,3)                                         [NGACCEL=1]
    #   coolprev (scalar)                                    [FREEZECOOL=1]
    # → N² + 6N + 1  doubles per cell per coolant
    cool_pop_bytes_per_cell = 0
    cool_pop_detail = []
    for n, ntemp, name, path, nlev_orig in eff:
        doubles = n*n + 6*n + 1
        b = doubles * DB
        cool_pop_bytes_per_cell += b
        cool_pop_detail.append((name, nlev_orig, n, doubles, b))
    cool_pop_total = cool_pop_bytes_per_cell * ncells

    # ------------------------------------------------------------------ #
    #  COMPONENT 2: ray-tracing path arrays (per cell)                    #
    # ------------------------------------------------------------------ #
    # initialization.F90 lines 48-65 (config.mk → preprocessor flags):
    #   RAYTHEIA=0 → neither -DRAYTHEIA nor -DRAYTHEIA_MO:
    #     epray(nrays)                         → nrays dp        per cell
    #     epoint(3, nrays, maxpoints+1)        → 3*nrays*(mp+1) integers  per cell
    #     projected(nrays, maxpoints+1)        → nrays*(mp+1) dp per cell
    #   RAYTHEIA=1 → -DRAYTHEIA + -DRAYTHEIA_MO  (memory-optimised, recommended):
    #     epray, plength, projected stored GLOBALLY once (not per cell)
    #   RAYTHEIA=2 → -DRAYTHEIA only (no MO flag):
    #     epray(nrays)                         → nrays dp        per cell
    #     length(nrays, maxpoints+1)           → nrays*(mp+1) dp per cell
    #     projected(nrays, maxpoints+1)        → nrays*(mp+1) dp per cell
    #   Always per cell regardless of mode:
    #     raytype(nrays)                       → nrays integers
    #     AV, rad_surface, NH2, NHD, NCO, NC, NS  (7 arrays × nrays dp)

    ray_scalar_bytes = (7 * nrays * DB           # AV + rad_surface + 5 columns
                        + nrays * DI)            # raytype

    if raytheia == 0:
        # config.mk RAYTHEIA=0: per-cell epoint + projected
        # epoint is declared real(kind=dp) in modules.F90 (pdr_node%epoint),
        # not integer, despite holding coordinate/index-like data.
        ray_path_bytes = (nrays * DB                          # epray
                          + 3 * nrays * (maxpoints+1) * DB    # epoint (real(dp))
                          + nrays * (maxpoints+1) * DB)       # projected
    elif raytheia == 1:
        # config.mk RAYTHEIA=1 = RAYTHEIA_MO: path arrays stored globally, not per cell
        ray_path_bytes = 0
    else:
        # config.mk RAYTHEIA=2: per-cell length (doubles) + projected
        ray_path_bytes = (nrays * DB                          # epray
                          + 2 * nrays * (maxpoints+1) * DB)  # length + projected

    ray_bytes_per_cell = ray_scalar_bytes + ray_path_bytes
    ray_total = ray_bytes_per_cell * ncells

    # RAYTHEIA_MO (raytheia=1) global path arrays: allocated once, not per cell
    if raytheia == 1:
        ray_global = (nrays * DB                          # epray
                      + 2 * nrays * (maxpoints+1) * DB)  # plength + projected
    else:
        ray_global = 0

    # ------------------------------------------------------------------ #
    #  COMPONENT 3: chemical abundances (per cell)                        #
    # ------------------------------------------------------------------ #
    chem_bytes_per_cell = nspec * DB
    chem_total = chem_bytes_per_cell * ncells

    # ------------------------------------------------------------------ #
    #  COMPONENT 4: other per-cell scalars + Fortran struct/descriptor     #
    #  overhead (pdr_node itself, and pdr_excit per coolant per cell)      #
    # ------------------------------------------------------------------ #
    # Every POINTER/ALLOCATABLE array field inside pdr_node (epoint, AV,
    # abundance, cooling, heating, column_N*, ...) and pdr_excit (pop, line,
    # solution, relativechange, pophist) costs an array descriptor in the
    # struct itself, on top of the backing heap memory counted elsewhere in
    # this function. probe_struct_overhead() measures this exactly; if it
    # could not run (no gfortran / sources not found), fall back to the old
    # rough guess for pdr_node scalars and 0 (i.e. underestimate) for the
    # pdr_excit descriptor overhead.
    node_bytes, excit_bytes = struct_overhead
    if node_bytes is None:
        node_bytes = 12 * DB + 4 * DI   # rough fallback: ~12 dp + 4 int scalars
    if excit_bytes is None:
        excit_bytes = 0                 # rough fallback: descriptor overhead ignored

    n_coolants = len(eff)
    other_total = node_bytes * ncells
    struct_overhead_total = excit_bytes * n_coolants * ncells

    # ------------------------------------------------------------------ #
    #  COMPONENT 5: static LAMDA data (once, not per cell)               #
    # ------------------------------------------------------------------ #
    # energies(N), weights(N): 2N dp
    # A_COEFFS(N,N), B_COEFFS(N,N), frequencies(N,N): 3N² dp
    # temperatures(7, NTEMP): 7*NTEMP dp
    # H_COL+HP_COL+EL_COL+HE_COL+H2_COL+PH2_COL+OH2_COL: 7*N²*NTEMP dp
    lamda_total = 0
    lamda_detail = []
    for n, ntemp, name, path, nlev_orig in eff:
        doubles = 2*n + 3*n*n + 7*ntemp + 7*n*n*ntemp
        b = doubles * DB
        lamda_total += b
        lamda_detail.append((name, n, ntemp, doubles, b))

    # ------------------------------------------------------------------ #
    #  COMPONENT 5b: CHEMANALYSIS temp_rate(nreac, pdr_ptot) (once)        #
    # ------------------------------------------------------------------ #
    # allocations.F90: allocate(temp_rate(1:nreac,1:pdr_ptot)) under
    # #ifdef CHEMANALYSIS -- every reaction rate for every cell, kept for
    # the whole run. Not per-cell-in-pdr_node, but scales with ncells and is
    # frequently the single largest allocation once CHEMANALYSIS=1 in
    # config.mk (it says "not recommended for 3D models" for this reason).
    chemanalysis_total = 0
    if chemanalysis_nreac:
        chemanalysis_total = chemanalysis_nreac * ncells * DB

    # ------------------------------------------------------------------ #
    #  COMPONENT 5c: cooling(:)/heating(:) arrays, incl. the             #
    #  re-ALLOCATE-on-POINTER leak in 3DPDR.F90 (see constant comment)    #
    # ------------------------------------------------------------------ #
    coolheat_bytes_per_cell = (COOLHEAT_BYTES_PER_CELL_FIXED if coolheat_bug_fixed
                               else COOLHEAT_BYTES_PER_CELL_BUGGY)
    coolheat_total = coolheat_bytes_per_cell * ncells

    # ------------------------------------------------------------------ #
    #  COMPONENT 6: temporary evalpop (peak, during cooling evaluation)   #
    # ------------------------------------------------------------------ #
    # coolingfunctions.F90: evalpop(0:nrays-1, 0:maxpoints, 1:N) per coolant
    # PRIVATE per OMP thread → n_threads copies exist simultaneously
    evalpop_per_thread = sum((nrays * (maxpoints+1) * n) * DB
                             for n, ntemp, name, path, nlev_orig in eff)
    evalpop_total = evalpop_per_thread * nthreads

    # ------------------------------------------------------------------ #
    #  Totals                                                             #
    # ------------------------------------------------------------------ #
    grand_persistent = (cool_pop_total + ray_total + chem_total
                        + other_total + lamda_total + ray_global
                        + struct_overhead_total + chemanalysis_total
                        + coolheat_total)

    return {
        'ncells'         : ncells,
        'nrays'          : nrays,
        'maxpoints'      : maxpoints,
        'mp_note'        : mp_note,
        'cool_pop_detail': cool_pop_detail,
        'cool_pop_total' : cool_pop_total,
        'ray_bytes_per_cell': ray_bytes_per_cell,
        'ray_path_bytes' : ray_path_bytes,
        'ray_scalar_bytes': ray_scalar_bytes,
        'ray_total'      : ray_total,
        'ray_global'     : ray_global,
        'chem_total'     : chem_total,
        'other_total'    : other_total,
        'struct_overhead_total': struct_overhead_total,
        'node_bytes'     : node_bytes,
        'excit_bytes'    : excit_bytes,
        'n_coolants'     : n_coolants,
        'chemanalysis_total': chemanalysis_total,
        'coolheat_total' : coolheat_total,
        'coolheat_bytes_per_cell': coolheat_bytes_per_cell,
        'lamda_total'    : lamda_total,
        'lamda_detail'   : lamda_detail,
        'evalpop_total'  : evalpop_total,
        'evalpop_per_thread': evalpop_per_thread,
        'grand_persistent': grand_persistent,
    }

# ---------------------------------------------------------------------------
# HEALPix comparison table
# ---------------------------------------------------------------------------
def healpix_table(coolants, res, lmax, raytheia, nthreads, nspec, max_level=4,
                  struct_overhead=(None, None), chemanalysis_nreac=None,
                  coolheat_bug_fixed=False):
    rows = []
    for lv in range(max_level + 1):
        r = compute_ram(coolants, res, lv, lmax, raytheia, nthreads, nspec,
                        struct_overhead, chemanalysis_nreac, coolheat_bug_fixed)
        rows.append((lv, 12 * 4**lv, r['grand_persistent'], r['ray_total'],
                     r['cool_pop_total'], r['evalpop_total']))
    return rows

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Estimate 3D-PDR RAM consumption from params.dat",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    ap.add_argument('--params',   default='params.dat',
                    help='params.dat file (default: params.dat)')
    ap.add_argument('--res',      type=int, required=True,
                    help='Grid resolution per side, N (model has N³ cells)')
    ap.add_argument('--healpix',  type=int, default=None,
                    help='HEALPix level override (default: read from params.dat)')
    ap.add_argument('--lmax',     type=int, default=0,
                    help='Cap max energy levels at N (0 = no cap, default: 0)')
    ap.add_argument('--raytheia', type=int, default=None, choices=[0, 1, 2],
                    help='config.mk RAYTHEIA value: 0=legacy, 1=RAYTHEIA_MO (recommended), '
                         '2=RAYTHEIA without MO (default: read from config.mk, else 0)')
    ap.add_argument('--threads',  type=int, default=1,
                    help='OpenMP threads for peak evalpop estimate (default: 1)')
    ap.add_argument('--nspec',    type=int, default=None,
                    help='Chemical species count (default: derived from config.mk NETWORK)')
    ap.add_argument('--chemdir',  default=None,
                    help='Path to chemfiles/ directory (default: auto-detect from params location)')
    ap.add_argument('--srcdir',   default=None,
                    help='Path to src/ directory, for config.mk and the struct-overhead probe '
                         '(default: <params dir>/src)')
    ap.add_argument('--configmk', default=None,
                    help='Path to config.mk (default: <srcdir>/config.mk)')
    ap.add_argument('--chemanalysis', type=int, default=None, choices=[0, 1],
                    help='Override CHEMANALYSIS flag (default: read from config.mk). '
                         'When 1, adds the temp_rate(nreac,pdr_ptot) global array.')
    ap.add_argument('--no-probe', action='store_true',
                    help="Don't compile a struct-overhead probe with gfortran; use the "
                         "rough fallback estimate instead (faster, less accurate)")
    ap.add_argument('--coolheat-fixed', type=int, default=None, choices=[0, 1],
                    help='Override whether the cooling(:)/heating(:) re-ALLOCATE-on-POINTER '
                         'leak in 3DPDR.F90 (~line 174) is present (default: auto-detect by '
                         'inspecting 3DPDR.F90 under --srcdir)')
    args = ap.parse_args()

    # ---- locate files ----
    params_path = os.path.abspath(args.params)
    if not os.path.isfile(params_path):
        sys.exit(f"ERROR: params file not found: {params_path}")

    run_dir  = os.path.dirname(params_path)
    chemdir  = args.chemdir if args.chemdir else os.path.join(run_dir, 'chemfiles')
    srcdir   = args.srcdir if args.srcdir else os.path.join(run_dir, 'src')
    configmk_path = args.configmk if args.configmk else os.path.join(srcdir, 'config.mk')

    # ---- parse config.mk ----
    config = parse_config_mk(configmk_path)
    if not config:
        print(f"WARNING: could not read config.mk at {configmk_path}; "
              f"RAYTHEIA/CHEMANALYSIS/struct-overhead defaults will be used.", file=sys.stderr)

    # ---- parse params.dat ----
    info = parse_params(params_path)
    healpix_level = args.healpix if args.healpix is not None else info['healpix_level']
    if healpix_level is None:
        print("WARNING: Could not read HEALPix level from params.dat; assuming 0.", file=sys.stderr)
        healpix_level = 0

    # ---- resolve RAYTHEIA mode ----
    if args.raytheia is not None:
        raytheia = args.raytheia
        raytheia_src = "user-supplied"
    elif 'RAYTHEIA' in config:
        raytheia = int(config['RAYTHEIA'])
        raytheia_src = f"from {configmk_path}"
    else:
        raytheia = 0
        raytheia_src = "default (no config.mk found)"

    # ---- nspec / nreac from config.mk NETWORK (fixes: species_*.dat never
    #      matched actual species_*.d filenames, so this always silently
    #      fell back to a hardcoded 33) ----
    network = config.get('NETWORK')
    suffix = NETWORK_SUFFIX.get(network) if network else None
    nspec = args.nspec
    nspec_src = "user-supplied"
    nreac = None
    nreac_src = None
    if suffix:
        species_path = os.path.join(chemdir, f'species_{suffix}')
        rates_path   = os.path.join(chemdir, f'rates_{suffix}')
        if nspec is None and os.path.isfile(species_path):
            nspec = count_datfile_lines(species_path)
            nspec_src = f"NETWORK={network} -> {os.path.basename(species_path)}"
        if os.path.isfile(rates_path):
            nreac = count_datfile_lines(rates_path)
            nreac_src = f"NETWORK={network} -> {os.path.basename(rates_path)}"
    if nspec is None:
        det = detect_nspec(chemdir)
        if det:
            nspec, nspec_fname = det
            nspec_src = f"auto-detected from {nspec_fname}"
        else:
            nspec = 33
            nspec_src = "default (REDUCED network, no config.mk/NETWORK match)"

    # ---- CHEMANALYSIS: adds temp_rate(nreac,pdr_ptot), see allocations.F90 ----
    if args.chemanalysis is not None:
        chemanalysis_on = bool(args.chemanalysis)
        chemanalysis_src = "user-supplied"
    else:
        chemanalysis_on = flag_on(config, 'CHEMANALYSIS', default=False)
        chemanalysis_src = f"from {configmk_path}" if 'CHEMANALYSIS' in config else "default (off)"
    chemanalysis_nreac = nreac if (chemanalysis_on and nreac) else None
    if chemanalysis_on and not nreac:
        print("WARNING: CHEMANALYSIS is on but nreac could not be determined "
              "(unknown NETWORK / rates file not found); temp_rate array omitted.",
              file=sys.stderr)

    # ---- struct/descriptor overhead probe ----
    if args.no_probe:
        struct_overhead = (None, None)
        struct_note = "skipped (--no-probe)"
    else:
        probe_config = dict(config)
        probe_config['RAYTHEIA'] = str(raytheia)
        node_bytes, excit_bytes, struct_note = probe_struct_overhead(srcdir, probe_config)
        struct_overhead = (node_bytes, excit_bytes)
        if node_bytes is None:
            print(f"WARNING: struct-overhead probe unavailable ({struct_note}); "
                  f"pdr_node/pdr_excit descriptor overhead will be underestimated. "
                  f"Install gfortran or check --srcdir.", file=sys.stderr)

    # ---- cooling(:)/heating(:) leak: auto-detect unless overridden ----
    if args.coolheat_fixed is not None:
        coolheat_fixed = bool(args.coolheat_fixed)
        coolheat_src = "user-supplied"
    else:
        coolheat_fixed, coolheat_src = detect_coolheat_bug_fixed(srcdir)

    # ---- load coolant LAMDA headers ----
    coolants = []
    missing  = []
    for bare in info['coolant_files']:
        path = os.path.join(chemdir, bare)
        if not os.path.isfile(path):
            # try without chemfiles/ prefix in case the path is already absolute
            path2 = os.path.join(run_dir, bare)
            if os.path.isfile(path2):
                path = path2
            else:
                missing.append(bare)
                continue
        try:
            nlev, ntemp, name = read_lamda_header(path)
            coolants.append((nlev, ntemp, name.strip(), path, bare))
        except Exception as ex:
            print(f"WARNING: could not read {path}: {ex}", file=sys.stderr)

    if missing:
        print(f"WARNING: LAMDA files not found: {missing}", file=sys.stderr)

    if not coolants:
        sys.exit("ERROR: No coolant LAMDA files could be read.")

    # ---- reshape to (nlev, ntemp, name, path) without bare filename ----
    coolants_4 = [(nlev, ntemp, name, path) for nlev, ntemp, name, path, bare in coolants]

    # ---- run main estimate ----
    r = compute_ram(coolants_4, args.res, healpix_level, args.lmax,
                    raytheia, args.threads, nspec,
                    struct_overhead, chemanalysis_nreac, coolheat_fixed)

    RAYTHEIA_LABELS = {0: "legacy (no grid; epoint per cell; maxpoints=1000)",
                       1: "RAYTHEIA_MO (config.mk RAYTHEIA=1; path arrays global; recommended)",
                       2: "RAYTHEIA no-MO (config.mk RAYTHEIA=2; length per cell; maxpoints=3N−2)"}

    # ------------------------------------------------------------------ #
    #  Print report                                                       #
    # ------------------------------------------------------------------ #
    W = 72
    print('=' * W)
    print(f"  3D-PDR RAM Estimator")
    print('=' * W)
    print(f"  params.dat   : {params_path}")
    print(f"  chemfiles/   : {chemdir}")
    print(f"  config.mk    : {configmk_path}{'' if config else ' (not found)'}")
    print(f"  Resolution   : {args.res}³ = {r['ncells']:,} cells")
    print(f"  HEALPix level: {healpix_level}  →  {r['nrays']} rays per cell")
    print(f"  RAYTHEIA mode: {raytheia}  ({RAYTHEIA_LABELS[raytheia]})  [{raytheia_src}]")
    print(f"  maxpoints    : {r['maxpoints']}  ({r['mp_note']})")
    print(f"  OMP threads  : {args.threads}  (for peak evalpop)")
    print(f"  nspec        : {nspec}  ({nspec_src})")
    print(f"  nreac        : {nreac if nreac is not None else 'unknown'}  ({nreac_src or 'no NETWORK match'})")
    print(f"  CHEMANALYSIS : {'on' if chemanalysis_on else 'off'}  ({chemanalysis_src})")
    print(f"  struct probe : {struct_note}")
    print(f"  cooling/heat : {'fixed' if coolheat_fixed else 'BUGGY (leak)'}  ({coolheat_src})")
    lmax_str = f"-lmax={args.lmax} active" if args.lmax > 0 else "off (all levels used)"
    print(f"  LEVMAX       : {lmax_str}")
    print()

    # ---- coolant table ----
    print('-' * W)
    print(f"  Coolant files  ({len(coolants_4)} total)")
    print('-' * W)
    hdr = f"  {'Name':<12} {'File levels':>12} {'Used levels':>12} {'Cell arrays':>14} {'Cell kB':>9}"
    print(hdr)
    print(f"  {'-'*12} {'-'*12} {'-'*12} {'-'*14} {'-'*9}")
    for name, nlev_orig, n_used, doubles, b in r['cool_pop_detail']:
        trunc = f"→ {n_used}" if n_used < nlev_orig else "    "
        print(f"  {name:<12} {nlev_orig:>12} {n_used:>10} {trunc}  "
              f"{doubles:>14,}  {b/1024:>8.2f}")
    print(f"  {'TOTAL':<12} {'':>12} {'':>12} {r['cool_pop_total']//r['ncells']:>14,}  "
          f"{r['cool_pop_total']/r['ncells']/1024:>8.2f}")
    print()

    # ---- per-component breakdown ----
    total = r['grand_persistent']
    components = [
        ("Coolant pop arrays",   r['cool_pop_total'],   "(line(N,N)+pop+solution+relchange+pophist)"),
        ("Ray-tracing per cell", r['ray_total'],         f"(epoint/length+projected+AV+columns; {r['nrays']} rays)"),
        ("Chemistry abundances", r['chem_total'],        f"({nspec} species per cell)"),
        ("Other per-cell data",  r['other_total'],       f"(pdr_node struct: {r['node_bytes']} B/cell)"),
        ("Coolant struct overhead", r['struct_overhead_total'],
         f"(pdr_excit descriptors: {r['excit_bytes']} B x {r['n_coolants']} coolants/cell)"),
        ("Static LAMDA data",    r['lamda_total'],       "(A,B,freq,collision tables; loaded once)"),
    ]
    if r['ray_global']:
        components.append(("Ray arrays (global)", r['ray_global'],
                           "(RAYTHEIA_MO: epray+plength+projected once)"))
    if r['chemanalysis_total']:
        components.append(("CHEMANALYSIS temp_rate", r['chemanalysis_total'],
                           f"(temp_rate({nreac},pdr_ptot); not recommended for 3D)"))
    coolheat_display_note = (f"({r['coolheat_bytes_per_cell']} B/cell; 3DPDR.F90 leak fixed)"
                             if coolheat_fixed else
                             f"({r['coolheat_bytes_per_cell']} B/cell; BUGGY re-alloc-on-POINTER "
                             f"leak in 3DPDR.F90 ~line 174 -- see --coolheat-fixed)")
    components.append(("Cooling/heating arrays", r['coolheat_total'], coolheat_display_note))

    print('-' * W)
    print(f"  Persistent RAM breakdown  (all {args.res}³ cells)")
    print('-' * W)
    w_name = 24
    for label, val, note in components:
        frac = val / total if total > 0 else 0
        print(f"  {label:<{w_name}} {fmt_bytes(val):>10}  {bar(frac)}  {frac*100:4.1f}%")
        print(f"  {' '*w_name}   {note}")

    print('-' * W)
    print(f"  {'TOTAL persistent RAM':<{w_name}} {fmt_bytes(total):>10}")
    print()

    # ---- peak temporary RAM ----
    print('-' * W)
    print(f"  Peak temporary RAM (during cooling-function evaluation)")
    print('-' * W)
    print(f"  evalpop(nrays, maxpoints, NLEV) × {args.threads} thread(s)")
    print(f"  = {fmt_bytes(r['evalpop_per_thread'])} per thread  →  "
          f"{fmt_bytes(r['evalpop_total'])} peak")
    print()

    # ---- grand total ----
    grand = r['grand_persistent'] + r['evalpop_total']
    print('=' * W)
    print(f"  ESTIMATED TOTAL PEAK  :  {fmt_bytes(grand)}")
    print(f"  (persistent {fmt_bytes(r['grand_persistent'])}  +  "
          f"peak temporary {fmt_bytes(r['evalpop_total'])})")
    print()
    print(f"  Note: excludes OS, MPI buffers, Fortran runtime overhead,")
    print(f"  output file I/O buffers, and CVODE solver workspace (~few GB).")
    print('=' * W)
    print()

    # ---- HEALPix level comparison table ----
    print('-' * W)
    print(f"  HEALPix level comparison  (res={args.res}, RAYTHEIA={raytheia},"
          f" lmax={'off' if args.lmax==0 else args.lmax}, threads={args.threads})")
    print('-' * W)
    print(f"  {'Level':>6}  {'Rays':>6}  {'Total RAM':>12}  {'ΔvsL0':>10}  "
          f"{'Ray arrays':>12}  {'Coolant pop':>12}  {'evalpop':>10}")
    print(f"  {'-'*6}  {'-'*6}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*10}")

    hrows = healpix_table(coolants_4, args.res, args.lmax, raytheia,
                          args.threads, nspec, max_level=5,
                          struct_overhead=struct_overhead,
                          chemanalysis_nreac=chemanalysis_nreac,
                          coolheat_bug_fixed=coolheat_fixed)
    base_total = hrows[0][2] + hrows[0][5]   # grand total at level 0
    for lv, nray, persistent, ray_tot, cool_tot, evalpop_tot in hrows:
        grand_lv = persistent + evalpop_tot
        delta = grand_lv / base_total if base_total > 0 else 1.0
        marker = " ←" if lv == healpix_level else ""
        print(f"  {lv:>6}  {nray:>6}  {fmt_bytes(grand_lv):>12}  "
              f"{delta:>9.2f}×  {fmt_bytes(ray_tot):>12}  "
              f"{fmt_bytes(cool_tot):>12}  {fmt_bytes(evalpop_tot):>10}{marker}")

    print()
    print(f"  Coolant pop is independent of HEALPix level.")
    print(f"  Ray arrays and evalpop both scale as 12×4^L per cell.")
    print('=' * W)

if __name__ == '__main__':
    main()
