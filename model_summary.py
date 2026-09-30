#!/usr/bin/env python3
"""Summarise a 3D-PDR model: parameters and elemental abundances (X/H).

Usage: python3 model_summary.py PREFIX [--cell N] [--top K]
  PREFIX includes the directory (e.g. sims/Freeze3D). Reads PREFIX.species,
  PREFIX.params and ONE cell of PREFIX.pdr.h5 if present (needs h5py), otherwise
  ONE line of PREFIX.pdr.fin (streamed, the file is never loaded in full).

Files (written by src/writeoutputs.F90):
  .species : "index  name" per line (no abundances/masses)
  .params  : row1 = G0 [Draine], CR field choice 'L/H/U' (CRATTENUATION=1) or zeta [s^-1]
             (zeta*1.3e-17 as written), dust-to-gas/metallicity, v_turb [km/s]
             (row1 is empty for CRATTENUATION=2)
             row2 = nx ny nz ; row3 = box lengths [pc]
  .pdr.fin : 3D: idx x y z Tgas Tdust etype rho UVfield ab(1..nspec)   (9 leading cols)
             1D: idx x AV Tgas Tdust etype rho UVfield ab(1..nspec)   (8 leading cols)
"""
import sys, re, os, argparse
import numpy as np
from itertools import islice

ELEM2 = ['He', 'Mg', 'Si', 'Fe', 'Na', 'Cl', 'Ne', 'Ar', 'Li', 'Ca', 'Al', 'Ti', 'Zn', 'Br', 'Be', 'Cr', 'Mn', 'Ni', 'Co']
ELEM2 = [e for e in ELEM2 if e != 'Co']  # 'Co' ambiguous with C+o; never appears (CO is C+O)
ELEM1 = set('HDCNOSPFKB')  # D handled as deuterium (own element label 'D')
NONSPECIES = {'CRP', 'PHOTON', 'CRPHOT', 'FREEZE', 'DESORB', 'CRDESORB', 'PAH', 'PAH-', 'PAH+', 'PAH0',
              'G0', 'G-', 'G+', 'GRAIN', 'GRAIN-', 'GRAIN0', 'GRAIN+', 'ELECTR', 'E-', 'e-', 'e'}
GRAINLIKE = {'PAH', 'PAH-', 'PAH+', 'PAH0', 'G0', 'G-', 'G+', 'GRAIN', 'GRAIN-', 'GRAIN0', 'GRAIN+'}
ELECTRONS = {'e-', 'E-', 'ELECTR', 'e', 'ELECTRON'}


def parse_species(name):
    """Return (kind, {element: count}) ; kind in gas/ice/electron/grain; None comp if unparseable."""
    s = name.strip()
    if s in ELECTRONS:
        return 'electron', {}
    if s in GRAINLIKE:
        return 'grain', {}
    kind = 'gas'
    if s.startswith('#') or s.startswith('@'):
        kind, s = 'ice', s[1:]
    m = re.match(r'^(.*?)_[sS]$', s)          # H2O_s style ices
    if m:
        kind, s = 'ice', m.group(1)
    s = re.sub(r'[+-]+$', '', s)              # charges
    s = re.sub(r'^[op]-?(?=H2|D2|H3)', '', s)  # ortho/para labels (oH2, p-H2)
    comp = {}
    i = 0
    while i < len(s):
        c = s[i]
        # isotopes: leading digits before element (13C) -> treat as own label
        m = re.match(r'(\d+)([A-Z][a-z]?)', s[i:])
        if m and m.group(2) in ELEM2 + list(ELEM1):
            el = m.group(1) + m.group(2); i += m.end(); label = el
        elif s[i:i+2] in ELEM2:
            label = s[i:i+2]; i += 2
        elif c in ELEM1:
            label = c; i += 1
            m2 = re.match(r'(\d{2})(?!\d)', s[i:])   # suffix isotope: C13, O18
            if m2 and c in 'CNOS' and int(m2.group(1)) >= 13 and s[i:i+2] in ('13','15','17','18','34'):
                label = c + m2.group(1); i += 2
        else:
            return kind, None
        m = re.match(r'\d+', s[i:])
        n = 1
        if m:
            n = int(m.group()); i += m.end()
        comp[label] = comp.get(label, 0) + n
    if not comp:
        return kind, None
    return kind, comp


def read_species(path):
    names = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if len(t) >= 2 and t[0].isdigit():
                names.append(t[1])
            else:
                names.append(t[0])
    return names


def read_params(path):
    with open(path) as f:
        l = [f.readline() for _ in range(3)]
    r1, r2, r3 = l[0].split(), l[1].split(), l[2].split()
    p = {}
    if len(r1) < 4:   # CRATTENUATION=2 writes no first row: rows shift up by one
        p['G0'] = None
        r2, r3 = r1, r2
    else:
        p['G0'] = float(r1[0])
        if r1[1].isalpha():
            p['cr'] = 'attenuated CR model "%s" (column-dependent zeta, CRATTENUATION=1)' % r1[1]
            p['zeta'] = None
        else:
            p['cr'] = None
            p['zeta'] = float(r1[1])
        p['dtg'] = float(r1[2]); p['vturb'] = float(r1[3])
    p['res'] = tuple(int(x) for x in r2[:3])
    p['box'] = tuple(float(x) for x in r3[:3])
    return p


def read_cell_h5(path, nspec, cell):
    """Read one cell from PREFIX.pdr.h5 (written by fin2h5). Datasets: abundance001..N (one 3D array per
    species, order = .species), x,y,z,tgas,tdust,etype,rho,uv, av000..; h5py shape is (nz,ny,nx).
    .pdr.fin line order is x slowest / z fastest, so fin line N <-> h5 [iz,iy,ix] with
    (ix,iy,iz)=unravel(N-1,(nx,ny,nz)). Returns None if h5py is unavailable."""
    try:
        import h5py
    except ImportError:
        print('NOTE: h5py not available; falling back to .pdr.fin')
        return None
    with h5py.File(path, 'r') as f:
        nz, ny, nx = f['rho'].shape
        if cell < 1 or cell > nx * ny * nz:
            sys.exit('cell %d out of range 1..%d' % (cell, nx * ny * nz))
        ix, iy, iz = np.unravel_index(cell - 1, (nx, ny, nz))
        pos = (int(iz), int(iy), int(ix))
        nab = sum(1 for k in f if k.startswith('abundance'))
        if nab != nspec:
            sys.exit('%s has %d abundance datasets but .species has %d species' % (path, nab, nspec))
        ab = np.array([f['abundance%03d' % (i + 1)][pos] for i in range(nspec)])
        info = dict(x=f['x'][pos], y=f['y'][pos], z=f['z'][pos], Tgas=f['tgas'][pos], Tdust=f['tdust'][pos],
                    etype=int(f['etype'][pos]), rho=f['rho'][pos], UV=f['uv'][pos], idx=cell)
    return info, ab, 9


def read_cell(path, nspec, cell):
    with open(path) as f:
        line = next(islice(f, cell - 1, cell), None)
    if line is None:
        sys.exit('cell %d beyond end of %s' % (cell, path))
    t = line.split()
    lead = len(t) - nspec
    if lead not in (8, 9):
        sys.exit('pdr.fin has %d columns for %d species (leading=%d): unexpected layout' % (len(t), nspec, lead))
    vals = [float(x) for x in t]
    hdr = vals[:lead]
    if lead == 9:
        info = dict(x=hdr[1], y=hdr[2], z=hdr[3], Tgas=hdr[4], Tdust=hdr[5], etype=int(hdr[6]), rho=hdr[7], UV=hdr[8])
    else:
        info = dict(x=hdr[1], AV=hdr[2], Tgas=hdr[3], Tdust=hdr[4], etype=int(hdr[5]), rho=hdr[6], UV=hdr[7])
    info['idx'] = int(hdr[0])
    return info, np.array(vals[lead:]), lead


def fmt_top(contrib, k):
    tot = sum(v for _, v in contrib)
    contrib = sorted(contrib, key=lambda x: -x[1])[:k]
    return ', '.join('%s %.0f%%' % (n, 100 * v / tot) for n, v in contrib if tot > 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix')
    ap.add_argument('--cell', type=int, default=1,
                    help='1-based cell index, i.e. line number in .pdr.fin (default 1)')
    ap.add_argument('--top', type=int, default=3, help='number of main carriers listed per element (default 3)')
    a = ap.parse_args()
    pre = a.prefix
    for suf in ('.pdr.fin', '.species', '.params'):
        if pre.endswith(suf):
            pre = pre[:-len(suf)]
    names = read_species(pre + '.species')
    p = read_params(pre + '.params') if os.path.exists(pre + '.params') else None
    res = None
    src = pre + '.pdr.fin'
    if os.path.exists(pre + '.pdr.h5'):
        res = read_cell_h5(pre + '.pdr.h5', len(names), a.cell)
        if res is not None:
            src = pre + '.pdr.h5'
    if res is None:
        res = read_cell(pre + '.pdr.fin', len(names), a.cell)
    info, ab, lead = res

    print('=' * 78)
    print('3D-PDR model summary : %s' % pre)
    print('=' * 78)
    if p and p['G0'] is None:
        print('FUV/CR/Z/v_turb  : not recorded in .params (CRATTENUATION=2)')
    if p and p['G0'] is not None:
        print('FUV field G0     : %.4g  (Draine units)' % p['G0'])
        if p['zeta'] is not None:
            print('CR ionization    : %.4g s^-1 (zeta)' % p['zeta'])
        else:
            print('CR ionization    : %s' % p['cr'])
        print('Dust-to-gas      : %.4g  (x Solar; also scales A_V/N_H)' % p['dtg'])
        print('v_turb           : %.4g km/s' % p['vturb'])
    if p:
        print('Resolution       : %d x %d x %d cells%s' % (p['res'] + (('  (1D/no grid)' if p['res'][0] == 0 else ''),)))
        print('Box size         : %.4g x %.4g x %.4g  pc' % p['box'])
    print('Source           : %s' % src)
    print('Species          : %d in .species (%s layout)' % (len(names), '3D' if lead == 9 else '1D'))
    print('Cell used        : cell %d (idx %d): n=%.4g cm^-3, Tgas=%.4g K, Tdust=%.4g K, UV=%.4g' %
          (a.cell, info['idx'], info['rho'], info['Tgas'], info['Tdust'], info['UV']))
    print('-' * 78)

    gas, ice = {}, {}
    contrib_g, contrib_i = {}, {}
    unparsed, skipped = [], []
    xe = 0.0
    for n, x in zip(names, ab):
        kind, comp = parse_species(n)
        if kind == 'electron':
            xe += x; continue
        if kind == 'grain' or n in NONSPECIES:
            skipped.append(n); continue
        if comp is None:
            unparsed.append(n); continue
        tgt, ct = (ice, contrib_i) if kind == 'ice' else (gas, contrib_g)
        for el, c in comp.items():
            tgt[el] = tgt.get(el, 0.0) + c * x
            ct.setdefault(el, []).append((n, c * x))
    tot = {}
    for d in (gas, ice):
        for el, v in d.items():
            tot[el] = tot.get(el, 0.0) + v
    for w in unparsed:
        print('WARNING: could not parse species "%s" (excluded from totals)' % w)
    if skipped:
        print('Listed, not counted (grain/PAH-like): %s' % ', '.join(skipped))
    d = dict(zip(names, ab))
    for k in ('H2', 'H', 'CO', 'C+', 'e-'):
        if k in d:
            print('  x(%s) = %.4e' % (k, d[k]), end='')
    print('\n  net charge check: sum(ions)-sum(anions)-x_e = %.3e' % (
        sum(x * (n.count('+') - n.count('-')) for n, x in zip(names, ab) if not n.startswith('#') and n not in ELECTRONS) - xe))
    print('-' * 78)
    H = tot.get('H', 0.0)
    if H <= 0:
        sys.exit('total H is zero?')
    has_ice = bool(ice)
    print('Total H (nuclei per n_H, should be ~1) = %.6g   [gas %.6g%s]' % (
        H, gas.get('H', 0), ', ice %.4g' % ice.get('H', 0) if has_ice else ''))
    print('Elemental abundances relative to total H (gas+ice), with gas/ice split:')
    hdr = '%-5s %11s %8s   %11s %11s   %s' % ('El', 'X/H', '12+log', 'gas X/H', 'ice X/H' if has_ice else '', 'main carriers')
    print(hdr)
    for el in sorted(tot, key=lambda e: (e != 'H', -tot[e])):
        r = tot[el] / H
        lg = 12 + np.log10(r) if r > 0 else float('nan')
        car = fmt_top(contrib_g.get(el, []) + contrib_i.get(el, []), a.top)
        print('%-5s %11.4e %8.3f   %11.4e %11s   %s' % (el, r, lg, gas.get(el, 0) / H,
              ('%11.4e' % (ice.get(el, 0) / H)) if has_ice else '', car))
    print('=' * 78)


if __name__ == '__main__':
    main()
