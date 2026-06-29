"""Compile and run 3D-PDR, streaming the combined output line by line.

Kept deliberately small: each function yields log lines as they appear so the
UI can show a live console. Callers decide how to display/store them.
"""
from __future__ import annotations

import os
import queue
import subprocess
import threading
from collections.abc import Iterator
from pathlib import Path

from . import paths


def _stream(cmd: list[str], cwd: Path, env: dict | None = None) -> Iterator[str]:
    """Run ``cmd`` in ``cwd``, yielding combined stdout/stderr lines, then a
    final ``"__EXIT__ <code>"`` sentinel line so callers know the result."""
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=env,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line.rstrip("\n")
    proc.wait()
    yield f"__EXIT__ {proc.returncode}"


def compile_code(clean: bool = True) -> Iterator[str]:
    """Compile 3D-PDR in ``src/`` via the makefile. Streams build output.

    Uses ``make clean; make`` (semicolon, not ``&&``) on purpose: the makefile's
    ``clean`` target runs ``rm *.mod *.o ../restart.bin`` without ``-f``, so on an
    already-clean tree it exits non-zero. We must not let that harmless failure
    abort the build, while still surfacing a *real* build error (the final exit
    code is ``make``'s). Runs under ``bash -lc`` so the login profile's
    ``SUNDIALS_DIR`` (and compiler env) is available.
    """
    src = paths.src_dir()
    shell_cmd = "make clean; make" if clean else "make"
    yield f"$ cd {src} && {shell_cmd}"
    yield from _stream(["bash", "-lc", shell_cmd], cwd=src)


def run_model() -> Iterator[str]:
    """Execute the compiled ./3DPDR binary from the root. Streams its output."""
    root = paths.pdr_root()
    exe = paths.executable()
    if not exe.exists():
        yield f"ERROR: executable not found at {exe}. Compile the code first."
        yield "__EXIT__ 127"
        return
    yield f"$ cd {root} && ./3DPDR"
    yield from _stream(["./3DPDR"], cwd=root)


def rttool_needs_build() -> bool:
    """True if the RTtool binary is missing or older than any RT-tool source."""
    exe = paths.rt_executable()
    if not exe.exists():
        return True
    src_dir = paths.rt_tool_dir()
    if not src_dir.is_dir():
        return False
    exe_mtime = exe.stat().st_mtime
    return any(p.stat().st_mtime > exe_mtime
               for p in src_dir.glob("*.F90"))


def compile_rttool() -> Iterator[str]:
    """Compile RTtool via its makefile (builds and moves the binary to root)."""
    d = paths.rt_tool_dir()
    yield f"$ cd {d} && make"
    yield from _stream(["bash", "-lc", "make"], cwd=d)


def run_rttool() -> Iterator[str]:
    """Run ./RTtool from the root (reads paramsRT.dat). Streams its output."""
    root = paths.pdr_root()
    exe = paths.rt_executable()
    if not exe.exists():
        yield f"ERROR: RTtool not found at {exe}. Compile it first."
        yield "__EXIT__ 127"
        return
    yield f"$ cd {root} && ./RTtool"
    yield from _stream(["./RTtool"], cwd=root)


def parse_exit(line: str) -> int | None:
    """Return the exit code if ``line`` is the sentinel, else None."""
    if line.startswith("__EXIT__"):
        try:
            return int(line.split()[1])
        except (IndexError, ValueError):
            return -1
    return None


# ===========================================================================
# Stoppable background run
# ===========================================================================
class BackgroundProcess:
    """Run a command in the background so the UI stays responsive and the run can
    be stopped.

    A daemon thread drains the process's combined stdout/stderr into a queue; the
    Streamlit page polls :meth:`drain` (e.g. from an auto-refreshing fragment) to
    append new lines, calls :meth:`is_running` to drive the Run/Stop buttons, and
    :meth:`stop` to terminate it. Nothing here touches Streamlit, so it is safe to
    keep an instance in ``st.session_state`` across reruns.
    """

    def __init__(self, cmd: list[str], cwd: Path, env: dict | None = None):
        self.cmd = cmd
        self.cwd = cwd
        self.env = env
        self.proc: subprocess.Popen | None = None
        self.lines: list[str] = []
        self.returncode: int | None = None
        self._q: queue.Queue = queue.Queue()
        self._done = threading.Event()
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        self.lines.append(f"$ cd {self.cwd} && {' '.join(self.cmd)}")
        # new session/process group so stop() can signal the whole tree
        self.proc = subprocess.Popen(
            self.cmd, cwd=str(self.cwd),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, env=self.env, start_new_session=True,
        )
        self._thread = threading.Thread(target=self._reader, daemon=True)
        self._thread.start()

    def _reader(self) -> None:
        assert self.proc is not None and self.proc.stdout is not None
        for line in self.proc.stdout:
            self._q.put(line.rstrip("\n"))
        self.proc.wait()
        self.returncode = self.proc.returncode
        self._done.set()

    def drain(self) -> list[str]:
        """Move any newly-available output into ``self.lines``; return the new lines."""
        new: list[str] = []
        while True:
            try:
                new.append(self._q.get_nowait())
            except queue.Empty:
                break
        self.lines.extend(new)
        return new

    def is_running(self) -> bool:
        return self.proc is not None and not self._done.is_set()

    def stop(self) -> None:
        """Terminate the run (SIGTERM the process group, then SIGKILL if needed)."""
        if self.proc is None or self._done.is_set():
            return
        try:
            os.killpg(os.getpgid(self.proc.pid), 15)  # SIGTERM the group
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(os.getpgid(self.proc.pid), 9)  # SIGKILL
        except ProcessLookupError:
            pass
        self.lines.append("__stopped by user__")


def model_process() -> BackgroundProcess | None:
    """A stoppable BackgroundProcess for the compiled ./3DPDR, or None if missing."""
    exe = paths.executable()
    if not exe.exists():
        return None
    return BackgroundProcess(["./3DPDR"], cwd=paths.pdr_root())
