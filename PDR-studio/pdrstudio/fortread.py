"""Robust reading of 3D-PDR / RT-tool ASCII output tables.

Fortran's ``E``/``ES`` edit descriptors drop the ``E`` when the exponent needs
three digits (|exp| >= 100), writing e.g. ``1.0805662-164`` instead of
``1.0805662E-164``. Fortran reads this back fine, but numpy/Python ``float()``
cannot. We repair those tokens before parsing.
"""
from __future__ import annotations

import io
import re
from pathlib import Path

import numpy as np

# A digit, then a sign and exactly three exponent digits, with no 'E' between
# (a properly written exponent has 'E' before the sign, so '(\d)' won't match).
_BADEXP = re.compile(r"(\d)([+-]\d{3})(?![0-9])")


def fix_fortran_exponents(text: str) -> str:
    """Insert the missing ``E`` into dropped-exponent tokens."""
    return _BADEXP.sub(r"\1E\2", text)


def loadtxt(path, **kwargs) -> np.ndarray:
    """``numpy.loadtxt`` that first repairs Fortran's 3-digit-exponent tokens."""
    text = fix_fortran_exponents(Path(path).read_text())
    return np.loadtxt(io.StringIO(text), **kwargs)


def to_float(token: str) -> float:
    """``float()`` that tolerates a dropped-exponent token."""
    return float(fix_fortran_exponents(token))
