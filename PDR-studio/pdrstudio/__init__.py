"""PDR-studio backend package.

Thin, well-separated helpers that the Streamlit UI calls. Each module owns one
concern (paths, config.mk, params.dat, chemistry files, compiling/running,
loading outputs, debug reports) so new features can be added without touching
the rest.
"""

__version__ = "0.1.0"
