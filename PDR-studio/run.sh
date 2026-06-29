#!/usr/bin/env bash
# Launch PDR-studio. Creates the venv on first run if missing.
export PDRSTUDIO_LLM_GGUF="$HOME/data/models/Llama-3.2-3B-Instruct-Q4_K_M.gguf"
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d .venv ]; then
    echo "Creating virtual environment (.venv) ..."
    python3 -m venv .venv
    .venv/bin/pip install --upgrade pip
    .venv/bin/pip install -r requirements.txt
fi

exec .venv/bin/streamlit run app.py "$@"
