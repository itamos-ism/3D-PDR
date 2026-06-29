#!/usr/bin/env bash
# Launch PDR-studio so others on your intranet can reach it at http://<your-ip>:8501
#
# SECURITY: anyone who can reach that address gets the full app — which runs/compiles
# 3D-PDR and executes the Python you type in the "Custom plot" tabs ON THIS MACHINE.
# Only use this on a trusted internal network.
set -euo pipefail
cd "$(dirname "$0")"

PORT="${PORT:-8501}"
IP="$(hostname -I 2>/dev/null | awk '{print $1}')"

# Optional: enable offline AI summaries (see run.sh).
# export PDRSTUDIO_LLM_GGUF="$HOME/models/Llama-3.2-3B-Instruct-Q4_K_M.gguf"

if [ ! -d .venv ]; then
    echo "Creating virtual environment (.venv) ..."
    python3 -m venv .venv
    .venv/bin/pip install --upgrade pip
    .venv/bin/pip install -r requirements.txt
fi

echo "============================================================"
echo " Share this with colleagues on the intranet:"
echo "     http://${IP:-<your-ip>}:${PORT}"
echo " (Ctrl-C to stop. If they can't connect, open the firewall:"
echo "     sudo ufw allow ${PORT}/tcp )"
echo "============================================================"

exec .venv/bin/streamlit run app.py \
    --server.address=0.0.0.0 \
    --server.port="${PORT}" \
    --server.headless=true \
    --browser.gatherUsageStats=false "$@"
