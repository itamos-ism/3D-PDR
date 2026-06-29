"""Optional, fully-offline local-LLM backend for natural-language pathway summaries.

No internet is used at run time. Two backends are auto-detected, in order:

1. **Ollama** — a locally-running server (default ``http://localhost:11434``).
   Install Ollama once, ``ollama pull llama3.2`` (or any small instruct model),
   and it is detected automatically. Talked to over stdlib HTTP (no extra deps).
2. **llama.cpp** — a local GGUF model via the optional ``llama-cpp-python``
   package, pointed to by the ``PDRSTUDIO_LLM_GGUF`` environment variable.

If neither is present the feature degrades gracefully: :func:`available` returns
False and the page shows setup instructions while still rendering the
deterministic pathway table.

The model only *phrases* a digest that is computed deterministically from the
model output (see ``chemanalysis.pathway_evidence``); it is instructed to use
only the reactions/numbers provided, so it does not invent chemistry.
"""
from __future__ import annotations

import functools
import json
import os
import urllib.request

OLLAMA_HOST = os.environ.get("OLLAMA_HOST", "http://localhost:11434").rstrip("/")
GGUF_PATH = os.environ.get("PDRSTUDIO_LLM_GGUF", "")

# Small instruct models we prefer if several are installed (name prefixes).
_PREFERRED = ["llama3.2", "qwen2.5", "phi3", "phi", "gemma2", "mistral", "llama3", "llama"]

SYSTEM = (
    "You are an assistant for photodissociation-region (PDR) astrochemistry. You "
    "are given pre-computed dominant formation and destruction pathways of ONE "
    "species across a 1D cloud, ordered by depth (visual extinction Av; higher Av "
    "means deeper and more shielded, usually denser and colder). Write a concise, "
    "scientifically accurate summary (4-6 sentences) of how the species is FORMED "
    "and DESTROYED, and how the dominant pathways change from the cloud edge to "
    "the interior. Mention only the most important reactions; ignore minor "
    "contributors. Use ONLY the reactions and numbers provided — do not invent "
    "reactions, rate values, or species. Refer to reactions in arrow notation."
)


# ---------------------------------------------------------------------------
# backend detection
# ---------------------------------------------------------------------------
def _ollama_models() -> list[str]:
    try:
        with urllib.request.urlopen(f"{OLLAMA_HOST}/api/tags", timeout=1.5) as r:
            data = json.load(r)
        return [m["name"] for m in data.get("models", [])]
    except Exception:
        return []


def _pick(models: list[str]) -> str | None:
    for pref in _PREFERRED:
        for m in models:
            if m.split(":")[0].lower().startswith(pref):
                return m
    return models[0] if models else None


def _gguf_ready() -> bool:
    if not (GGUF_PATH and os.path.exists(GGUF_PATH)):
        return False
    try:
        import llama_cpp  # noqa: F401
        return True
    except Exception:
        return False


def available() -> bool:
    """True if any local LLM backend can be used right now."""
    return bool(_ollama_models()) or _gguf_ready()


def backend_name() -> str | None:
    models = _ollama_models()
    if models:
        return f"Ollama · {_pick(models)}"
    if _gguf_ready():
        return f"llama.cpp · {os.path.basename(GGUF_PATH)}"
    return None


def setup_hint() -> str:
    return (
        "No local LLM was detected, so the AI summary is disabled (the deterministic "
        "pathway table below still works). To enable fully-offline summaries, pick one:\n\n"
        "**Option A — Ollama (recommended).** Install **[Ollama](https://ollama.com)**, "
        "then in a terminal run:\n"
        "```bash\n"
        "ollama pull llama3.2      # or qwen2.5 / phi3 — any small instruct model\n"
        "```\n"
        "Ollama runs a local server on `http://localhost:11434` and PDR-studio picks it "
        "up automatically (no internet needed afterwards, no extra Python packages).\n\n"
        "**Option B — llama.cpp (a single GGUF file).** Install the binding into the "
        "app's `.venv` and point to a downloaded model:\n"
        "```bash\n"
        "pip install llama-cpp-python\n"
        "export PDRSTUDIO_LLM_GGUF=/path/to/model.gguf\n"
        "```\n"
        "Then restart PDR-studio. Both options run entirely on your machine."
    )


# ---------------------------------------------------------------------------
# generation
# ---------------------------------------------------------------------------
def _ollama_generate(model: str, prompt: str, system: str, max_tokens: int) -> str:
    body = json.dumps({
        "model": model, "prompt": prompt, "system": system, "stream": False,
        "options": {"temperature": 0.2, "num_predict": max_tokens},
    }).encode()
    req = urllib.request.Request(f"{OLLAMA_HOST}/api/generate", data=body,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r).get("response", "").strip()


@functools.lru_cache(maxsize=1)
def _llama(path: str):
    from llama_cpp import Llama
    return Llama(model_path=path, n_ctx=4096, verbose=False)


def _llamacpp_generate(prompt: str, system: str, max_tokens: int) -> str:
    llm = _llama(GGUF_PATH)
    msgs = ([{"role": "system", "content": system}] if system else []) \
        + [{"role": "user", "content": prompt}]
    out = llm.create_chat_completion(messages=msgs, temperature=0.2, max_tokens=max_tokens)
    return out["choices"][0]["message"]["content"].strip()


def generate(prompt: str, system: str = SYSTEM, max_tokens: int = 400) -> str:
    """Run the prompt on the first available local backend."""
    models = _ollama_models()
    if models:
        return _ollama_generate(_pick(models), prompt, system, max_tokens)
    if _gguf_ready():
        return _llamacpp_generate(prompt, system, max_tokens)
    raise RuntimeError("No local LLM backend available.")


# ---------------------------------------------------------------------------
# prompt assembly + high-level entry point
# ---------------------------------------------------------------------------
def build_prompt(ev: dict) -> str:
    """Render a pathway-evidence dict (chemanalysis.pathway_evidence) into the
    compact, factual context the model summarises."""
    sp = ev["species"]
    dens = (f"{ev['nH_range'][0]:.2e} cm-3 (uniform)" if ev["uniform_density"]
            else f"{ev['nH_range'][0]:.2e}–{ev['nH_range'][1]:.2e} cm-3")
    ab = ev["abundance"]
    L = [
        f"Species: {sp}",
        f"Cloud: {ev['npoints']} depth points; Av {ev['av_range'][0]:.3g}–"
        f"{ev['av_range'][1]:.3g} mag; n_H {dens}; "
        f"Tgas {ev['Tgas_range'][1]:.3g}→{ev['Tgas_range'][0]:.3g} K (edge→interior).",
        f"Abundance n({sp})/n_H: {ab['edge']:.2e} at the edge, peak {ab['max']:.2e} "
        f"near Av={ab['peak_Av']:.3g}, {ab['deep']:.2e} in the interior.",
        "",
        "Dominant FORMATION pathways across the column:",
    ]
    for r in ev["formation_top"]:
        L.append(f"  - {r['reaction']}  (~{r['typical_pct']}% of formation where it "
                 f"dominates; top route over ~{r['dominant_fraction'] * 100:.0f}% of the "
                 f"column, Av {r['av_lo']:.3g}–{r['av_hi']:.3g})")
    L.append("Dominant DESTRUCTION pathways across the column:")
    for r in ev["destruction_top"]:
        L.append(f"  - {r['reaction']}  (~{r['typical_pct']}% of destruction where it "
                 f"dominates; top route over ~{r['dominant_fraction'] * 100:.0f}% of the "
                 f"column, Av {r['av_lo']:.3g}–{r['av_hi']:.3g})")
    L += ["", "Snapshots at representative depths:"]
    for s in ev["snapshots"]:
        f = "; ".join(f"{x['reaction']} ({x['pct']}%)" for x in s["formation"]) or "—"
        d = "; ".join(f"{x['reaction']} ({x['pct']}%)" for x in s["destruction"]) or "—"
        L.append(f"  [{s['label']}] Av={s['Av']:.3g}, n_H={s['nH']:.2e}, "
                 f"Tgas={s['Tgas']:.3g} K, n({sp})/n_H={s['abundance']:.2e}")
        L.append(f"      formation:   {f}")
        L.append(f"      destruction: {d}")
    return "\n".join(L)


def summarize(evidence: dict) -> str:
    """Build the prompt from evidence and run it on the local backend."""
    if not evidence:
        return ""
    return generate(build_prompt(evidence))
