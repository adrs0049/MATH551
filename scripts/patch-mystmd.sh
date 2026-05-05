#!/usr/bin/env bash
# Patch mystmd for:
#   1. Algorithm/property proof environments in LaTeX/PDF export
#   2. Sequential HTML rendering to reduce memory on CI
#
# Run after installing dependencies:
#   pip install -r requirements.txt && bash scripts/patch-mystmd.sh

set -euo pipefail

# ── Locate the bundled CJS files ──────────────────────────────────────

JB_CJS=$(python3 -c "
import jupyter_book, os
print(os.path.join(os.path.dirname(jupyter_book.__file__), 'dist', 'jupyter-book.cjs'))
" 2>/dev/null) || true

MYST_CJS=$(python3 -c "
import mystmd_py, os
print(os.path.join(os.path.dirname(mystmd_py.__file__), 'myst.cjs'))
" 2>/dev/null) || true

# ── Patch 1: Algorithm/property environments ──────────────────────────

patch_proof_envs() {
    local CJS="$1"
    local NAME="$2"

    if [ ! -f "$CJS" ]; then
        echo "  SKIP: $NAME not found"
        return 0
    fi

    if grep -q 'case "algorithm":' "$CJS" 2>/dev/null; then
        echo "  OK: $NAME already patched"
        return 0
    fi

    sed -i '/return "corollary";/a\    case "algorithm":\n      return "algorithm";\n    case "property":\n      return "property";' "$CJS"

    if grep -q 'case "algorithm":' "$CJS" 2>/dev/null; then
        echo "  PATCHED: $NAME"
    else
        echo "  WARNING: $NAME patch may not have applied"
    fi
}

echo "Patch 1: Algorithm/property proof environments"
[ -n "$JB_CJS" ] && patch_proof_envs "$JB_CJS" "jupyter-book"
[ -n "$MYST_CJS" ] && patch_proof_envs "$MYST_CJS" "mystmd"

# ── Patch 2: SIGKILL the Remix dev server on appServer.stop() ────────
# Stops Remix's exit handler from deleting public/build/ before mystmd's
# final copySync reads it. Without this, copySync ENOENTs and the artifact
# either fails or — worse — the snapshot self-heal restores prod-build
# hashes that don't match the dev-rendered HTML, so every asset 404s.

echo ""
echo "Patch 2: SIGKILL kill signal for dev server"
patched=false
if [ -n "$MYST_CJS" ]; then
    python3 scripts/patch-kill-signal.py "$MYST_CJS" && patched=true
fi
if [ -n "$JB_CJS" ]; then
    python3 scripts/patch-kill-signal.py "$JB_CJS" && patched=true
fi
if [ "$patched" = false ]; then
    echo "  SKIP: no CJS files found to patch"
fi

# ── Patch 3: Sequential HTML rendering ────────────────────────────────
# mystmd's default Promise.all + p-limit(5) holds 5 routes' worth of
# Remix-rendered React trees in memory at once and OOMs the 16GB CI
# runner around page 50. Sequential rendering keeps peak memory low
# enough to fit. (Earlier variant also restarted the Remix server every
# 15 pages, but Remix v1.17 takes ~30s to rebuild on restart — longer
# than the fetch retry window — so post-restart fetches fail. The
# no-restart variant here lets the server run continuously; combined
# with the SIGKILL patch above it produces a complete artifact.)

echo ""
echo "Patch 3: Sequential HTML rendering"
patched=false
if [ -n "$MYST_CJS" ]; then
    python3 scripts/patch-html-sequential.py "$MYST_CJS" && patched=true
fi
if [ -n "$JB_CJS" ]; then
    python3 scripts/patch-html-sequential.py "$JB_CJS" && patched=true
fi
if [ "$patched" = false ]; then
    echo "  SKIP: no CJS files found to patch"
fi

echo ""
echo "All patches applied."
