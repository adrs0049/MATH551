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

# ── Patch 2: Sequential HTML rendering (DISABLED) ─────────────────────
# Disabled with mystmd 1.9.0 / jupyter-book 2.1.5 upgrade. The original
# variant restarted the Remix server every 15 pages to free memory, but
# Remix v1.17 takes ~30s to rebuild on restart — longer than the fetch
# retry window — so post-restart fetches fail and the post-fetch
# copySync step never runs, producing artifacts with no
# _build/html/build/_assets/. The sequential-without-restart variant
# also hits the ENOENT race more often than mystmd's default concurrent
# rendering (longer wall time = more chances for Remix's dev cleanup to
# delete public/build files before copySync reads them). The workflow's
# snapshot + self-heal step handles the ENOENT race correctly, so we
# rely on that and let mystmd run its default Promise.all + p-limit(5).
# Re-enable if 73-page builds OOM on the CI runner.

echo ""
echo "All patches applied."
