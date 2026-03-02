#!/usr/bin/env bash
# Patch myst-to-tex's kindToEnvironment to support algorithm and property
# proof environments in LaTeX/PDF export.
#
# Run after installing dependencies:
#   pip install -r requirements.txt && bash scripts/patch-mystmd.sh
#
# This adds 'algorithm' and 'property' to the proof environment switch
# statement in myst-to-tex, which is bundled inside jupyter-book and mystmd.
#
# Needs to be re-run after upgrading jupyter-book or mystmd.

set -euo pipefail

patch_file() {
    local CJS="$1"
    local NAME="$2"

    if [ ! -f "$CJS" ]; then
        echo "  SKIP: $NAME not found at $CJS"
        return 0
    fi

    if grep -q 'case "algorithm":' "$CJS" 2>/dev/null; then
        echo "  OK: $NAME already patched"
        return 0
    fi

    # Insert algorithm and property cases after 'return "corollary";'
    sed -i '/return "corollary";/a\    case "algorithm":\n      return "algorithm";\n    case "property":\n      return "property";' "$CJS"

    if grep -q 'case "algorithm":' "$CJS" 2>/dev/null; then
        echo "  PATCHED: $NAME — algorithm and property environments enabled"
    else
        echo "  WARNING: $NAME patch may not have applied correctly"
        return 1
    fi
}

echo "Patching myst-to-tex for algorithm/property support..."

# Patch jupyter-book's bundled copy
JB_CJS=$(python3 -c "
import jupyter_book, os
print(os.path.join(os.path.dirname(jupyter_book.__file__), 'dist', 'jupyter-book.cjs'))
" 2>/dev/null) || true
if [ -n "$JB_CJS" ]; then
    patch_file "$JB_CJS" "jupyter-book"
fi

# Patch mystmd's bundled copy (used by 'myst build' directly)
MYST_CJS=$(python3 -c "
import mystmd_py, os
print(os.path.join(os.path.dirname(mystmd_py.__file__), 'myst.cjs'))
" 2>/dev/null) || true
if [ -n "$MYST_CJS" ]; then
    patch_file "$MYST_CJS" "mystmd"
fi

echo "Done."
