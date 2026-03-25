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

echo "Done with algorithm/property patch."

# ── Patch 2: Reduce HTML static build concurrency to 1 ────────────────
# The default p-limit(5) causes 5 pages to render simultaneously,
# which can OOM on GitHub runners (16GB). We reduce to 1 and add
# a for-loop instead of Promise.all to allow GC between pages.

patch_html_concurrency() {
    local CJS="$1"
    local NAME="$2"

    if [ ! -f "$CJS" ]; then
        echo "  SKIP: $NAME not found at $CJS"
        return 0
    fi

    if grep -q 'sequential HTML rendering' "$CJS" 2>/dev/null; then
        echo "  OK: $NAME HTML concurrency already patched"
        return 0
    fi

    # Replace Promise.all(routes.map(...)) with sequential for-loop
    python3 -c "
import re

with open('$CJS', 'r') as f:
    content = f.read()

old = 'await Promise.all(routes.map(async (route) => limitConnections(async () => {'
new = '''// Patched: sequential HTML rendering to reduce memory pressure
  for (const route of routes) {'''

content = content.replace(old, new, 1)

# Close the limitConnections wrapper - replace the matching '})));' with '}'
# Find the specific pattern after the route processing block
content = content.replace(
    '}\\n  })));\\n  appServer.stop();',
    '}\\n  }\\n  appServer.stop();',
    1
)

# If the above didn't match, try a more targeted approach
if 'sequential HTML rendering' in content and '})));' in content:
    # Replace the first }))) after the sequential comment
    idx = content.index('sequential HTML rendering')
    rest = content[idx:]
    rest = rest.replace('})));', '}', 1)
    content = content[:idx] + rest

with open('$CJS', 'w') as f:
    f.write(content)
" 2>/dev/null

    if grep -q 'sequential HTML rendering' "$CJS" 2>/dev/null; then
        echo "  PATCHED: $NAME — sequential HTML rendering enabled"
    else
        echo "  WARNING: $NAME HTML concurrency patch may not have applied"
    fi
}

echo ""
echo "Patching HTML build concurrency..."

if [ -n "$MYST_CJS" ]; then
    patch_html_concurrency "$MYST_CJS" "mystmd"
fi

echo "Done with algorithm/property patch."

# ── Patch 2: Sequential HTML rendering to reduce memory ───────────────
# The default parallel fetch (p-limit 5) causes Remix to accumulate
# memory for all pages, OOMing on 16GB GitHub runners.
echo ""
echo "Patching HTML build for sequential rendering..."

if [ -n "$MYST_CJS" ]; then
    python3 scripts/patch-html-sequential.py "$MYST_CJS"
fi

echo "Done."
