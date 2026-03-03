#!/usr/bin/env bash
#
# Execute notebooks referenced in myst.yml (the TOC) so their outputs are cached.
# Usage:
#   ./scripts/execute-notebooks.sh              # execute TOC notebooks
#   ./scripts/execute-notebooks.sh --all        # execute all notebooks
#   ./scripts/execute-notebooks.sh notebooks/stiffness.ipynb  # execute one
#
set -uo pipefail

TIMEOUT=${NOTEBOOK_TIMEOUT:-300}
MYSTyml="myst.yml"

if [ $# -gt 0 ] && [ "$1" = "--all" ]; then
    # All notebooks in the notebooks/ directory
    mapfile -t notebooks < <(find notebooks -name '*.ipynb' -not -path '*/.ipynb_checkpoints/*' | sort)
elif [ $# -gt 0 ]; then
    # Explicit list
    notebooks=("$@")
else
    # Only notebooks referenced in myst.yml
    mapfile -t notebooks < <(grep '\.ipynb' "$MYSTyml" | sed 's/.*file: *//' | sed 's/ *$//' | sort)
fi

total=${#notebooks[@]}
passed=0
failed=0
failed_list=()

echo "Executing $total notebooks (timeout=${TIMEOUT}s per notebook)"
echo "=============================================="

for nb in "${notebooks[@]}"; do
    echo -n "  $nb ... "

    if jupyter execute "$nb" \
        --inplace \
        --timeout="$TIMEOUT" \
        2>/dev/null; then
        echo "OK"
        passed=$((passed + 1))
    else
        echo "FAILED"
        failed=$((failed + 1))
        failed_list+=("$nb")
    fi
done

echo "=============================================="
echo "Results: $passed passed, $failed failed (of $total)"

if [ ${#failed_list[@]} -gt 0 ]; then
    echo ""
    echo "Failed notebooks:"
    for nb in "${failed_list[@]}"; do
        echo "  - $nb"
    done
    echo ""
    echo "WARNING: Some notebooks failed but continuing (outputs may be stale)."
fi
