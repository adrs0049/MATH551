#!/usr/bin/env bash
#
# Execute all notebooks in-place so their outputs are cached.
# Usage:
#   ./scripts/execute-notebooks.sh              # execute all
#   ./scripts/execute-notebooks.sh notebooks/stiffness.ipynb  # execute one
#
set -uo pipefail

TIMEOUT=${NOTEBOOK_TIMEOUT:-300}

if [ $# -gt 0 ]; then
    notebooks=("$@")
else
    mapfile -t notebooks < <(find notebooks -name '*.ipynb' -not -path '*/.ipynb_checkpoints/*' | sort)
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
    exit 1
fi
