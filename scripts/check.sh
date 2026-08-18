#!/usr/bin/env bash
# Run what CI runs, before pushing.
#
# CI for this repo is two jobs (.github/workflows/tests.yml): `pytest -q` across a
# Python matrix, and `mkdocs build --strict`. The second is the one that catches a docs
# page nobody linked, and it is quick -- so `--fast` runs the checks that take seconds
# and skips the full suite, for the loop between edits.
#
#   scripts/check.sh          # everything CI runs
#   scripts/check.sh --fast   # docs + contract tests only (a few seconds)
set -uo pipefail
cd "$(dirname "$0")/.."

FAST=0
[[ "${1:-}" == "--fast" ]] && FAST=1
FAILED=()

run() {                       # run <label> <command...>
  local label=$1; shift
  printf '\n\033[1m== %s ==\033[0m\n' "$label"
  if "$@"; then
    printf '\033[32mOK\033[0m  %s\n' "$label"
  else
    printf '\033[31mFAILED\033[0m  %s\n' "$label"
    FAILED+=("$label")
  fi
}

# The docs contract: every command in docs/commands.md, every documented flag real, every
# docs page in the mkdocs nav. Seconds, and it catches most of what the docs job would.
run "docs contract tests" python -m pytest tests/test_docs_contract.py -q

if command -v mkdocs >/dev/null 2>&1; then
  run "mkdocs build --strict" mkdocs build --strict
  rm -rf site
else
  printf '\n\033[33mSKIPPED\033[0m  mkdocs build --strict (pip install -r docs/requirements.txt)\n'
fi

if [[ $FAST -eq 0 ]]; then
  run "pytest" python -m pytest -q
fi

printf '\n'
if [[ ${#FAILED[@]} -eq 0 ]]; then
  printf '\033[32mall checks passed\033[0m\n'
else
  printf '\033[31m%d check(s) failed:\033[0m %s\n' "${#FAILED[@]}" "${FAILED[*]}"
  exit 1
fi
