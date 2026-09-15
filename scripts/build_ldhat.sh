#!/usr/bin/env bash
# Build LDhat for `ld_recombination --emit ldhat`.
#
# LDhat (McVean lab; https://github.com/auton1/LDhat) is the recommended recombination-rate
# estimator for Plasmodium: pyrho's fused-LASSO over-smooths Pf SNP density and collapses
# each chromosome to ~one rate, while LDhat's `interval` rjMCMC resolves gene-scale
# structure at the same sample size. There is no bioconda build for osx-arm64, and two
# things bite when building from source, both handled here:
#
#   1. LDhat is old C. Modern clang defaults to C23, which turns its K&R-style
#      non-prototype calls into hard ERRORS. We build with `-std=gnu89` and silence the
#      two relevant warnings.
#   2. On Apple Silicon a NATIVE arm64 LDhat binary bus-errors at run time (misaligned
#      access in the old pointer code). The fix is to build for x86_64 and run under
#      Rosetta. This script does that automatically on Darwin/arm64 -- and builds ALL of
#      interval/lkgen/stat the same way, because `ld_recombination` runs them under one
#      `arch -x86_64` prefix and a mixed-architecture install breaks that.
#
# Usage:
#   scripts/build_ldhat.sh                 # clone + build + install to ~/.local/share/ldhat
#   scripts/build_ldhat.sh --prefix DIR    # install the binaries + lookup table into DIR
#   scripts/build_ldhat.sh --src DIR       # build an existing LDhat checkout instead of cloning
#   scripts/build_ldhat.sh --ref TAG       # git ref to clone (default: master)
#
# Then point the subcommand at it (the script prints this at the end):
#   export LDHAT_DIR=<prefix>
#   plasgenomicsutils ld_recombination ... --emit ldhat            # uses $LDHAT_DIR
#   plasgenomicsutils ld_recombination ... --emit ldhat --ldhat-dir <prefix>
set -uo pipefail

PREFIX="${LDHAT_DIR:-$HOME/.local/share/ldhat}"
SRC=""
REPO="https://github.com/auton1/LDhat.git"
REF="master"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX=$2; shift 2 ;;
    --src)    SRC=$2; shift 2 ;;
    --repo)   REPO=$2; shift 2 ;;
    --ref)    REF=$2; shift 2 ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

die() { printf '\033[31mbuild_ldhat: %s\033[0m\n' "$*" >&2; exit 1; }
say() { printf '\033[1m== %s ==\033[0m\n' "$*"; }

# --- toolchain: native everywhere except Apple Silicon, where we cross to x86_64 --------
CFLAGS="-O2 -std=gnu89 -Wno-implicit-function-declaration -Wno-deprecated-non-prototype"
ARCH_PREFIX=()
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "arm64" ]]; then
  arch -x86_64 true 2>/dev/null \
    || die "Rosetta is required to run x86_64 LDhat on Apple Silicon. Install it with:
           softwareupdate --install-rosetta --agree-to-license"
  CC="clang -arch x86_64"
  ARCH_PREFIX=(arch -x86_64)
  say "Apple Silicon: building LDhat for x86_64 (runs under Rosetta)"
else
  CC="${CC:-cc}"
  say "building LDhat natively with: $CC"
fi

# --- get the source ---------------------------------------------------------------------
CLEAN_SRC=0
if [[ -z "$SRC" ]]; then
  command -v git >/dev/null || die "git not found (needed to clone $REPO; or pass --src)"
  SRC=$(mktemp -d)/LDhat
  CLEAN_SRC=1
  say "cloning $REPO@$REF"
  git clone --depth 1 --branch "$REF" "$REPO" "$SRC" 2>/dev/null \
    || git clone --depth 1 "$REPO" "$SRC" || die "git clone failed"
fi
[[ -f "$SRC/makefile" || -f "$SRC/Makefile" ]] || die "no makefile in $SRC -- not an LDhat checkout"
[[ -f "$SRC/interval.c" || -f "$SRC/pair_int.c" ]] || die "$SRC does not look like LDhat source"

# --- build interval + lkgen + stat ------------------------------------------------------
say "compiling interval, lkgen, stat"
( cd "$SRC" && make clean >/dev/null 2>&1
  "${ARCH_PREFIX[@]}" make CC="$CC" CFLAGS="$CFLAGS" interval lkgen stat ) \
  || die "make failed (see output above)"
for b in interval lkgen stat; do
  [[ -x "$SRC/$b" ]] || die "$b was not built"
done

# --- install ----------------------------------------------------------------------------
say "installing to $PREFIX"
mkdir -p "$PREFIX/lk_files" || die "cannot create $PREFIX"
cp -f "$SRC/interval" "$SRC/lkgen" "$SRC/stat" "$PREFIX/"
# the bundled precomputed lookup table (n=120, theta=0.001) that --ldhat-lk defaults to
if [[ -f "$SRC/lk_files/lk_n120_t0.001.gz" ]]; then
  cp -f "$SRC/lk_files/lk_n120_t0.001.gz" "$PREFIX/lk_files/"
else
  printf '\033[33mnote: no lk_files/lk_n120_t0.001.gz in the source; pass --ldhat-lk pointing\n      at a lookup table whose sequence count >= your panel size.\033[0m\n'
fi

# --- verify: uniform architecture + lkgen actually runs ---------------------------------
say "verifying"
archs=$(for b in interval lkgen stat; do file "$PREFIX/$b" | grep -oE 'x86_64|arm64'; done | sort -u)
echo "  binary architecture(s): $archs"
[[ $(wc -l <<<"$archs") -eq 1 ]] || die "binaries are NOT a single architecture -- the arch prefix will break one of them"
if [[ -f "$PREFIX/lk_files/lk_n120_t0.001.gz" ]]; then
  tmp=$(mktemp -d); ( cd "$tmp" && gunzip -c "$PREFIX/lk_files/lk_n120_t0.001.gz" > lk.txt \
    && "${ARCH_PREFIX[@]}" "$PREFIX/lkgen" -lk lk.txt -nseq 20 >/dev/null 2>&1 \
    && head -1 new_lk.txt | grep -q '^20 ' ) \
    && echo "  lkgen smoke test: OK" || die "lkgen smoke test failed -- the binaries may not run on this machine"
  rm -rf "$tmp"
fi
[[ $CLEAN_SRC -eq 1 ]] && rm -rf "$(dirname "$SRC")"

cat <<EOF

$(printf '\033[32mLDhat installed to %s\033[0m' "$PREFIX")

Point the subcommand at it:
  export LDHAT_DIR="$PREFIX"          # add to your shell profile to make it permanent
  plasgenomicsutils ld_recombination --vcf CALLSET.bcf --fws-table fws.tsv \\
      --region Pf3D7_07_v3 --maf 0.02 --emit ldhat --out-prefix recomb_map

Or per-run without the env var:  ... --emit ldhat --ldhat-dir "$PREFIX"
EOF
