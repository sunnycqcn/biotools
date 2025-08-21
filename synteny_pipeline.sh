#!/bin/bash
# ============================
# Flexible Synteny Pipeline (jcvi)
# 支持任意数量基因组
# 默认跑所有 pairwise 组合；也可以用 -pairs 指定组合
#
# 使用示例（默认跑所有组合）:
#   ./synteny_pipeline.sh \
#       -gff1 genome1.gff -gff2 genome2.gff -gff3 genome3.gff \
#       -cds1 genome1.cds.fa -cds2 genome2.cds.fa -cds3 genome3.cds.fa \
#       -prefix G1 G2 G3
#
# 使用示例（只跑指定组合）:
#   ./synteny_pipeline.sh ... -prefix G1 G2 G3 -pairs "G1:G2,G1:G3"
# ============================

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
  case $1 in
    -gff*) gffs+=("$2"); shift 2 ;;
    -cds*) cdss+=("$2"); shift 2 ;;
    -prefix) prefixes=(); shift
             while [[ $# -gt 0 && $1 != -* ]]; do
               prefixes+=("$1"); shift
             done ;;
    -pairs) user_pairs="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# --- Check input ---
n=${#prefixes[@]}
if [[ ${#gffs[@]} -ne $n ]] || [[ ${#cdss[@]} -ne $n ]]; then
  echo "Error: prefix / gff / cds 数量必须相等!"
  exit 1
fi

# --- Activate conda env ---
echo "==> Activating conda env ..."
source /isilon/projects/J-002571_clubroot_src_chen/fuf/miniconda/bin/activate \
  /isilon/projects/J-002571_clubroot_src_chen/fuf/miniconda/envs/jcvi

# --- Step 1: GFF -> BED ---
echo "==> Converting GFF3 to BED ..."
for i in "${!gffs[@]}"; do
  prefix=${prefixes[$i]}
  gff=${gffs[$i]}
  python -m jcvi.formats.gff bed --type=mRNA "$gff" -o "${prefix}.bed"
done

# --- Step 2: CDS format ---
echo "==> Formatting CDS fasta ..."
for i in "${!cdss[@]}"; do
  prefix=${prefixes[$i]}
  cds=${cdss[$i]}
  python -m jcvi.formats.fasta format "$cds" "${prefix}.cds"
done

# --- Step 3: Generate pairs ---
pairs=()
if [[ -n "$user_pairs" ]]; then
  # 用户指定组合，例如 "DAR:ECD10,DAR:G3"
  IFS=',' read -r -a pair_list <<< "$user_pairs"
  for p in "${pair_list[@]}"; do
    pairs+=("$p")
  done
else
  # 默认：全 pairwise
  for ((i=0; i<n; i++)); do
    for ((j=i+1; j<n; j++)); do
      pairs+=("${prefixes[$i]}:${prefixes[$j]}")
    done
  done
fi

# --- Step 4: Run synteny for each pair ---
for pair in "${pairs[@]}"; do
  p1=${pair%%:*}
  p2=${pair##*:}
  echo "==> Running pairwise synteny: $p1 vs $p2"

  # Ortholog search
  python -m jcvi.compara.catalog ortholog "$p1" "$p2" --cscore=.99 --no_strip_names --notex

  # Dotplot
  python -m jcvi.graphics.dotplot "$p1.$p2.anchors" --notex

  # Depth check
  python -m jcvi.compara.synteny depth --histogram "$p1.$p2.anchors"

  # Simplified anchors
  python -m jcvi.compara.synteny screen --minspan=0 --simple \
    "$p1.$p2.anchors" "$p1.$p2.anchors.simple"
done

echo "==> Pipeline finished successfully!"
