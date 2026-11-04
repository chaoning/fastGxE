#!/usr/bin/env bash
# Checks that the three-step GxE workflow agrees with the single-command mode, and
# that the results do not depend on how the genotypes are stored or the work divided.
#
#   test/test_steps.sh <fastgxe2.0 binary> <new work directory>
#
# Needs plink 1.9, found on PATH or through $PLINK. Set THREADS to change the number
# of threads (default 4).
set -uo pipefail

FASTGXE=$(readlink -f "$1")
WORK=$2
REPO=$(cd "$(dirname "$0")/.." && pwd)
EX=$REPO/example
PLINK=${PLINK:-plink}
THREADS=${THREADS:-4}
SEED=2026
failures=0

mkdir "$WORK" || exit 1
WORK=$(readlink -f "$WORK")
cd "$WORK" || exit 1

cat > compare_tables.py <<'PY'
import sys

def agree(a, b):
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    return abs(x - y) <= 1e-5 * max(abs(x), abs(y))

with open(sys.argv[1]) as f:
    first = f.read().splitlines()
with open(sys.argv[2]) as f:
    second = f.read().splitlines()
if len(first) != len(second):
    sys.exit(1)
for a, b in zip(first, second):
    fields_a, fields_b = a.split('\t'), b.split('\t')
    if len(fields_a) != len(fields_b) or not all(map(agree, fields_a, fields_b)):
        sys.exit(1)
PY

run() { "$FASTGXE" "$@" --threads "$THREADS" >> run.log 2>&1 || { echo "command failed: fastgxe2.0 $*"; exit 1; }; }
same() {
    if cmp -s "$1" "$2"; then echo "ok    $3"; else echo "FAIL  $3"; failures=$((failures + 1)); fi
}
rejected() {
    local label=$1; shift
    if "$FASTGXE" "$@" >> run.log 2>&1; then echo "FAIL  $label"; failures=$((failures + 1)); else echo "ok    $label"; fi
}
# Numbers must agree to a relative 1e-5. Results are only byte-identical for the same
# SNP blocks: other blocks change the last bits of the score, and near its null mean
# (p around 0.46) saddle() turns that into differences of about 1e-6 in p_multi.
close() {
    if python3 "$WORK/compare_tables.py" "$1" "$2"; then
        echo "ok    $3"
    else
        echo "FAIL  $3"; failures=$((failures + 1))
    fi
}
merge() { awk 'FNR>1 || NR==1' "$@"; }
no_order() { cut -f2- "$1"; }

echo "== preparing genotypes and GRM"
# Keep the file order of the chromosomes, so the split genome lists SNPs in the same order.
chromosomes=$(awk '{print $1}' "$EX/test.bim" | uniq)
for c in $chromosomes; do
    "$PLINK" --bfile "$EX/test" --chr "$c" --keep-allele-order --make-bed --out "chr$c" >> run.log 2>&1
    echo "chr$c" >> chr_list.txt
done
run --make-grm --bfile "$EX/test" --out grm
run --process-grm --group --grm grm --cut-value 0.05

MODEL="--grm grm --data $EX/test_simu_pheno.txt --trait pheno --env-int E1:E40"

echo "== three steps against the single command"
run --test-gxe $MODEL --bfile "$EX/test" --seed $SEED --out one
run --step 1 $MODEL --out s1
run --step 2 --null s1 --bfile "$EX/test" --seed $SEED --out s2
run --step 3 --null s1 --gamma s2 --bfile "$EX/test" --out s3
same one.var s1.var ".var"
same one.h2 s1.h2 ".h2"
same one.gamma s2.gamma ".gamma"
same one.assoc s3.assoc ".assoc"
same one.perE.assoc s3.perE.assoc ".perE.assoc"

echo "== gamma from split genotypes"
run --step 2 --null s1 --mbfile chr_list.txt --seed $SEED --out s2_split
same s2.gamma s2_split.gamma "--mbfile draws the same SNPs as --bfile"

echo "== random panel is topped up when SNPs fail QC"
run --step 2 --null s1 --bfile "$EX/test" --maf 0.3 --seed $SEED --keep-random --out s2_strict
run --step 2 --null s1 --mbfile chr_list.txt --maf 0.3 --seed $SEED --out s2_strict_split
same s2_strict.gamma s2_strict_split.gamma "--mbfile tops up with the same SNPs as --bfile"
drawn=$(grep 'Random panel' run.log | tail -1 | sed 's/.*used, \([0-9]*\) drawn/\1/')
if [ "$drawn" -gt 1000 ]; then echo "ok    QC failures made it draw more ($drawn SNPs)"; else echo "FAIL  drew only $drawn SNPs, so nothing was topped up"; failures=$((failures + 1)); fi
kept=$(awk -F'\t' 'NR>1 && $2>0.3 && $2<0.7 {n++; seen[$1]=1} END {print n+0, length(seen)}' s2_strict.main.random.assoc)
if [ "$kept" = "1000 1000" ]; then echo "ok    panel holds 1000 distinct SNPs that all pass QC"; else echo "FAIL  panel: $kept"; failures=$((failures + 1)); fi

echo "== scan by chromosome"
files=(); perE_files=()
while read -r fileset; do
    run --step 3 --null s1 --gamma s2 --bfile "$fileset" --out "s3.$fileset"
    files+=("s3.$fileset.assoc"); perE_files+=("s3.$fileset.perE.assoc")
done < chr_list.txt
merge "${files[@]}" > by_chr.assoc
merge "${perE_files[@]}" > by_chr.perE.assoc
close <(no_order s3.assoc) <(no_order by_chr.assoc) ".assoc by chromosome (without order)"
close <(no_order s3.perE.assoc) <(no_order by_chr.perE.assoc) ".perE.assoc by chromosome (without order)"

echo "== --split-task"
for i in 1 2 3 4 5 6 7; do
    run --step 3 --null s1 --gamma s2 --bfile "$EX/test" --split-task 7 $i --out s3_split
done
merge s3_split.7_{1..7}.assoc > split.assoc
merge s3_split.7_{1..7}.perE.assoc > split.perE.assoc
close s3.assoc split.assoc ".assoc in 7 parts"
close s3.perE.assoc split.perE.assoc ".perE.assoc in 7 parts"

for i in 1 2 3 4 5; do
    run --step 3 --null s1 --gamma s2 --mbfile chr_list.txt --split-task 5 $i --out s3_msplit
done
merge s3_msplit.5_{1..5}.assoc > msplit.assoc
close s3.assoc msplit.assoc ".assoc with --mbfile, in 5 parts across chromosomes"

for i in 1 2 3; do
    run --step 3 --null s1 --gamma s2 --bfile chr1 --split-task 3 $i --out s3_chr1
done
merge s3_chr1.3_{1..3}.assoc > chr1_split.assoc
close s3.chr1.assoc chr1_split.assoc ".assoc of one chromosome in 3 parts"

echo "== --exact-cut 0 needs no V^-1"
run --test-gxe $MODEL --bfile "$EX/test" --seed $SEED --exact-cut 0 --out one_noexact
run --step 3 --null s1 --gamma s2 --bfile "$EX/test" --exact-cut 0 --out s3_noexact
same one_noexact.assoc s3_noexact.assoc ".assoc with --exact-cut 0"

# Cut the null-model file just before V^-1 (see fastgxe/null_model_io.cpp for the layout).
vinv_offset=$(python3 - <<'PY'
import struct
with open('s1.null.bin', 'rb') as f:
    header = f.read(48)
    num_id, num_env, num_varcom = struct.unpack('<qqq', header[24:48])
    for _ in range(1 + num_id + num_env):
        size = struct.unpack('<q', f.read(8))[0]
        f.seek(size, 1)
    print(f.tell() + 8 * (num_varcom + num_id + num_id * num_env))
PY
)
head -c "$vinv_offset" s1.null.bin > cut.null.bin
run --step 3 --null s1 --gamma s2 --bfile chr22 --exact-cut 0 --out s3_full
run --step 3 --null cut --gamma s2 --bfile chr22 --exact-cut 0 --out s3_cut
same s3_full.assoc s3_cut.assoc "--exact-cut 0 does not read V^-1"
rejected "a null model without V^-1 is rejected when V^-1 is needed" --step 2 --null cut --bfile chr22 --out x

echo "== .fam order differs between filesets (simulated samples)"
python3 - <<'PY'
import numpy as np
rng = np.random.default_rng(1)
n, m, q = 300, 4000, 3
ids = [f"S{i:04d}" for i in range(n)]
dosage = rng.binomial(2, rng.uniform(0.1, 0.9, m), (n, m))
code = np.array([0b11, 0b10, 0b00])  # .bed codes of dosage 0, 1 and 2

def write_fileset(prefix, samples, snps, chrom):
    with open(prefix + ".fam", "w") as f:
        for s in samples:
            f.write(f"{ids[s]} {ids[s]} 0 0 0 -9\n")
    with open(prefix + ".bim", "w") as f:
        for j in snps:
            f.write(f"{chrom} snp{j} 0 {j + 1} A G\n")
    bits = code[dosage[np.ix_(samples, snps)]].T  # SNPs x samples
    padded = np.zeros((bits.shape[0], (len(samples) + 3) // 4 * 4), dtype=np.uint8)
    padded[:, :len(samples)] = bits
    packed = padded[:, ::4] | padded[:, 1::4] << 2 | padded[:, 2::4] << 4 | padded[:, 3::4] << 6
    with open(prefix + ".bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01]) + packed.tobytes())

everyone, first, second = np.arange(n), np.arange(m // 2), np.arange(m // 2, m)
write_fileset("syn", everyone, np.arange(m), 1)
write_fileset("syn_a", everyone, first, 1)
write_fileset("syn_b", everyone, second, 2)
write_fileset("syn_b_reordered", rng.permutation(n), second, 2)
write_fileset("syn_b_incomplete", everyone[:-1], second, 2)

env = rng.standard_normal((n, q))
beta = rng.normal(0, 0.02, m)
y = (dosage - dosage.mean(axis=0)) @ beta + env @ np.array([0.3, -0.2, 0.1]) + rng.standard_normal(n)
with open("syn_pheno.txt", "w") as f:
    f.write("iid y " + " ".join(f"E{j + 1}" for j in range(q)) + "\n")
    for i in range(n):
        f.write(f"{ids[i]} {y[i]:.6f} " + " ".join(f"{v:.6f}" for v in env[i]) + "\n")
PY
printf 'syn_a\nsyn_b\n' > syn_list.txt
printf 'syn_a\nsyn_b_reordered\n' > syn_reordered_list.txt
printf 'syn_a\nsyn_b_incomplete\n' > syn_incomplete_list.txt
run --make-grm --bfile syn --out syn_grm
run --process-grm --group --grm syn_grm --cut-value 0.05
run --step 1 --grm syn_grm --data syn_pheno.txt --trait y --env-int E1:E3 --out syn1
run --step 2 --null syn1 --mbfile syn_list.txt --num-random-snp 500 --seed $SEED --out syn2
run --step 2 --null syn1 --mbfile syn_reordered_list.txt --num-random-snp 500 --seed $SEED --out syn2r
same syn2.gamma syn2r.gamma ".gamma with a re-ordered .fam"
run --step 3 --null syn1 --gamma syn2 --mbfile syn_list.txt --out syn3
run --step 3 --null syn1 --gamma syn2 --mbfile syn_reordered_list.txt --out syn3r
same syn3.assoc syn3r.assoc ".assoc with a re-ordered .fam"
rejected "a sample missing from one .fam" --step 2 --null syn1 --mbfile syn_incomplete_list.txt --out x

echo "== random panel limits (simulated samples)"
run --step 2 --null syn1 --mbfile syn_list.txt --num-random-snp 100000 --seed $SEED --out syn2_all
used=$(awk '$1=="n_random_used"{print $2}' syn2_all.gamma)
requested=$(awk '$1=="n_random_requested"{print $2}' syn2_all.gamma)
if [ "$used" -gt 0 ] && [ "$used" -le 4000 ] && [ "$requested" = 100000 ]; then echo "ok    asking for more SNPs than exist uses them all ($used)"; else echo "FAIL  used=$used requested=$requested"; failures=$((failures + 1)); fi
rejected "no SNP passes QC" --step 2 --null syn1 --mbfile syn_list.txt --maf 0.6 --out x

echo "== --maf treats both alleles alike (simulated samples)"
python3 - <<'PY'
import numpy as np
rng = np.random.default_rng(2)
n, m = 300, 600
rare_snp = np.arange(m) % 3 == 0
dosage = rng.binomial(2, np.where(rare_snp, 0.002, 0.3), (n, m))
code = np.array([0b11, 0b10, 0b00])  # .bed codes of dosage 0, 1 and 2

def write_fileset(prefix, values):
    with open(prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"S{i:04d} S{i:04d} 0 0 0 -9\n")
    with open(prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"1 snp{j} 0 {j + 1} A G\n")
    bits = code[values].T  # SNPs x samples
    padded = np.zeros((m, (n + 3) // 4 * 4), dtype=np.uint8)
    padded[:, :n] = bits
    packed = padded[:, ::4] | padded[:, 1::4] << 2 | padded[:, 2::4] << 4 | padded[:, 3::4] << 6
    with open(prefix + ".bed", "wb") as f:
        f.write(bytes([0x6c, 0x1b, 0x01]) + packed.tobytes())

write_fileset("maf_a", dosage)
write_fileset("maf_b", 2 - dosage)  # every SNP counted through its other allele
PY
run --make-grm --bfile maf_a --out maf_a
run --make-grm --bfile maf_b --out maf_b
same maf_a.grm.N maf_b.grm.N "--maf keeps the same SNPs whichever allele is counted"
python3 - > maf_check.txt <<'PY'
import numpy as np
first = np.fromfile("maf_a.grm.bin", dtype="<f4")
second = np.fromfile("maf_b.grm.bin", dtype="<f4")
print("same" if np.allclose(first, second, rtol=1e-5, atol=1e-5) else "different")
PY
if [ "$(cat maf_check.txt)" = same ]; then echo "ok    the GRM does not depend on which allele is counted"; else echo "FAIL  GRM depends on the allele"; failures=$((failures + 1)); fi

echo "== GRM files"
num_id=$(wc -l < grm.grm.id)
if [ "$(stat -c %s grm.grm.bin)" = $((num_id * (num_id + 1) / 2 * 4)) ]; then echo "ok    .grm.bin is float32"; else echo "FAIL  .grm.bin size"; failures=$((failures + 1)); fi
if [ "$(head -c 8 grm.grm.sp.bin)" = FGMASPF4 ]; then echo "ok    .grm.sp.bin has its format header"; else echo "FAIL  .grm.sp.bin header"; failures=$((failures + 1)); fi

# The same GRM in the float64 layouts of version 1.0.
python3 - <<'PY'
import numpy as np
np.fromfile("grm.grm.bin", dtype="<f4").astype("<f8").tofile("legacy.grm.bin")
with open("grm.grm.sp.bin", "rb") as f:
    f.read(8)
    records = np.frombuffer(f.read(), dtype=[("i", "<i8"), ("j", "<i8"), ("v", "<f4")])
legacy = np.zeros(len(records), dtype=[("i", "<i8"), ("j", "<i8"), ("v", "<f8")])
legacy["i"], legacy["j"], legacy["v"] = records["i"], records["j"], records["v"]
legacy.tofile("legacy.grm.sp.bin")
PY
for suffix in id N group group.size; do cp "grm.grm.$suffix" "legacy.grm.$suffix"; done

run --step 1 --grm legacy --data "$EX/test_simu_pheno.txt" --trait pheno --env-int E1:E40 --out legacy_s1
same s1.null.bin legacy_s1.null.bin "a float64 GRM gives the same null model"
run --process-grm --reformat --grm legacy --out-fmt 3 --out dense_float32
same grm.grm.bin dense_float32.grm.bin "--reformat converts a float64 .grm.bin"
run --process-grm --reformat --sparse --grm legacy --out-fmt 3 --out sparse_float32
same grm.grm.sp.bin sparse_float32.grm.sp.bin "--reformat --sparse converts a float64 .grm.sp.bin"

mkdir regroup
cp legacy.grm.bin legacy.grm.id regroup/
run --process-grm --group --grm regroup/legacy --cut-value 0.05
same grm.grm.sp.bin regroup/legacy.grm.sp.bin "--group reads a float64 .grm.bin"
same grm.grm.group regroup/legacy.grm.group "--group finds the same groups"

run --process-grm --make-epis --prefix1 grm --prefix2 grm --out epis_float32
run --process-grm --make-epis --prefix1 grm --prefix2 legacy --out epis_mixed
same epis_float32.grm.bin epis_mixed.grm.bin "--make-epis takes one GRM of each width"
run --process-grm --out-close --grm grm --cut-value 0.05 --out close_dense
run --process-grm --out-close --grm legacy --cut-value 0.05 --out close_dense_legacy
same close_dense.grm.pruned_iid close_dense_legacy.grm.pruned_iid "--out-close reads both .grm.bin widths"
run --process-grm --out-close --sparse --grm grm --cut-value 0.05 --out close_sparse
run --process-grm --out-close --sparse --grm legacy --cut-value 0.05 --out close_sparse_legacy
same close_sparse.grm.pruned_iid close_sparse_legacy.grm.pruned_iid "--out-close --sparse reads both widths"

run --make-grm --bfile "$EX/test" --npart 2 1 --out part
run --make-grm --bfile "$EX/test" --npart 2 2 --out part
run --process-grm --merge --grm part --npart 2
python3 - > part_check.txt <<'PY'
import numpy as np
whole = np.fromfile("grm.grm.bin", dtype="<f4")
merged = np.fromfile("part.grm.bin", dtype="<f4")
print("same" if whole.size == merged.size and np.allclose(whole, merged, rtol=1e-6, atol=1e-6) else "different")
PY
if [ "$(cat part_check.txt)" = same ]; then echo "ok    a GRM built in parts matches the whole"; else echo "FAIL  GRM built in parts"; failures=$((failures + 1)); fi

# A two-column .grm.id (FID IID), as written by fastGMA and GCTA. Simulated samples.
python3 - <<'PY'
with open("syn_two.grm.id", "w") as f:
    for i in range(300):
        f.write(f"S{i:04d} S{i:04d}\n")
PY
for suffix in sp.bin N group group.size; do cp "syn_grm.grm.$suffix" "syn_two.grm.$suffix"; done
run --step 1 --grm syn_two --data syn_pheno.txt --trait y --env-int E1:E3 --out syn_two_s1
same syn1.null.bin syn_two_s1.null.bin "a two-column .grm.id gives the same null model"

head -c 1000 grm.grm.bin > short.grm.bin
cp grm.grm.id short.grm.id
rejected "a .grm.bin of the wrong size" --process-grm --group --grm short --cut-value 0.05
for suffix in id group group.size; do
    cp "grm.grm.$suffix" "partial.grm.$suffix"
    cp "grm.grm.$suffix" "noheader.grm.$suffix"
done
cp grm.grm.sp.bin partial.grm.sp.bin
printf 'x' >> partial.grm.sp.bin
rejected "a partial record in .grm.sp.bin" --step 1 --grm partial --data "$EX/test_simu_pheno.txt" --trait pheno --env-int E1:E40 --out x
tail -c +9 grm.grm.sp.bin > noheader.grm.sp.bin
rejected "float32 records without their header" --step 1 --grm noheader --data "$EX/test_simu_pheno.txt" --trait pheno --env-int E1:E40 --out x

echo "== errors"
run --step 1 --grm grm --data "$EX/test_simu_pheno.txt" --trait pheno --env-int E1:E20 --out other
rejected "gamma of another null model" --step 3 --null other --gamma s2 --bfile chr22 --out x
rejected "phenotype option in step 3" --step 3 --null s1 --gamma s2 --bfile chr22 --data x --out x
rejected "--split-task in step 2" --step 2 --null s1 --bfile chr22 --split-task 2 1 --out x
rejected "--bfile together with --mbfile" --step 2 --null s1 --bfile chr22 --mbfile chr_list.txt --out x
rejected "no genotypes" --step 2 --null s1 --out x
rejected "--split-task with part above total" --step 3 --null s1 --gamma s2 --bfile chr22 --split-task 3 7 --out x
rejected "--snp-range with --mbfile" --step 3 --null s1 --gamma s2 --mbfile chr_list.txt --snp-range a b --out x
rejected "unknown step" --step 4 --null s1 --out x

echo
if [ "$failures" -eq 0 ]; then echo "all checks passed"; else echo "$failures check(s) failed"; fi
exit "$failures"
