cd /net/zootopia/disk1/chaon/WORK/fastgxe/example/

fastgxe --make-grm --bfile test --code-type 2 --out test

fastgxe --process-grm --group --grm test.agrm --cut-value 0.05 --out test.agrm

fastgxe --test-gxe --grm test --bfile test --data test_simu_pheno.txt --trait pheno \
--env-int E1:E40 --out test_gxe




