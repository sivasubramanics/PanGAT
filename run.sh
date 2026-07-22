kmc_tools simple test_data/output/01.kmc/Elliot_HAP1 test_data/output/01.kmc/Elliot_HAP2 intersect Elliot -ocright
kmc_tools transform Elliot dump Elliot.tsv
for i in test_data/input/genotype_data/*;do p=$(pwd)/${i};s=$(basename ${i});for f in ${p}/*;do echo $s,$f;done;done > test_data/input/samples.csv