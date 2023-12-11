#!/bin/bash


#thread_num=(0 1 2 4 8 16 24 32 48)

for j in 1
do
    rm -f /mnt/mem/apex.data
    LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libjemalloc.so.2 ./build/pmdk/src/PMDK/src/nondebug/libpmemobj.so.1 ./build/pmdk/src/PMDK/src/nondebug/libpmem.so.1" numactl --cpunodebind=0 --membind=0 ./build/benchmark \
        --keys_file=$1 \
        --keys_file_type=$2 \
        --keys_type=$3 \
        --init_num_keys=10000000 \
        --workload_keys=$4 \
        --total_num_keys=$5 \
        --operation=$6 \
        --insert_frac=${10} \
        --lookup_distribution=uniform \
        --theta=0.99 \
        --using_epoch=$7 \
        --thread_num=${12} \
 	    --poisoning_frac=${11} \
        --index=$8 \
        --sort_bulkload=$9
done
