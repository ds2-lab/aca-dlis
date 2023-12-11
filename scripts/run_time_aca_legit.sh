#!/usr/bin/env bash

run_exps(){

    worload_size_list=( 
        10000000
    )
        
        
    lookup_perc_list=(
        0.5
        0.9
    )

    output_file="/home/ubuntu/data/attack_out/out/perf_aca_greedy_legit.out"
    
    #legitimate logic written in whitebox.cpp file
    cpp_bin=`pwd`/build/time_aca_whitebox

    file_list=(
        #"/home/ubuntu/data/datasets/longlat-200M.bin.data"
        #"/home/ubuntu/data/datasets/longitudes-200M.bin.data"

        "/home/ubuntu/data/datasets/ycsb-200M.bin.data"
        "/home/ubuntu/data/datasets/lognormal-190M.bin.data"
    )

    init_num_keys=10000000
    sample_run=5
    data_type=0

    
    for i in ${!file_list[@]}; do
        in_file=${file_list[$i]}

        if [[ "$in_file" == *"longitudes"* ]]; then
            data_type=0
        elif [[ "$in_file" == *"longlat"* ]]; then
            data_type=1
        elif [[ "$in_file" == *"ycsb"* ]]; then
            data_type=2
        else
            data_type=3
        fi

            for size in ${worload_size_list[@]}; do
                for lookperc in ${lookup_perc_list[@]}; do
                    for (( j=1; j<=$sample_run; j++)); do
                        $cpp_bin --key_file=$in_file --dataset_type=$data_type --init_num_keys=$init_num_keys --workload_size=$size --c=0 --budget=0 --lookup_perc=$lookperc >> $output_file
                    done
                done
            done

    done
}

run_exps 
