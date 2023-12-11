#!/bin/bash

file_list=(
    "./datasets/longitudes-200M.bin.data"
    "./datasets/longlat-200M.bin.data"
    "./datasets/ycsb-200M.bin.data"
    "./datasets/lognormal-190M.bin.data"
)

budget_list=(0.005 0.01 0.015 0.02)

thread_num=(1 2 4 8)

insert_frac=(0.1)

sample_run=5


for i in ${!file_list[@]}; do
    in_file=${file_list[$i]}

    for j in ${insert_frac[@]}; do
        for k in ${thread_num[@]}; do
            for a in ${budget_list[@]}; do 
                for (( b=1; b<=$sample_run; b++)); do
                    if [[ "$in_file" == *"longitudes"* ]] || [[ "$in_file" == *"longlat"* ]]; then
                        ./run.sh $in_file binary double 10000000 200000000 mixed 0 apex 1 $j $a $k
                    else
                        ./run.sh $in_file binary int 10000000 200000000 mixed 0 apex 1 $j $a $k
                    fi
                done
            done
        done
    done
done

