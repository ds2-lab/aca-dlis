#!/usr/bin/env bash

# gen_job $cpp_bin $in_file $init_num_keys $size $c $b $lookperc
gen_job(){
    echo "$1 --key_file=$2 --dataset_type=$3 --init_num_keys=$4 --workload_size=$5 --c=$6 --budget=$7 --lookup_perc=$8 --sub_key_file=$9" >> $job_description
}

local(){
    /home/ubuntu/sofs/parallel-20230122/src/parallel \
        --group --resume-failed -j 85% --load 80% --progress \
        --retries 3 --joblog $1 -a $2 >> $3
}

remote(){

    /home/ubuntu/sofs/parallel-20230122/src/parallel --sshloginfile `pwd`/nodelists \
        --group --resume-failed -j 90% --load 85% --progress \
        --retries 3 --joblog $1 -a $2 > $3
}


run_exps(){

    worload_size_list=( 
        10000000
        #20000000
        #50000000
    )
        
    
    c_list=(
        200
        400
        600
    )

    b_list=(
        #0.05
        #0.1
        #0.15
        #0.2
        0.5
    )
    
    lookup_perc_list=(
        #0.5
        0.9
    )


    output_file="/home/ubuntu/data/attack_out/out/perf_aca_greedy_graybox.out"
    job_description="/home/ubuntu/data/attack_out/jobdir/perf_aca.des"
    job_log="/home/ubuntu/data/attack_out/logdir/perf_aca.log"
    
    cpp_bin=`pwd`/build/time_aca_graybox

    file_list=(
        #"/home/ubuntu/data/datasets/longlat-200M.bin.data"
        #"/home/ubuntu/data/datasets/longitudes-200M.bin.data"
        #"/home/ubuntu/data/datasets/ycsb-200M.bin.data"
        "/home/ubuntu/data/datasets/lognormal-190M.bin.data"

    )

    sub_in_file="/home/ubuntu/data/attack_out/kde/lognormal_tophat_50000000_1.5"

    init_num_keys=10000000
    sample_run=5
    data_type=0



    #config here to enable baseline experiments
    set_baseline=0


    if [ -f $job_description ]; then
        rm $job_description
    fi

    if [ -f $job_log ]; then
        rm $job_log
    fi

    
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


        #baseline: legitimate workloads
        if [[ "$set_baseline" -eq 1 ]]; then

            for size in ${worload_size_list[@]}; do
                for lookperc in ${lookup_perc_list[@]}; do
                    for (( j=1; j<=$sample_run; j++)); do
                        gen_job $cpp_bin $in_file $data_type $init_num_keys $size 0 0 $lookperc
                    done
                done
            done
        fi


        #workload with poisoning inserts

        for size in  ${worload_size_list[@]}; do
            for c in ${c_list[@]}; do
                for b in ${b_list[@]}; do
                    for lookperc in ${lookup_perc_list[@]}; do
                        for (( j=0; j<$sample_run; j++)); do
                            gen_job $cpp_bin $in_file $data_type $init_num_keys $size $c $b $lookperc ${sub_in_file}_${j}
                        done
                    done
                done
            done
        done
    done


    local $job_log $job_description $output_file



    #remote $job_log $job_description $output_file

    #./run_perf_aca_blackbox.sh

    #./shutdown_nodes.sh

    sudo poweroff
    
}

run_exps 
