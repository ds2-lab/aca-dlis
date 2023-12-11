#!/usr/bin/env bash


#gen_job $cpp_double_bin $double_file $output_file 1 $size $perc $a i "binary"
gen_job(){
    echo "$1 --keys_file=$2 --dataset_type=$3 --budget=$4 --total_num_keys=$5 --num_action=$6 --keys_file_type=$7 --init_perc=$8 --setting=1" >> $job_description
}

local(){
    /home/ubuntu/sofs/parallel-20230122/src/parallel --timeout 4000% \
        --group --resume-failed -j 90% --load 80% --progress \
        --retries 3 --joblog $1 -a $2 > $3

}

remote(){

    /home/ubuntu/sofs/parallel-20230122/src/parallel --sshloginfile `pwd`/nodelists --timeout 4000% \
        --group --resume-failed -j 90% --load 80% --progress \
        --retries 3 --joblog $1 -a $2 > $3

}


run_exps(){

    size_list=( 
        #500000
        #5000000
        #20000000
        50000000
        #100000000
        #150000000
    )
    
    action_list=(
        1
        2
        4
        6
    )

    init_list=(
        #10
        0
        #40
        #60
        #50
    )

    
    budget_list=(
        5
        10
        15
        20
        25
        30
    )

    

    output_file="/home/ubuntu/data/attack_out/out/whitebox_double_b.out"
    job_description="/home/ubuntu/data/attack_out/jobdir/whitebox.des"
    job_log="/home/ubuntu/data/attack_out/logdir/whitebox.log"
    
    cpp_bin=`pwd`/build/space_aca_dn

    file_list=(
        "/home/ubuntu/data/datasets/longlat-200M.bin.data"
        "/home/ubuntu/data/datasets/longitudes-200M.bin.data"

        #"/home/ubuntu/data/datasets/ycsb-200M.bin.data"
        #"/home/ubuntu/data/datasets/lognormal-190M.bin.data"
    )




    if [ -f $output_file ]; then
        rm $output_file
    fi

    if [ -f $job_description ]; then
        rm $job_description
    fi

    if [ -f $job_log ]; then
        rm $job_log
    fi


    for i in ${!file_list[@]}; do
        in_file=${file_list[$i]}
        for size in  ${size_list[@]}; do
            #attack with diff budget config
            for init_perc in ${init_list[@]}; do
                for perc in ${budget_list[@]}; do
                    for a in ${action_list[@]}; do
                        for j in {1..3..1}; do
                            if [[ "$in_file" == *"longitudes"* ]]; then
                                gen_job $cpp_bin $in_file $i $perc $size $a "binary" $init_perc 
                            elif [[ "$in_file" == *"longlat"* ]]; then
                                gen_job $cpp_bin $in_file $i $perc $size $a "binary" $init_perc
                            elif [[ "$in_file" == *"ycsb"* ]]; then
                                gen_job $cpp_bin $in_file $i $perc $size $a "binary" $init_perc
                            else
                                gen_job $cpp_bin $in_file $i $perc $size $a "binary" $init_perc
                            fi
                        done
                    done
                done
            done
        done
    done



    local $job_log $job_description $output_file

    #remote $job_log $job_description $output_file
    
    #./shutdown_nodes.sh

    #poweroff
    
}

run_exps 
