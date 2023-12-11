#!/usr/bin/env bash

#gen_job $cpp_bin $in_file $sub_in_file $i $size $perc 200000000 $a "binary" $h $k
gen_job(){
    echo "$1 --key_file=$2 --sub_key_file=$3 --dataset_type=$4 --total_num_keys=$5 --budget=$6 --num_action=$7 --keys_file_type=$8 --h=$9 --kde_type=${10} --init_perc=${11} --setting=2" >> $job_description
}

local(){
    /home/ubuntu/sofs/parallel-20230122/src/parallel --timeout 3000% \
        --group --resume-failed -j 90% --load 85% --progress \
        --retries 3 --joblog $1 -a $2 > $3

}

remote(){

    /home/ubuntu/sofs/parallel-20230122/src/parallel --sshloginfile `pwd`/nodelists --timeout 3000% \
        --group --resume-failed -j 90% --load 85% --progress \
        --retries 3 --joblog $1 -a $2 > $3

}


run_exps(){

    size_list=( 
        #5000000
        #20000000
        50000000
        100000000
        150000000
    )
    

    action_list=(
        #1
        2
        #4
        #6
        #8
    )
    
    budget_list=(
        5
        10
        15
        20
        25
        30
    )

    h_list=(
        1.5
        #0.5
        #1
    )

    output_file="/home/ubuntu/data/attack_out/out/graybox.out"
    job_description="/home/ubuntu/data/attack_out/jobdir/graybox.des"
    job_log="/home/ubuntu/data/attack_out/logdir/graybox.log"
    
    cpp_bin=`pwd`/build/space_aca_dn

    file_list=(
        #"/home/ubuntu/data/datasets/longlat-200M.bin.data"
        #"/home/ubuntu/data/datasets/longitudes-200M.bin.data"

        "/home/ubuntu/data/datasets/ycsb-200M.bin.data"
        "/home/ubuntu/data/datasets/lognormal-190M.bin.data"
    )

    kde_list=(
        #"gaussian"
        "tophat"  
    )
    sub_size=5



    if [ -f $output_file ]; then
        rm $output_file
    fi

    if [ -f $job_description ]; then
        rm $job_description
    fi

    if [ -f $job_log ]; then
        rm $job_log
    fi

    kde_dir="/home/ubuntu/data/attack_out/kde"


    for i in ${!file_list[@]}; do
        in_file=${file_list[$i]}

        for size in  ${size_list[@]}; do
            for h in ${h_list[@]}; do
                
            #attack with diff budget config
                for perc in ${budget_list[@]}; do
                    for a in ${action_list[@]}; do
                        for k in ${kde_list[@]}; do
                            for (( s=0; s<$sub_size; s++)); do
                                if [[ "$in_file" == *"longitudes"* ]]; then
                                    sub_in_file=${kde_dir}/longitude_${k}_${size}_${h}_${s}
                                    gen_job $cpp_bin $in_file $sub_in_file $i $size $perc $a "binary" $h $k 50
                                
                                elif [[ "$in_file" == *"longlat"* ]]; then
                                        sub_in_file=${kde_dir}/longlat_${k}_${size}_${h}_${s}
                                    gen_job $cpp_bin $in_file $sub_in_file $i $size $perc $a "binary" $h $k 50
                                
                                elif [[ "$in_file" == *"ycsb"* ]]; then
                                    sub_in_file=${kde_dir}/ycsb_${k}_${size}_${h}_${s}
                                    gen_job $cpp_bin $in_file $sub_in_file $i $size $perc $a "binary" $h $k 50
                                
                                else
                                    sub_in_file=${kde_dir}/lognormal_${k}_${size}_${h}_${s}
                                    gen_job $cpp_bin $in_file $sub_in_file $i $size $perc $a "binary" $h $k 50
                                fi
                            done
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
