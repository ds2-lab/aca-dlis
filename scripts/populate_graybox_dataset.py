import sys
from sklearn.neighbors import KernelDensity
import numpy as np
import pandas as pd
import os 
import time


def main():

    print('====python script====')

    double_type_datasets = {
            'longitude': '/home/ubuntu/data/datasets/longitudes-200M.bin.data',
            'longlat': '/home/ubuntu/data/datasets/longlat-200M.bin.data'
            }



    int64_type_datasets = {
            'ycsb': '/home/ubuntu/data/datasets/ycsb-200M.bin.data',
            'lognormal': '/home/ubuntu/data/datasets/lognormal-190M.bin.data'
            }

    sub_size = 1
    total_size = 200000000

    dataset_size = [
            50000000
            #100000000,
            #150000000
            ]
    hs = [0.5, 1]

    kernel_type = ['tophat']

    double_dt = np.dtype(np.float64)
    int64_dt = np.dtype(np.uint64)

    output_dir = '/home/ubuntu/data/attack_out/kde/'


    for key in double_type_datasets:

        dataset_path = double_type_datasets[key]
        X = np.fromfile(dataset_path, dtype=double_dt).reshape(-1,1)
        print('====finish loading file for', key, '=====')
        target_dataset = X[:total_size]
        
        for k in kernel_type:
            for h in hs:
                kde = KernelDensity(kernel=k, bandwidth=h).fit(target_dataset)

                for size in dataset_size:
                    for s in range(sub_size):
                        print('pupolating data for', key, 'at', 'kernel='+k,'size='+ str(size), 'h='+str(h), 'sub='+str(s))

                        sample_dataset = kde.sample(n_samples=size).flatten()

                        if len(sample_dataset) == size:
                            out_path = output_dir + key + '_' + k + '_' + str(size) +'_' + str(h) + '_' + str(s)
                            sample_dataset.astype(double_dt).tofile(out_path)
                            print('done pupolating data for', key, 'at', 'kernel='+k,'size='+ str(size), 'h='+str(h), 'sub='+str(s))



    for key in int64_type_datasets:

        dataset_path = int64_type_datasets[key]
        X = np.fromfile(dataset_path, dtype=int64_dt).reshape(-1,1)
        print('====finish loading file for', key, '=====')

        if key == 'lognormal':
            total_size = 190000000
            
        target_dataset = X[:total_size]



        for k in kernel_type:
            for h in hs:
                kde = KernelDensity(kernel=k, bandwidth=h).fit(target_dataset)

                for size in dataset_size:

                    for s in range(sub_size):
                        print('pupolating data for', key, 'at', 'kernel='+k,'size='+ str(size), 'h='+str(h), 'sub='+str(s))

                        sample_dataset = kde.sample(n_samples=size).flatten()
                       
                                               
                        if len(sample_dataset) == size:
                            out_path = output_dir + key + '_' + k + '_' + str(size) +'_' + str(h) + '_' + str(s)
                            sample_dataset.astype(int64_dt).tofile(out_path)
                            print('done pupolating data for', key, 'at', 'kernel='+k,'size='+ str(size), 'h='+str(h), 'sub='+str(s) )



if __name__ == "__main__":
    main()
