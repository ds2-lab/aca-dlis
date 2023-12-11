// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

#include "zipf.h"
#include "general_zipf.h"
#include "../util/uniform.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <experimental/algorithm>

template <class T>
bool load_binary_data(T data[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if (!is.is_open()) {
        return false;
    }
    is.read(reinterpret_cast<char*>(data), std::streamsize(length * sizeof(T)));
    is.close();
    return true;
}

template <class T>
bool load_text_data(T array[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str());
    if (!is.is_open()) {
        return false;
    }
    int i = 0;
    std::string str;
    while (std::getline(is, str) && i < length) {
        std::istringstream ss(str);
        ss >> array[i];
        i++;
    }
    is.close();
    return true;
}

template <class T>
T* get_search_keys(T array[], int num_keys, int num_searches) {
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, num_keys - 1);
    auto* keys = new T[num_searches];
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(gen);
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_zipf(T array[], int num_keys, int num_searches) {
    auto* keys = new T[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys);
    for (int i = 0; i < num_searches; i++) {
        int pos = zipf_gen.nextValue();
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_zipf_with_theta(T array[], int num_keys, int num_searches, double theta) {
    auto* keys = new T[num_searches];
    std::default_random_engine generator_;
    zipfian_int_distribution<int> dis(0, num_keys - 1, theta);
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(generator_);
        keys[i] = array[pos];
    }
    std::random_shuffle(&keys[0], &keys[num_keys - 1]);
    return keys;
}

template <class T>
void inject_poisoning_inserts(T* keys, uint64_t poisoning_inserts, uint64_t consecutive_num, uint64_t base_num_inserts, char* 
    workload, int mixed_generate_num, int dataset_type){

    int num_locations = 1;
    consecutive_num = 10000;

    std::cout<<"consecutive_num="<<consecutive_num<<"\n";
    num_locations = static_cast<uint64_t> (poisoning_inserts / consecutive_num);
    std::cout<<"num_locations="<<num_locations<<"\n";
    
    int estimates_cnt = 0;
    std::vector<T> estimates;
    typedef std::pair<int, T> OPT;
    OPT *mixed_workload = reinterpret_cast<OPT *>(workload);

    std::experimental::sample(keys, keys + base_num_inserts, std::back_inserter(estimates),
    				num_locations, std::mt19937{std::random_device{}()});

    std::cout<<"dataset_type="<<dataset_type<<"\n";
    

    T diff = 0;
    if(dataset_type == 0 ) {
        diff = 0.0000000000001;
    } else if (dataset_type == 1) {
        diff = 0.00000000001;
    } else {
        diff = 1;
    }
    std::cout<<"diff="<<diff<<"\n";


    T consecutive_key_start = estimates[estimates_cnt] + diff;
    
    int poisoning_insert_left = 0;





/*

    std::cout.precision(19);
    std::cout<< std::scientific;
    
    std::cout<<"poisoning_insert_left="<<poisoning_insert_left
        <<", poisoning_inserts="<<poisoning_inserts
        <<", num_locations="<<num_locations
        <<", estimates[estimates_cnt]="<<estimates[estimates_cnt]
        <<", mixed_generate_num="<<mixed_generate_num
        <<", diff="<<diff<<"\n";

    

    
    for(int i = 0; i < mixed_generate_num; i++) {
        if(mixed_workload[i].first == 0) {
            consecutive_key_start = mixed_workload[i].second;
               break;
        }
     
    }*/




    for(int i = 0; i < mixed_generate_num; i++) {

        if(poisoning_insert_left >= poisoning_inserts) {
            break;
        }

        if(mixed_workload[i].first == 0) {
            mixed_workload[i].second = consecutive_key_start; 
            consecutive_key_start -= diff;
            poisoning_insert_left++;

            if(poisoning_insert_left % consecutive_num == 0 && num_locations > 1) {
                estimates_cnt++;
                consecutive_key_start = estimates[estimates_cnt];
            }

        }

    }





/*
    for(int i = 0; i < mixed_generate_num; i++) {
    	std::cout<<"["<<mixed_workload[i].first << "," << mixed_workload[i].second<<"], ";
    }*/

    std::cout<<"\n";


}

template <class T>
char* generate_workload(std::string& operation_, uint64_t& generate_num, T* keys, uint64_t base_insert_num, double
                        insert_frac_, int dataset_type, uint64_t batch_size, double poisoning_frac_, 
                            int num_lookups_per_batch, std::string& distribution = "uniform", double theta = 0.99){
    char *workload = nullptr;


    if(operation_ == "mixed"){
        typedef std::pair<int, T> OP;
        unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937_64 g2 (seed1);
        uint64_t u64Random = g2();
        UniformRandom rng(u64Random);
        uint32_t random;
        uint32_t insert_sign = static_cast<uint32_t>(insert_frac_ * 100);

        uint64_t my_inserts = 0;
        uint64_t my_searches = 0;
        uint64_t mixed_generate_num = 0;
        mixed_generate_num = generate_num;
        // generate_num only means the approximate number of inserts in mixed workload 
        if(insert_sign != 0){
            
            // uncomment here to calculate the total workload divided by the insert_frac_
            //mixed_generate_num = static_cast<uint64_t>(generate_num / insert_frac_);
            

            //using multiplication to calculate the total workload size 
            mixed_generate_num = static_cast<uint64_t>(generate_num);
        }

        auto mixed_workload = new OP[mixed_generate_num];
        T* lookup_keys = nullptr;
        if(num_lookups_per_batch != 0){
            lookup_keys = get_search_keys(keys, base_insert_num, num_lookups_per_batch);
        }


        for(uint64_t i = 0; i < mixed_generate_num; ++i){
            random = rng.next_uint32() % 100;
            if (random < insert_sign) { /*insert*/
                mixed_workload[i].first = 0;
                mixed_workload[i].second = keys[base_insert_num + my_inserts];
                my_inserts++;
                if(my_inserts == generate_num) break;
            } else { /*get*/
                mixed_workload[i].first = 1;
                mixed_workload[i].second = lookup_keys[my_searches % num_lookups_per_batch];
                my_searches++;
                if (((my_searches % num_lookups_per_batch) == 0) && (num_lookups_per_batch != 0))
                {
                    delete [] lookup_keys;
                    lookup_keys = get_search_keys(keys, base_insert_num, num_lookups_per_batch);
                }
            }
        }

        if(lookup_keys != nullptr){
            delete [] lookup_keys; 
        }


        generate_num = my_inserts + my_searches;
        auto poison_insert_num = static_cast<uint64_t> (poisoning_frac_ * generate_num);



        if(poison_insert_num != 0) {
            inject_poisoning_inserts(keys, poison_insert_num, batch_size, base_insert_num, 
            reinterpret_cast<char*> (mixed_workload), generate_num, dataset_type);
        }

        std::cout << "Generated mixed operations = " << generate_num << ": #inserts=" 
            << my_inserts << "; #searches=" << my_searches <<"; #poisoning inserts="
            <<poison_insert_num<< std::endl;
        workload = reinterpret_cast<char*>(mixed_workload);
    } else if ((operation_ == "insert") || (operation_ == "debug") || (operation_ == "var_insert")) {
        // Reuse the existing key array
        auto my_workload = new T[generate_num];
        memcpy(my_workload, keys + base_insert_num, sizeof(T) * generate_num);
        workload = reinterpret_cast<char*>(my_workload);


        //inject_poisoning_inserts(keys, generate_num, batch_size, base_insert_num, workload, generate_num, dataset_type);


        if(distribution == "zipf"){
            std::cout << "Only support zipf for update/search/range-query" << std::endl;
            exit(-1);
        }
    } else if (operation_ == "erase"){
        // Need to avoid the duplicate, and ensure these keys already existis
        auto my_workload = new T[generate_num];
        double erase_frac = generate_num / static_cast<double>(base_insert_num);
        uint32_t erase_sign = static_cast<uint32_t>(erase_frac * 10);
        uint64_t erase_idx = 0;

        for(uint64_t i = 0; i < base_insert_num; i += 10){
            if((base_insert_num - i - 10) < (generate_num - erase_idx)){
                uint64_t end = i + generate_num - erase_idx;
                for (uint64_t j = i; j < end; ++j)
                {
                    my_workload[erase_idx++] = keys[j];
                }
                break;
            }

            // For every 10 keys, delete erase_sign keys
            for (uint64_t j = 0; j < erase_sign; ++j)
            {
                my_workload[erase_idx++] = keys[i + j];
            }
        }

        if (erase_idx != generate_num)
        {
            LOG_FATAL("Generate wrong number of delete keys");
        }

        std::random_shuffle(&my_workload[0], &my_workload[generate_num - 1]);
        workload = reinterpret_cast<char*>(my_workload); 
    } else if (operation_ == "range_debug"){
        auto my_workload = new T[generate_num];
        memcpy(my_workload, keys, sizeof(T) * generate_num);
        std::sort(my_workload, my_workload + generate_num,
                [](auto const& a, auto const& b) { return a < b; });
        workload = reinterpret_cast<char*>(my_workload);
    } else {
        // Update/search/range-query operation
        T *lookup_keys = nullptr;
        if(distribution == "zipf"){
            std::cout << "generate zipf workload with theta = " << theta << std::endl;
            lookup_keys = get_search_keys_zipf_with_theta(keys, base_insert_num, generate_num, theta);
        }else{
            lookup_keys = get_search_keys(keys, base_insert_num, generate_num);
        }

        workload = reinterpret_cast<char*>(lookup_keys);
    }

    return workload;
}
