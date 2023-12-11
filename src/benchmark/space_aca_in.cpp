// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

/*
 * Simple benchmark that runs a mixture of point lookups and inserts on ALEX.
 */


#include <algorithm>
#include <experimental/algorithm>
#include <vector>
#include <iomanip>
#include <iostream>
#include <random>
#include <chrono> //For system_clock
#include <bits/stdc++.h>
#include <map>
#include <numeric>
#include <tuple>
#include <stdlib.h>
#include <sstream>
#include <iterator>

#include "absl/strings/str_format.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/algorithms/knapsack_solver.h"

#include "flags.h"
#include "utils.h"
#include "alex_nodes.h"
#include "attack.h"                             /*Rui:our attack logics*/
#include "alex.h"
#include "leap_util.h"

// Modify these if running your own workload
//#define KEY_TYPE double
//#define PAYLOAD_TYPE double


#define KEY_TYPE unsigned long long 
#define PAYLOAD_TYPE unsigned long long

//sample szie for lookup

using namespace std;

/*
 * Required flags:
 * --keys_file              path to the file that contains keys
 * --keys_file_type         file type of keys_file (options: binary or text)
 * --init_num_keys          number of keys to bulk load with
 * --total_num_keys         total number of keys in the keys file
 * --batch_size             number of operations (lookup or insert) per batch
 *
 * Optional flags:
 * --insert_frac            fraction of operations that are inserts (instead of
 * lookups)
 * --lookup_distribution    lookup keys distribution (options: uniform or zipf)
 * --time_limit             time limit, in minutes
 * --print_batch_stats      whether to output stats for each batch
 */

template <class T>
KEY_TYPE find_duplicate_blackbox(T keys[], int init_num_keys, alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index);

KEY_TYPE find_duplicate_whitebox(alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, int dataset);

KEY_TYPE find_duplicate_graybox(alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, int dataset);

void get_dn_with_smallest_sparsity(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, 
        std::vector<std::pair<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*, std::pair<double, double>>>& ks_ss);

int main(int argc, char* argv[]) {




    //read flags
    auto flags = parse_flags(argc, argv);
    auto total_num_keys = stoi(get_required(flags, "total_num_keys")); 
    auto key_file = get_required(flags, "key_file");
    auto init_num_keys = stoi(get_required(flags, "init_num_keys"));
    auto num_duplicates = stoi(get_required(flags, "num_duplicates"));
    auto setting = stoi(get_required(flags, "setting"));
    auto dataset = stoi(get_required(flags, "dataset"));


    auto key_type = get_with_default(flags, "key_type", "binary");
    auto sub_total_num_keys = stoi(get_with_default(flags, "sub_total_num_keys", "0"));
    auto sub_key_file = get_with_default(flags, "sub_key_file", "/");


    // key arrays
    auto keys = new KEY_TYPE[total_num_keys];

    // read original dataset
    if (key_type == "binary") {
        load_binary_data(keys, total_num_keys, key_file);
    } else if (key_type == "text") {
        load_text_data(keys, total_num_keys, key_file);
    } else {
        std::cerr << "--key_type must be either 'binary' or 'text'"
            << std::endl;
        return 1;
    }


    // calcualte legitimal key index
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[init_num_keys];
    std::mt19937_64 gen_payload(std::random_device{}());
    std::mt19937_64 gen_key(std::random_device{}());


    for(int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index; 

    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });
    index.bulk_load(values, init_num_keys);


    KEY_TYPE duplicate_key;



    if(setting == 0) {
        duplicate_key = find_duplicate_blackbox(keys, init_num_keys, &index);

        //duplicate_key = -17237.8630338678849512;
        //duplicate_key = -17381.7294313226666418;
        //duplicate_key = -66.7098473016988862;
        //duplicate_key = -49.7777477515247782;
        //duplicate_key = 3.12773;
        //duplicate_key = 154993691;
    } else if(setting == 1) {
        duplicate_key = find_duplicate_whitebox(&index, dataset);
        //duplicate_key = -17237.8630338678849512;

    } else if(setting == 2){

        if(sub_total_num_keys == 0 || sub_key_file == "/") {
            std::cout<<"sub_total_num_keys or sub_key_file is wrong\n";
            return 1;
        }

        auto sub_keys = new KEY_TYPE[sub_total_num_keys];

        if (key_type == "binary") {
            load_binary_data(sub_keys, sub_total_num_keys, sub_key_file);
        } else if (key_type == "text") {
            load_text_data(sub_keys, sub_total_num_keys, sub_key_file);
        } else {
            std::cerr << "--key_type must be either 'binary' or 'text'"
                << std::endl;
            return 1;
        }

        alex::Alex<KEY_TYPE, PAYLOAD_TYPE> sub_index;
        auto sub_values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[init_num_keys];

        for(int i = 0; i < init_num_keys; i++) {
            sub_values[i].first = sub_keys[i];
            sub_values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }

        std::sort(sub_values, sub_values + init_num_keys,
                [](auto const& a, auto const& b) { return a.first < b.first; });
        sub_index.bulk_load(sub_values, init_num_keys);

        duplicate_key = find_duplicate_whitebox(&sub_index, dataset);

        //duplicate_key = -30.6549;

        delete[] sub_keys;
        delete[] sub_values;

    } else {
        std::cout<<"unrecognized setting flag\n";
        return 1;
    }

    ExpeData b;

    KEY_TYPE diff = 1;

    if(dataset == 0) {
        diff = 0.0000000000001;
    } else if(dataset == 1) {
        diff = 0.00000000001;
    }     

    for(int i = 0; i < num_duplicates; i++) {
        b.dbg_ = i;
        index.insert(duplicate_key, static_cast<PAYLOAD_TYPE>(gen_payload()), &b);
    }

    delete[] keys;
    delete[] values;
    return 0;

}


template <class T>
KEY_TYPE find_duplicate_blackbox(T keys[], int init_num_keys, alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index){

    //generate duplicates
    int key_range_sample_size = 1000;
    std::vector<KEY_TYPE> estimates;
    estimates.clear();
    std::experimental::sample(keys, keys + init_num_keys, std::back_inserter(estimates),
            key_range_sample_size, std::mt19937{std::random_device{}()});

    KEY_TYPE min_key = estimates[0];
    KEY_TYPE max_key = estimates[estimates.size() - 1];

    std::uniform_real_distribution<> dis(min_key, max_key);
    //std::uniform_int_distribution<unsigned long long> dis(min_key, max_key);

    std::mt19937 gen(std::random_device{}());


    KEY_TYPE guess_key = dis(gen);

    std::cout<< "index min_key: " << index->get_min_key_pub() << " index max_key:"<< index->get_max_key_pub() <<" sampled min_key: "<<min_key << " sampled max_key:"<<max_key<<" guess_key:"<<guess_key <<"\n";

    return guess_key;

}

KEY_TYPE get_right_ks(alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node) {
        
        int bucketID = parent->model_.predict(node->first_key());
        int dup = static_cast<int>(node->duplication_factor_);
        int repeats = 1 << dup;
        int start_bucketID = bucketID - (bucketID %repeats); 
        int num_buckets = 1 << 1;

        int end_bucketID = start_bucketID + num_buckets;
        int mid_bucketID = start_bucketID + num_buckets / 2; 

        KEY_TYPE right_kspace = (mid_bucketID - parent->model_.b_) / parent->model_.a_;

        return right_kspace;
}


void get_dn_with_smallest_sparsity(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, 
        std::vector<std::pair< alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*, std::pair<double, double>>>& ks_ss) {

     for(auto node: dns) {
        
         if(node->num_keys_ <= 0) {
            continue;
        }

        auto parent = index->get_parent(node);
        int bucketID = parent->model_.predict(node->first_key());
        int dup = static_cast<int>(node->duplication_factor_);
        int repeats = 1 << dup;
        int start_bucketID = bucketID - (bucketID %repeats); 
        int num_buckets = 1 << 1;

        int end_bucketID = start_bucketID + num_buckets;
        int mid_bucketID = start_bucketID + num_buckets / 2; 

        KEY_TYPE right_kspace = (mid_bucketID - parent->model_.b_) / parent->model_.a_;

        //left
        
        int cnt = 0;

        double right_ns = 1;
        double left_ns = 1;

        int node_cnt = 0;


        for(int i = 0; i < node->data_capacity_; i++) {
            if(node->check_exists(i) && node->get_key(i) > right_kspace) {
                left_ns = static_cast<double> (cnt) / node->num_keys_;
                cnt = 0;
                break;
            } else if(node->check_exists(i)  && node->get_key(i) < right_kspace) {
                cnt++;
            }
        }


        for(int i = node->data_capacity_ - 1; i >= 0; i--) {
            if(node->check_exists(i) && node->get_key(i) < right_kspace) {
                right_ns = static_cast<double> (cnt) / node->num_keys_;
                cnt = 0;
                break;
            } else if(node->check_exists(i)  && node->get_key(i) > right_kspace) {
                cnt++;
            }
        }

        ks_ss.push_back({node, {left_ns, right_ns}});
     }

}

KEY_TYPE find_duplicate_whitebox(alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, int dataset){
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* target_node;

    index->get_all_datanodes(dns);
    std::vector<std::pair<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*, std::pair<double, double>>> ks_ss;


    /*get_dn_with_smallest_sparsity(dns, index, ks_ss);

    std::sort(ks_ss.begin(),
            ks_ss.begin() + ks_ss.size(),
            [](auto &a, auto &b) {
                return a.second.first < b.second.first;
            });
    


    target_node = ks_ss[30000].first;*/
    

    //KEY_TYPE ks = get_right_ks(index->get_parent(target_node), target_node);


    std::sort(dns.begin(),
            dns.begin() + dns.size(),
            [](auto *a, auto *b) {

            //int num_buckets = 1 << a->duplication_factor;
            //int end_bucketID = start_bucketID + num_buckets;
            //int mid_bucketID = start_bucketID + num_buckets / 2;
            
            //return a->data_size() < b->data_size();

            //int shifts1 = a->shifts_per_insert();
            //int shifts2 = a->shifts_per_insert();

            std::pair<int, int> max1 = attack::get_max_continuous_key(a, 0, a->data_capacity_);
            std::pair<int, int> max2 = attack::get_max_continuous_key(b, 0, b->data_capacity_);

            return max1.first > max2.first; 
            //

            //return a->expected_avg_shifts_ > b->expected_avg_shifts_;

            });

    KEY_TYPE ret_key;


    //20000:403, 158
    //30000:403, 167
    //25000:413, 
    target_node = dns[0.5 * dns.size()];

    KEY_TYPE diff = 1;

    if(dataset == 0) {
        diff = 0.0000000000001;
    } else if(dataset == 1) {
        diff = 0.00000000001;
    } 

    auto max_idx = attack::get_max_continuous_key(target_node, 0, target_node->data_capacity_);

    ret_key = target_node->get_key(max_idx.second) - diff;
    std::cout<< "index min_key: " << index->get_min_key_pub() << " index max_key:"<< index->get_max_key_pub() << ", segment length:"<< max_idx.first <<"\n";
    return ret_key;

}

KEY_TYPE find_duplicate_graybox(alex::Alex<KEY_TYPE, PAYLOAD_TYPE>* index, int dataset)
{
    KEY_TYPE ret_key;

    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* target_node;

    index->get_all_datanodes(dns);
    std::sort(dns.begin(),
            dns.begin() + dns.size(),
            [](auto const* a, auto const* b) {
            return a->data_size()  < b->data_size(); 
            });

    target_node = dns[0];
    //std::uniform_real_distribution<> dis(target_node->min_key_, target_node->max_key_);
    //std::uniform_int_distribution<unsigned long long> dis(target_node->min_key_, target_node->max_key_);
    //std::mt19937 gen(std::random_device{}());

    ret_key = target_node->first_key();
    //ret_key = dis(gen);
    std::cout<< "index min_key: " << index->get_min_key_pub() << " index max_key:"<< index->get_max_key_pub() <<"\n";

    return ret_key;

}
