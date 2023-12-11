// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

/*
 * Simple benchmark that runs a mixture of point lookups and inserts on ALEX.
 */


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
#include "alex.h"

// Modify these if running your own workload
#define KEY_TYPE double 
#define PAYLOAD_TYPE double

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






int main(int argc, char* argv[]) {


    //read flags
    auto flags = parse_flags(argc, argv);
    auto total_num_keys = stoi(get_required(flags, "total_num_keys"));
    auto init_num_keys = stoi(get_required(flags, "init_num_keys"));
    //auto budget = stod(get_required(flags, "budget"));
    std::string key_type = get_required(flags, "key_type");
    auto key_file = get_required(flags, "key_file");
    auto num_poisoning_keys = stoi(get_required(flags, "num_poisoning_keys"));
    auto k = 10000;

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

    std::mt19937_64 gen_payload(std::random_device{}());



    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[init_num_keys];
    for(int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }


    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> *index = new alex::Alex<KEY_TYPE, PAYLOAD_TYPE>;

    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index->bulk_load(values, init_num_keys);


    int key_range_sample_size = 1000;


    //for Longitudes/Longlat datasets
    double diff = 0.0000000000001;

    //uncomment for Lognormal and YCSB datasets
    //int diff = 1;

    int num_locations = num_poisoning_keys / k;

    std::vector<KEY_TYPE> estimates;
    estimates.clear();
    std::experimental::sample(keys, keys + total_num_keys, std::back_inserter(estimates),
            key_range_sample_size, std::mt19937{std::random_device{}()});

    std::sort(estimates.begin(), estimates.end(),
            [](auto const& a, auto const& b) { return a < b; });


    KEY_TYPE min_key = estimates[0];
    KEY_TYPE max_key = estimates[estimates.size() - 1];
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()


    std::uniform_real_distribution<> dis(min_key, max_key);
    //std::uniform_int_distribution<unsigned long long> dis(min_key, max_key);



    //std::cout << std::setprecision(19);
    //std::cout << std::scientific;

    int total_poison_inserts = 0;
    for(int i = 0; i < num_locations; i++) {

        //KEY_TYPE key_val = 1;

        KEY_TYPE key_val = dis(gen);

        int j = 0; 


        while( j < k) {


            /*
            std::cout<< index->stats_.num_model_nodes<<","
                << index->stats_.num_model_node_expansions<<","
                << index->stats_.num_data_nodes<<","
                << index->stats_.num_downward_splits<<","
                << index->stats_.num_sideways_splits<<","
                << index->stats_.num_keys <<","
                << index->model_size() + index->data_size()<<"\n";*/

            auto ret = index->insert(key_val, static_cast<PAYLOAD_TYPE>(gen_payload()));
            key_val -= diff;

            if(ret.second) {
                j++;
                total_poison_inserts++;
            }
        }



    }



    std::cout<< index->stats_.num_model_nodes<<","
        << index->stats_.num_model_node_expansions<<","
        << index->stats_.num_data_nodes<<","
        << index->stats_.num_downward_splits<<","
        << index->stats_.num_sideways_splits<<","
        << index->stats_.num_keys <<","
        << index->model_size() + index->data_size()<<"\n";


    delete index;
    delete[] keys;
    delete[] values;

    return 0;

}
