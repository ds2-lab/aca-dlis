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


int main(int argc, char* argv[]) {


    std::mt19937_64 gen_payload(std::random_device{}());
    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index; 


    for(int i = 0; i < 10; i++) {
        index.insert(i, static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    for(int i = 0; i < 20; i++) {
        index.insert(1, static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    cout<<"data nodes:"<<index.stats_.num_data_nodes;

    return 0;

}
