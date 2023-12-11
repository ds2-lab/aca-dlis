#pragma once


#include <vector>


struct FindBenchmarkStruc {
    int num_exp_searchs=0;
    int depth=-1;
};

struct Cost {

    double iters_cost = 0.0;
    double shifts_cost = 0.0;
    int occup = 0;
    double single_poison_key = 0.0;
    
};
struct GapInfo{
    int max = 0;
    int start_idx = 0;
    int end_idx = 0;
};

struct Stats{
    int64_t total_selected_size = 0;
    int64_t total_budget_needed = 0;
    long long total_expected_mem = 0;
    int64_t total_budget_to_trigger_all_dns = 0;
    int mck_num_keys = 0;
    long long total_mem = 0;
    int mck_solution_flag = 0;
    int data_node_size = 0;
    int exp_num_actions =0;
    int aut_num_actions =0;
    int num_data_nodes = 0;
    int num_mck_call=0;
};


struct ExpeData {
    int shift = 0;
    // 0: no split
    // split = 1: means split sideways
    // 2: split downwards
    int split = 0; 
    int depth = -1;
    int reserved_memo = 0;
    int running_memo = 0;
    int num_keys = -1;
    int keys_capacity = -1;
    int datanode_pos = -2;
    int dbg_=0;

    //0: expand and scale
    //1: expand and retrain
    //2: split sideway
    //3: split downward
    int action_flag = -1;
    int fanout_tree_flag = 0;
    double left_avg = 0.0;
    double right_avg = 0.0;
    double whole_avg = 0.0;
    int left_num_keys = 0;
    int right_num_keys = 0;

    int fail_reason = 0;
    double empirical_cost = 0.0;
    double actual_cost = 0.0;

    int i = 0;
    double a = 0.0;
    double b = 0.0;
    int boundary = 0;
    int capacity = 0;
    double expand_expected_shifts_cost = 0.0;
    double expand_expected_search_iterations_cost = 0.0;
    double split_expected_shifts_cost = 0.0;
    double split_expected_search_iterations_cost = 0.0;

    int expected_combin_shifts = 0;
    int expected_combin_search_iters = 0;
    int expected_combin_pos_occup = 0;
    int expected_left_shifts = 0;
    int expected_left_search_iters = 0;
    int expected_left_pos_occup = 0;
    int expected_right_shifts = 0;
    int expected_right_search_iters = 0;
    int expected_right_pos_occup = 0;

    int expected_seperated_insert_flag = 0;

    double single_poison_key = 0.0;

    double expand_expected_shifts = 0.0;
    double expand_expected_search_iterations = 0.0;
    double split_expected_shifts = 0.0;
    double split_expected_search_iterations = 0.0;

    double expand_traversal_cost = 0.0;
    double split_traversal_cost = 0.0;

    int action_keys = 0;
 
    std::vector<int> segments;
    std::vector<int> gaps;
    std::vector<double> keys;
};

void reset_stats(ExpeData* l) {
   l->expected_combin_search_iters = 0;
   l->expected_combin_shifts = 0;
   l->expected_left_search_iters = 0;
   l->expected_left_shifts = 0;
   l->expected_right_search_iters = 0;
   l->expected_right_shifts = 0;
   l->action_flag = 0;
   l->fail_reason = 0;
   
}

void resetExpeData(ExpeData* data) {

    data->shift=0;
    data->split=0;
    data->depth=-1;
    data->reserved_memo=0;
    data->running_memo=0;
    data->num_keys = -1;
    data->keys_capacity = -1;
    data->datanode_pos = -2;
    data->action_keys = 0;
};

void resetStats(Stats &s) {
    s.total_selected_size = 0;
    s.total_budget_needed = 0;
    s.total_expected_mem = 0;
    s.total_budget_to_trigger_all_dns=0;
}

