#include <iostream>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

#include <numeric>

#include "absl/strings/str_format.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/algorithms/knapsack_solver.h"
#include <iomanip>
#include <limits>
#include <type_traits>
#include <stdlib.h>


#include "alex.h"
#include "alex_fanout_tree.h"
#include "alex_base.h"
#include "leap_util.h"


//uncomment this for double key type 

//#define KEY_TYPE double 
//#define PAYLOAD_TYPE double


//uncomment this for double key type 
#define KEY_TYPE unsigned long long
#define PAYLOAD_TYPE unsigned long long
using namespace std;



namespace attack {



    /*Function decalrations*/


    bool is_integer(KEY_TYPE threshold);

    void dump_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<KEY_TYPE> &substitude_poisoning_keys);

    std::string turn_on_node_at_index(int target_index, int total_len);
    int launch_greedy_attack();

    int run_mck_solver(std::vector<int64_t> &extra_mems, 
            std::vector<int64_t> &budgets,
            std::vector<int> &picked_dn_indexes,
            int total_budget,
            int num_choice_allowed,
            Stats &s);

    void insert_left_keys(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_keys,
            std::vector<KEY_TYPE> &poison_keys,
            std::vector<double> &chosen_node_density,
            std::vector<KEY_TYPE> &substitude_poisoning_keys);

    void gen_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);

    void gen_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size);

    void  gen_int_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size);


    void gen_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            int budget);

    void gen_double_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);

    void gen_double_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size);

    int insert_keys_to_datanode(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int budget,
            int &num_actual_inserts, 
            int num_actions_chosen, 
            std::vector<KEY_TYPE> &poison_keys, 
            std::vector<KEY_TYPE> &substitude_poisoning_keys,
            int blackbox_enable);


    bool is_all_leaf(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns);

    void select_dn(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &selected_dns, 
            int large_num, int small_num, int large_threshold, int small_threshold);

    void compute_budgets(std::vector<int64_t> &budgets, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, 
            int num_actions_allowed);

    void run_knapsack_solver(std::vector<int64_t> &extra_memos,
            std::vector<int64_t> &budgets,
            int64_t capacity,
            std::vector<int> &target_dns,
            operations_research::KnapsackSolver &solver);

    int launch_batched_knapsack_attack(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_batch,
            int total_budget, 
            int num_actions_allowed,
            int key_distribution_flag,
            Stats &s,
            std::vector<double> &chosen_node_density,
            std::vector<KEY_TYPE> &substitude_poisoning_keys);

    void compute_extra_mems(std::vector<int64_t> &extra_mems, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, 
            int num_actions_allowed);

    int64_t compute_extra_mem(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node, int num_actions_allowed);

    int64_t compute_budget(std::vector<int64_t> &budgets,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node, 
            int nodeidx,
            int num_actions_allowed);


    std::pair<KEY_TYPE, KEY_TYPE> get_slot_key_boundary(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, int pos);

    int get_max_continuous_gap_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end);

    int get_mid_occupied_pos(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end);

    std::vector<int> get_pupolated_slots(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);
    std::vector<int> get_unpupolated_slots(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);


    std::pair<int,int> get_max_continuous_key(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end);

    void gen_keys_over_slots_evenly(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            std::vector<int> slots, 
            int keys_per_slot,
            KEY_TYPE poison_keys[], 
            int desired_key_size);

    std::pair<KEY_TYPE,KEY_TYPE> find_segment_neighbor_boundary(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);
    int is_populated(KEY_TYPE injected_key, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);

    bool should_examine_node(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);

    void collect_dn_stat_before_action(KEY_TYPE key, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, ExpeData *l);

    void dump_segments(std::ofstream &out, std::vector<int> &segments, int dn_idx);

    void dump_gaps(std::ofstream &out, std::vector<int> &gaps, int dn_idx);
    std::pair<KEY_TYPE, KEY_TYPE> get_datanode_key_boundaries(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE> *parent);

    void get_nodes_with_cost( std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<KEY_TYPE, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &node_cost_diff_pairs);

    void get_nodes_with_num_keys( std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &node_with_num_keys);

    void get_nodes_with_max_gap_size(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &get_nodes_with_max_gap_pairs);

    void get_node_mem_pairs();


    std::pair<KEY_TYPE, KEY_TYPE> get_node_key_range(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, const 
            alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent);


    void get_parents(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> &parents,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);

    void get_parents(std::vector<std::pair<KEY_TYPE, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &dns,
            std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> &parents,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);

    int get_left_partion_num_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent);
    int get_right_partion_num_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent);

    void dump_dn_stats(std::ofstream &dn_stat_out_file, ExpeData *l, 
            int i, int desired_key_size, int total_insert_keys,
            int data_capacity, std::pair<KEY_TYPE, KEY_TYPE> *model);

    void merge_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, ExpeData *l);

    void dump_keys(std::ofstream &out, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, int i);

    void dump_split_segment_distro(std::ofstream &out, ExpeData *l, int dn_idx);




    bool should_insert_left_range(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* dn, 
            alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            std::pair<KEY_TYPE, KEY_TYPE> &boundary,
            int counts[]);

    std::vector<int> get_node_segments(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);

    std::vector<int> get_node_gaps(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);



    void get_poison_key_boundaries_gaps(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start_pos, int end_pos, std::pair<KEY_TYPE, KEY_TYPE> &poison_key_boundaries);


    void get_poison_key_boundaries_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start_pos, int end_pos, std::pair<KEY_TYPE, KEY_TYPE> &poison_key_boundaries);

    void encode_actions(int num_actions,
            int num_nodes,
            std::string actions_chosen,
            int current_action,
            std::vector<std::string> &res);

    bool validate_keys(KEY_TYPE keys[], int size,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);

    int get_max_continuous_key_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end,std::pair<int,int> &pos);

    KEY_TYPE get_avg_continuous_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start, int end);


    std::pair<int, int> get_max_gap(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent,
            std::pair<int, int> &left_pos, std::pair<int, int> &right_pos);

    std::pair<KEY_TYPE,KEY_TYPE> get_avg_gap_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent);


    void dump_vector_to_file(std::string path, std::vector<int>);

    void datanode_to_string(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);



    void print_vector_elements(std::vector<int64_t> &vec, std::string name, int num_element);
    void print_vector_elements(std::vector<int> &vec, std::string name, int num_element);

    void print_array(KEY_TYPE keys[], int start ,int end);

    void print_error(std::string s);

    void to_string(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> nodes, 
            std::vector<int> extra_memory, 
            std::vector<int> budget,
            int num_actions_allowed,
            int num_nodes);
    void print_index_stats(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index);
    void print_datanode_stats(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node);


    void print_execution_stats(std::string s);




    int insert(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            KEY_TYPE keys[],
            int size, int node_id,
            ExpeData *l);

    int compute_node_mems_and_budgets_on_actions(
            std::vector<std::string> &actions_str,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int cur_capacity, 
            int num_actions,
            int cur_inserted_num_keys,
            std::vector<int> &extra_mems,
            std::vector<int> &budgets);

    int get_new_data_capacity(int prev_num_keys);

    KEY_TYPE get_new_expendsion_threshold(int new_data_capacity, int prev_num_keys);

    int round_up(double expandsion_threshould);

    int get_budget(int current_num_keys, int new_data_capacity);

    int get_extra_mem(int cur_capacity, int old_capacity);

    void compute_dns_mems_and_budgets_base_on_num_actions(std::vector<std::string> &actions_str, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns, 
            std::vector<int> &extra_mems, 
            std::vector<int> &budgets);


    void gen_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<int> &picked_dn_indexes,
            std::vector<int64_t> &budgets,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_actions_allowed,
            int key_distribution_flag,
            int total_budget);

    int insert_poison_keys(std::vector<KEY_TYPE> &poison_keys, 
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int budget);

    std::vector<std::string> actions_decoder(std::string str_in);

    std::pair<int,int> compute_mems_budgets(std::string actions, 
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node);
    void compute_node_memos_and_budgets(
            std::vector<std::string> &actions_str,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            std::vector<int> &extra_mems,
            std::vector<int> &budgets);



    /*Functions implemetations begins from this line*/

    void dump_vector_to_file(std::string out_file_path, std::vector<int> v)
    {
        std::ofstream fd;
        fd.open(out_file_path, std::fstream::trunc);
        int size = v.size();
        for(int i = 0; i < size; i++) {
            fd << v[i] << std::endl;
        }
        fd.close();

    }



    void gen_keys_with_norm_distribution(KEY_TYPE array[], int num_keys, KEY_TYPE mean, KEY_TYPE stddev)
    {
        std::random_device rd;
        std::default_random_engine eng(rd());
        std::normal_distribution<double> distribution(mean, stddev);

        for(int i = 0; i < num_keys; i++) {
            array[i] = distribution(eng);
        }
    }

    void gen_keys_with_lognorm_distribution(KEY_TYPE array[], int num_keys, KEY_TYPE mean, KEY_TYPE stddev)
    {

        std::random_device rd;
        std::default_random_engine eng(rd());
        std::lognormal_distribution<double> distribution(mean, stddev);

        for(int i = 0; i < num_keys; i++) {
            array[i] = distribution(eng);
        }
    }

    void gen_keys_with_uniform_distribution(KEY_TYPE array[], int num_keys, KEY_TYPE min, KEY_TYPE max)
    {

        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_real_distribution<double> distribution(min, max);

        for(int i = 0; i < num_keys; i++) {
            array[i] = distribution(eng);
        }
    }

    int launch_batched_knapsack_attack(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_batch,
            int total_budget, 
            int num_actions_allowed,
            int key_distribution_flag,
            Stats &s,
            std::vector<double> &chosen_node_density,
            std::vector<KEY_TYPE> &substitude_poisoning_keys)
    {
        std::vector<int64_t> extra_memos;
        std::vector<int> target_dns;
        std::vector<int64_t> budgets;
        std::vector<KEY_TYPE> poison_keys;
        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;


        int num_keys_per_batch = total_budget / num_batch;
        int num_keys_before_inserts = index.stats_.num_keys;
        int num_actual_inserts = 0;

        for(int i= 0; i < num_batch; i++) {
            dns.clear();
            extra_memos.clear();
            budgets.clear();
            target_dns.clear();
            poison_keys.clear();     

            index.get_all_datanodes(dns);

            compute_extra_mems(extra_memos, dns, num_actions_allowed);
            compute_budgets(budgets, dns, num_actions_allowed);

            if(i == 0) {
                for(auto b:budgets) {
                    s.total_budget_to_trigger_all_dns += b;
                }
            }

            if(i == num_batch - 1) {
                num_keys_per_batch = total_budget - num_actual_inserts;
            }

            run_mck_solver(extra_memos, budgets, target_dns, num_keys_per_batch, num_actions_allowed,s);
            //run_knapsack_solver(extra_memos, budgets, num_keys_per_batch, target_dns, solver);


            for(auto index:target_dns) {
                s.total_selected_size += dns[index]->data_size();
                s.total_budget_needed += budgets[index];
                s.total_expected_mem += extra_memos[index];
                chosen_node_density.push_back(dns[index]->num_keys_ / static_cast<double> (dns[index]->data_capacity_));
            }

            gen_poison_keys(poison_keys, target_dns, 
                    budgets, dns, index, num_actions_allowed, key_distribution_flag, num_keys_per_batch);
            num_actual_inserts += insert_poison_keys(poison_keys, index, num_keys_per_batch);

        }



        return num_actual_inserts;
    }

    void run_knapsack_solver(std::vector<int64_t> &extra_memos,
            std::vector<int64_t> &budgets,
            int64_t batch_budget,
            std::vector<int> &target_dns,
            operations_research::KnapsackSolver &solver) 
    {
        std::vector<int64_t> capacity = {batch_budget};
        solver.Init(extra_memos, {budgets}, capacity);
        int64_t computed_value = solver.Solve();

        for (std::size_t i = 0; i < extra_memos.size(); ++i) {
            if (solver.BestSolutionContains(i)) {
                target_dns.push_back(i);
            }
        }

    }

    void get_nodes_with_cost( std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<KEY_TYPE, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &node_cost_diff_pairs)
    {
        int size = nodes.size();
        for(int i = 0; i < size; i++) {
            node_cost_diff_pairs.push_back({nodes[i]->empirical_cost() - nodes[i]->cost_, nodes[i]});
        }
    }

    void get_nodes_with_num_keys( std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &node_with_num_keys)
    {
        int size = nodes.size();
        for(int i = 0; i < size; i++) {
            node_with_num_keys.push_back({nodes[i]->num_keys_, nodes[i]});
        }
    }

    int insert_poison_keys(std::vector<KEY_TYPE> &poison_keys, 
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int budget)
    {
        std::mt19937_64 gen_payload(std::random_device{}());
        int total_inserts = 0;
        for(int i = 0; i < poison_keys.size(); i++) {
            auto ret = index.insert(poison_keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
            total_inserts++;
        }
        return total_inserts;
    }


    bool is_integer(KEY_TYPE threshold) 
    {
        int a = threshold;

        KEY_TYPE zeros = threshold - a;

        if(zeros > 0) {
            return false;
        }

        return true;
    }

    std::vector<std::string> actions_decoder(std::string str_in)
    {
        std::vector<std::string> res;
        std::stringstream ss(str_in);

        while (ss.good()) {
            std::string substr;
            getline(ss, substr, ',');
            res.push_back(substr);
        }
        return res;

    }

    std::pair<int,int> compute_mems_budgets(std::string actions, 
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node)
    {   


        std::pair<int,int> res;
        /*
           std::vector<std::string> act_strs = actions_decoder(actions);

           std::vector<int> nodes; nodes varible holds total datanodes generated by this an actin
           nodes.push_back();


           int size = act_strs.size();
           for(int i = 0; i < size; i++) { 
           auto act_str = act_strs[i];
           char action = act_str[act_str.size() - 1];
           }*/

        return res;
    }

    void compute_node_memos_and_budgets(
            std::vector<std::string> &actions_str,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            std::vector<int> &extra_mems,
            std::vector<int> &budgets)
    {
        int length = actions_str.size();
        std::pair<int,int> res;

        for(int i = 0; i < length; i++) {
            res = compute_mems_budgets(actions_str[i], node);
            extra_mems.push_back(res.first);
            budgets.push_back(res.second);
        }

    }

    int compute_node_mems_and_budgets_on_actions(
            std::vector<std::string> &actions_str,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            int cur_capacity, 
            int num_actions,
            int cur_inserted_num_keys,
            std::vector<int> &extra_mems,
            std::vector<int> &budgets) 
    {
        //simulate the current fullness of a given datanode's capacity
        KEY_TYPE prev_expansion_threshold = 0.8 * cur_capacity;
        int prev_num_keys = 0;
        int prev_capacity = cur_capacity;

        prev_num_keys = round_up(prev_expansion_threshold);

        //calculating new stats base on current stats after num_actions
        KEY_TYPE expandsion_threshould = 0.0;
        int new_data_capacity = 0;
        int new_num_keys = 0;

        for(int i = 0; i < num_actions; i++) {
            new_data_capacity = get_new_data_capacity(prev_num_keys); 
            expandsion_threshould = get_new_expendsion_threshold(new_data_capacity, prev_num_keys);
            new_num_keys = round_up(expandsion_threshould);
            extra_mems.push_back(get_extra_mem(prev_capacity, new_data_capacity));
            budgets.push_back(get_budget(cur_inserted_num_keys, new_num_keys));

            //update prev stats
            prev_expansion_threshold = expandsion_threshould;
            prev_num_keys = new_num_keys;
            prev_capacity = new_data_capacity;
        }
        return new_data_capacity;
    }

    int get_new_data_capacity(int prev_num_keys) {
        return std::max(static_cast<int>(prev_num_keys / 0.6), prev_num_keys + 1);
    }

    KEY_TYPE get_new_expendsion_threshold(int new_data_capacity, int prev_num_keys){
        return  std::min(std::max(static_cast<KEY_TYPE>(new_data_capacity * 0.8), static_cast<KEY_TYPE>(prev_num_keys + 1)),
                static_cast<KEY_TYPE>(new_data_capacity));
    }

    int round_up(double expandsion_threshould) {

        if(is_integer(expandsion_threshould)){
            return static_cast<int> (expandsion_threshould);
        } 

        return static_cast<int>(expandsion_threshould + 1);
    }

    int get_budget(int current_num_keys, int num_keys_needed) {
        return num_keys_needed - current_num_keys;
    }



    int get_extra_mem(int cur_capacity, int new_capacity){
        int old_mem_bytes = cur_capacity * sizeof(KEY_TYPE) + 
            cur_capacity * sizeof(PAYLOAD_TYPE) + static_cast<size_t>(std::ceil(cur_capacity / 64.));

        int new_mem_bytes = new_capacity * sizeof(KEY_TYPE) + 
            new_capacity * sizeof(PAYLOAD_TYPE) + static_cast<size_t>(std::ceil(new_capacity / 64.));


        return  new_mem_bytes - old_mem_bytes;
    }

    void compute_dns_mems_and_budgets_base_on_num_actions(std::vector<std::string> &actions_str, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns, 
            std::vector<int> &extra_mems, 
            std::vector<int> &budgets)
    {

        for(auto node: dns){
            compute_node_memos_and_budgets(actions_str, node, extra_mems, budgets); 
        } 

    }


    int run_mck_solver(std::vector<int64_t> &extra_mems, 
            std::vector<int64_t> &budgets,
            std::vector<int> &picked_dn_indexes,
            int total_budget,
            int num_choice_allowed,
            Stats &s)
    {



        //print_execution_stats("===before run_mck_solver()===");
        const int num_items = static_cast<int>(extra_mems.size());
        std::vector<int> all_items(num_items);
        std::iota(all_items.begin(), all_items.end(), 0);

        const std::vector<int> bin_capacities = {{total_budget}}; //15%, 20%, 30%
        const int num_bins = static_cast<int>(bin_capacities.size());
        std::vector<int> all_bins(num_bins);
        std::iota(all_bins.begin(), all_bins.end(), 0);

        operations_research::sat::CpModelBuilder cp_model;


        std::map<int, operations_research::sat::BoolVar> x;
        for (int i : all_items) {
            x[i] = cp_model.NewBoolVar().WithName(absl::StrFormat("x_%d_%d", i));
        }

        // Constraints.
        // choose exactly 1 item out of 3 consecutive ones.
        // 0-2, 3-5, 6-8, 9-11, 12-14
        // 5, 3->15
        // 1, 3: 20, 40, 80
        for (int i = 0; i < num_items; i+=num_choice_allowed) {
            std::vector<operations_research::sat::BoolVar> copies;

            for(int j = i; j < i + num_choice_allowed; j++) {
                copies.push_back(x[j]);
            }
            cp_model.AddAtMostOne(copies);
        }

        operations_research::sat::LinearExpr bin_weight;

        for (int i : all_items) {
            bin_weight += x[i] * budgets[i];
        }
        cp_model.AddLessOrEqual(bin_weight, bin_capacities[0]);

        // Objective.
        // Maximize total value of packed items.
        operations_research::sat::LinearExpr objective;
        for (int i : all_items) {
            objective += x[i] * extra_mems[i];
        }
        cp_model.Maximize(objective);



        operations_research::sat::Model model;

        operations_research::sat::SatParameters parameters;
        parameters.set_max_time_in_seconds(400.0);
        //parameters.set_num_search_workers(6);
        model.Add(NewSatParameters(parameters));


        const operations_research::sat::CpSolverResponse response 
            = SolveCpModel(cp_model.Build(), &model);


        //response.status() == operations_research::sat::CpSolverStatus::FEASIBLE

        if (response.status() == operations_research::sat::CpSolverStatus::OPTIMAL || 
                response.status() == operations_research::sat::CpSolverStatus::FEASIBLE) {

            long long total_mem = 0.0;
            for (int b : all_bins) {
                int bin_budget = 0.0;
                long long bin_mem = 0.0;
                for (int i : all_items) {
                    if (SolutionIntegerValue(response, x[i]) > 0) {
                        picked_dn_indexes.push_back(i);
                        bin_budget += budgets[i];
                        bin_mem += extra_mems[i];   
                    }
                }

                total_mem += bin_mem;
            }

            s.total_expected_mem = total_mem;

            if(response.status() == operations_research::sat::CpSolverStatus::OPTIMAL) {
                s.mck_solution_flag = 1;

            } else {

                s.mck_solution_flag = 2;

            }

            return 0;

        } else {
            s.mck_solution_flag = 0;

            exit(EXIT_FAILURE);
        }

    }


    void print_vector_elements(std::vector<int64_t> &vec, std::string name, int num_element)
    {
        std::ostringstream buffer;
        std::copy(vec.begin(), vec.begin() + num_element - 1,
                std::ostream_iterator<int>(buffer, ", "));
        buffer << vec.back();

        std::cout<<name + "={" <<buffer.str()<<"}"<<" total_size="<<vec.size()<<std::endl;

    }

    void print_vector_elements(std::vector<int> &vec, std::string name, int num_element)
    {
        std::ostringstream buffer;
        std::copy(vec.begin(), vec.begin() + num_element - 1,
                std::ostream_iterator<int>(buffer, ", "));
        buffer << vec.back();

        std::cout<<name + "={" <<buffer.str()<<"}"<<" total_size="<<vec.size()<<std::endl;

    }

    void gen_double_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<int> &picked_dn_indexes,
            std::vector<int64_t> &budgets,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_actions_allowed,
            int key_distribution_flag,
            int total_budget)
    {
        //print_execution_stats("=====gen_double_type_poison_keys()====");
        //std::cout<<"total_budget="<<total_budget<<"\n";

        int node_idx;
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node;
        std::pair<KEY_TYPE, KEY_TYPE> node_key_boundary;
        std::random_device rd;
        std::default_random_engine eng(rd());
        //std::mt19937 eng(rd());

        int keys_left = total_budget;

        double mean = 0.0;
        double attack_key = 0.0;


        if(key_distribution_flag == 1){
            for(int i = 0; i < picked_dn_indexes.size(); i++) {
                node_idx = picked_dn_indexes[i] / num_actions_allowed;
                node = dns[node_idx]; 
                node_key_boundary = get_datanode_key_boundaries(node, 
                        index.get_parent(node));

                mean = (node_key_boundary.first + node_key_boundary.second) / 2;

                // simulate keys follow norm distribution
                std::normal_distribution<double> norm_distr(mean, 1.0);
                std::uniform_real_distribution<double> uniform_distr(node_key_boundary.first,
                        node_key_boundary.second);

                keys_left -= budgets[picked_dn_indexes[i]];
                for(int j = 0; j < budgets[picked_dn_indexes[i]]; j++) {
                    attack_key = norm_distr(eng);

                    if(node_key_boundary.first <= attack_key 
                            && attack_key < node_key_boundary.second){
                        poison_keys.push_back(attack_key);
                    } else {
                        poison_keys.push_back(uniform_distr(eng));
                    }
                }

                if(i == picked_dn_indexes.size() - 1 && keys_left > 0) {
                    for(int j = 0; j < keys_left; j++) {
                        if(node_key_boundary.first <= attack_key 
                                && attack_key < node_key_boundary.second){
                            poison_keys.push_back(attack_key);
                        } else {
                            poison_keys.push_back(uniform_distr(eng));
                        }

                    }
                }

            }

        }else{

            //attack keys follow uniform distribution
            for(int i = 0; i < picked_dn_indexes.size(); i++) {
                node_idx = picked_dn_indexes[i] / num_actions_allowed;
                node = dns[node_idx]; 
                node_key_boundary = get_datanode_key_boundaries(node, index.get_parent(node));
                std::uniform_real_distribution<double> distr(node_key_boundary.first, 
                        node_key_boundary.second);

                keys_left -= budgets[picked_dn_indexes[i]];
                for(int j = 0; j < budgets[picked_dn_indexes[i]]; j++) {
                    poison_keys.push_back(distr(eng));
                }

                if(i == picked_dn_indexes.size() - 1 && keys_left > 0) {
                    //use_left_over_keys_greedyly(index, keys_left, poison_keys);
                }   

            }



        } 


    }

    void gen_int_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<int> &picked_dn_indexes,
            std::vector<int64_t> &budgets,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_actions_allowed,
            int key_distribution_flag,
            int total_budget)
    {

        //print_execution_stats("===gen_int_type_poison_keys()===");
        int node_idx;
        std::pair<KEY_TYPE, KEY_TYPE> node_key_boundary;
        std::random_device rd;
        std::default_random_engine eng(rd());

        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node;

        int expense = 0;

        int64_t mean = 0;
        int64_t attack_key = 0;

        if(key_distribution_flag == 1){
            for(auto i:picked_dn_indexes) {
                //print_execution_stats("===gen_int_type_poison_keys()===1==");
                node_idx = i / num_actions_allowed;
                node = dns[node_idx]; 
                node_key_boundary = get_datanode_key_boundaries(node, 
                        index.get_parent(node));

                mean = (node_key_boundary.first + node_key_boundary.second) / 2;

                // simulate keys follow norm distribution
                std::binomial_distribution<int> norm_distr(mean, 1.0);
                std::uniform_int_distribution<int> uniform_distr(node_key_boundary.first,
                        node_key_boundary.second);
                expense += budgets[i];

                for(int j = 0; j < budgets[i]; j++) {
                    attack_key = uniform_distr(eng);

                    /*if(node_key_boundary.first <= attack_key 
                      && attack_key < node_key_boundary.second){
                      poison_keys.push_back(attack_key);
                      } else {
                      poison_keys.push_back(uniform_distr(eng));
                      }*/
                }

            }
        }else{

            //attack keys follow uniform distribution
            for(auto i:picked_dn_indexes) {
                node_idx = i / num_actions_allowed;
                node = dns[node_idx]; 
                node_key_boundary = get_datanode_key_boundaries(node, index.get_parent(node));
                std::uniform_int_distribution<int> distr(node_key_boundary.first, 
                        node_key_boundary.second);

                expense += budgets[i];

                for(int j = 0; j < budgets[i]; j++) {
                    poison_keys.push_back(distr(eng));
                }
            }


        } 


        //check if reach total budget
        if(expense < total_budget) {
            int keys_left = total_budget - expense;

            node = dns[dns.size() / 2];

            node_key_boundary = get_datanode_key_boundaries(node, index.get_parent(node));
            std::uniform_int_distribution<int> uniform_distr(node_key_boundary.first, 
                    node_key_boundary.second);
            for(int i = 0; i < keys_left; i++) {
                poison_keys.push_back(uniform_distr(eng));
            }
        }

    }

    void gen_double_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size)
    {


        //std::cout<<"gen_double_type_poison_keys() 4 paras\n";
        std::random_device rd;
        //int seed = 1;
        std::default_random_engine eng(rd());


        std::pair<double, double> node_key_boundary = get_datanode_key_boundaries(node, 
                index.get_parent(node));

        std::uniform_real_distribution<double> uniform_distr(node_key_boundary.first,
                node_key_boundary.second);


        if(size != 0 ) {

            for(int i = 0; i < size; i++){
                poison_keys.push_back(uniform_distr(eng));
            }

            return;
        }


        if(node->num_keys_ >= node->expansion_threshold_) {
            size = 1;
        } else {
            size = ceil(node->expansion_threshold_) - node->num_keys_ + 1;
        }

        for(int i = 0; i < size; i++){
            poison_keys.push_back(uniform_distr(eng));
        }

    }

    void gen_int_type_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size)
    {

        std::random_device rd;
        std::default_random_engine eng(rd());

        std::pair<KEY_TYPE, KEY_TYPE> node_key_boundary = get_datanode_key_boundaries(node, 
                index.get_parent(node));

        std::uniform_int_distribution<unsigned long long> uniform_distr(node_key_boundary.first,
                node_key_boundary.second);

        if(size != 0) {

            for(int i = 0; i < size; i++){
                poison_keys.push_back(uniform_distr(eng));
            }

            return;
        }


        if(node->num_keys_ >= node->expansion_threshold_) {
            size = 1;
        } else {
            size = ceil(node->expansion_threshold_) - node->num_keys_ + 1;
        }

        for(int i = 0; i < size; i++){
            poison_keys.push_back(uniform_distr(eng));
        }

    }


    void gen_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int size)
    {
        if(std::is_floating_point<KEY_TYPE>::value){
            gen_double_type_poison_keys(poison_keys,
                    index,
                    node,
                    size);
        }else{
            gen_int_type_poison_keys(poison_keys,
                    index,
                    node,
                    size);

        }

    }

    void gen_poison_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<int> &picked_dn_indexes,
            std::vector<int64_t> &budgets,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_actions_allowed,
            int key_distribution_flag,
            int total_budget)
    {

        if(std::is_floating_point<KEY_TYPE>::value){
            gen_double_type_poison_keys(poison_keys,
                    picked_dn_indexes,
                    budgets,
                    dns,
                    index,
                    num_actions_allowed, 
                    key_distribution_flag, 
                    total_budget);
        }else{


            gen_int_type_poison_keys(poison_keys,
                    picked_dn_indexes,
                    budgets,
                    dns,
                    index,
                    num_actions_allowed, 
                    key_distribution_flag, 
                    total_budget);

        }

    }



    void print_array(KEY_TYPE keys[], int start ,int end)
    {
        std::ostringstream os;
        os<<"{";
        for(int i = start; i < end - 2; i++) {
            os << keys[i] <<",";
        }
        os<<keys[end - 1]<<"}";

        std::string str(os.str());

        std::cout<<str<<"\n";
    }



    int64_t compute_extra_mem(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node, int num_actions_allowed)
    {

        int current_mem = node->data_size();
        int total_mem = current_mem;

        for(int i = 0; i < num_actions_allowed; i++) {
            total_mem = total_mem * 0.8 / 0.6;
        }
        return total_mem - current_mem;
    }

    void compute_extra_mems(std::vector<int64_t> &extra_mems, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, 
            int num_actions_allowed)
    {
        int size = dns.size();
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < num_actions_allowed; j++) {
                extra_mems.push_back(compute_extra_mem(dns[i], j + 1));
            }
        }

    }

    void compute_budgets(std::vector<int64_t> &budgets, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns, 
            int num_actions_allowed)
    {
        //print_execution_stats("inside compute_budgets()");
        int size = dns.size();
        for(int i = 0; i < size; i++) {
            for(int j = 0; j < num_actions_allowed; j++){
                budgets.push_back(compute_budget(budgets, dns[i], i, j + 1));
            }
        }
    }

    int64_t compute_budget(std::vector<int64_t> &budgets,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            int nodeidx,
            int num_actions_allowed)
    {
        if (num_actions_allowed == 1) {
            if(node->expansion_threshold_ <= node->num_keys_) {
                return 1;
            }

            return ceil(node->expansion_threshold_) - node->num_keys_ + 1;
        }

        return (budgets[budgets.size() - 1] + node->num_keys_ - 1) / 0.6 * 0.8 - node->num_keys_ + 1;

    }

    bool is_all_leaf(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns) {
        int size = dns.size();

        for(int i = 0; i < size; i++) {
            if(!dns[i]->is_leaf_){
                return false;
            }
        }

        return true;
    }



    void select_dn(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &selected_dns, 
            int large_num, int small_num, int large_threshold, int small_threshold)
    {
        int current_large = 0;
        int current_small = 0;
        int size;

        for(auto datanode: dns) {
            if(current_large + current_small == large_num + small_num) {
                break;
            }

            size = datanode->data_size();
            if(current_large < large_num && size > large_threshold) {
                selected_dns.push_back(datanode);
                current_large++;
            } else if (current_small < small_num && size < small_threshold){
                selected_dns.push_back(datanode);
                current_small++;
            }
        }
    }

    //this function insert keys to a datanode 
    //adjust where to insert keys after each action
    int insert_keys_to_datanode(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int budget,
            int &num_actual_inserts, 
            int num_actions_chosen, 
            std::vector<KEY_TYPE> &poison_keys,
            std::vector<KEY_TYPE> &substitude_poisoning_keys,
            int blackbox_enable)
    {


        int ret = 0;

        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *prev, *next;
        prev = node->prev_leaf_;
        next = node->next_leaf_;


        vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;



        // root is the datanode
        if(prev == nullptr && next == nullptr) {
            for(int j = 0; j < num_actions_chosen; j++) {
                dns.clear();
                index.get_all_datanodes(dns);

                for(auto node: dns){
                    poison_keys.clear();


                    gen_poison_keys(poison_keys, index, node, 0);
                    dump_keys(poison_keys, substitude_poisoning_keys);
                    ret += insert_poison_keys(poison_keys, index, poison_keys.size());
                }

            }

            return ret;

        }






        // root is not datanode
        // we have more than one datanodes
        
        if(prev != nullptr && next != nullptr) {
            for(int j = 0; j < num_actions_chosen; j++) {
                dns.clear();

                node = next->prev_leaf_;

                while(node != prev) {
                    dns.push_back(node);
                    node = node->prev_leaf_;
                }

                for(auto node: dns){
                    poison_keys.clear();
                    gen_poison_keys(poison_keys, index, node, 0);

                    //dump substitude keys into a vector to insert into target index
                    
                    dump_keys(poison_keys, substitude_poisoning_keys);
                    ret += insert_poison_keys(poison_keys, index, poison_keys.size());

                }
            }

        } else if(prev == nullptr) {

            for(int j = 0; j < num_actions_chosen; j++) {
                dns.clear();
                node = next->prev_leaf_;

                while(node != nullptr) {
                    dns.push_back(node);
                    node = node->prev_leaf_;
                }

                for(auto node: dns){
                    poison_keys.clear();
                    gen_poison_keys(poison_keys, index, node, 0);

                    //dump substitude keys into a vector to insert into target index
                    dump_keys(poison_keys, substitude_poisoning_keys);

                    ret += insert_poison_keys(poison_keys, index, poison_keys.size());

                }
            }
            
        } else if(next == nullptr){

            dns.clear();
            node = prev->next_leaf_;

            while(node != nullptr) {
                dns.push_back(node);
                node = node->next_leaf_;
            }

            for(auto node: dns){
                poison_keys.clear();
                gen_poison_keys(poison_keys, index, node, 0);

                //dump substitude keys into a vector to insert into target index
                dump_keys(poison_keys, substitude_poisoning_keys);

                ret += insert_poison_keys(poison_keys, index, poison_keys.size());

            }

        
        } 

        return ret;

    }
    void compute_budgets_performance(std::vector<int64_t> &inarr, 
            int c,
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns,
            int num_choice) {

        for(auto node: dns) {
            for(int i = 0; i < num_choice; i++) {
                int budget = 0;

                //int num_keys_need = ceil(node->expansion_threshold_) - node->num_keys_ + 1;
                if(i % num_choice == 0) {
                    budget = c;
                } else {
                    budget = ceil(node->expansion_threshold_) - node->num_keys_ + 1;
                }
                inarr.push_back(budget);
            }
        }
    }

    void compute_expected_latencies(std::vector<int64_t> &inarr, 
            std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns,
            int num_choice) {

        for(auto node: dns) {
            for(int i = 0; i < num_choice; i++) {
                if(i % num_choice == 0) {
                    inarr.push_back(0.7 * 100 * node->num_keys_);
                } else {
                    inarr.push_back(0.3 * 100 * node->num_keys_);
                }
            }
        }

    }

    
    int launch_performance_mck_attack(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, 
            int total_budget, 
            std::vector<int> &picked_dn_indexes,
            std::vector<int64_t> &budgets,
            vector<int64_t> &expected_latencies,
            int c) {
        
        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node;


        dns.clear();
        index.get_all_datanodes(dns);

        compute_budgets_performance(budgets, c, dns, 2);
        compute_expected_latencies(expected_latencies, dns, 2);

        if(budgets.size() == 0 || expected_latencies.size() == 0 ) {
            return 1;
        }


        Stats s;

        int flag = run_mck_solver(expected_latencies, budgets, picked_dn_indexes, total_budget, 2, s);

        if(picked_dn_indexes.size() == 0) {
            return 1;
        }
        
        return 0;
        
    }



    int launch_mck_attack(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int total_budget, 
            int num_actions_allowed,
            int key_distribution_flag,
            Stats &s,
            std::vector<double> &chosen_node_density,
            std::vector<KEY_TYPE> &substitude_poisoning_keys,
            int blackbox_enable) 
    {
        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;


        dns.clear();
        index.get_all_datanodes(dns);
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node;

        std::vector<KEY_TYPE> poison_keys;      /*holds the poisoning keys*/
        std::vector<int64_t> extra_mems;            /*vector contains the extra memories after x number actions*/
        std::vector<int64_t> budgets;               /*vector contains the budgets needed for the x number actions*/
        std::vector<int> picked_dn_indexes;     /*this vector contains the item index picked by the mck solver*/



        compute_budgets(budgets, dns, num_actions_allowed);
        //cout<<"compute_budgets()\n";
        compute_extra_mems(extra_mems, dns, num_actions_allowed);
        //cout<<"compute_extra_mems()\n";


        int flag = run_mck_solver(extra_mems, budgets, picked_dn_indexes, total_budget, num_actions_allowed, s);

        int num_actual_inserts = 0;
        int num_actions = 0;

        int before_mck_attack = index.stats_.num_keys;

        int sum_1 = 0;

        for(int i = 0; i < picked_dn_indexes.size(); i++) {
            sum_1 += budgets[picked_dn_indexes[i]];
        }

        int nodeidx, num, budget, num_actions_chosen;

        for(int i = 0; i < picked_dn_indexes.size(); i++) {
            nodeidx = picked_dn_indexes[i] / num_actions_allowed;
            node = dns[nodeidx];

            //chosen_node_density.push_back(node->num_keys_ / static_cast<double> (node->data_capacity_));
            num_actions_chosen = picked_dn_indexes[i] % num_actions_allowed + 1;
            num_actions += num_actions_chosen;
            budget = budgets[picked_dn_indexes[i]];
            
            
            if(blackbox_enable == 1) {
                gen_poison_keys(substitude_poisoning_keys, index, node, budget * 1.05);
            }


            num = insert_keys_to_datanode(node, index, budget, 
                    num_actual_inserts, num_actions_chosen, poison_keys, 
                    substitude_poisoning_keys,
                    blackbox_enable);
            
            num_actual_inserts += num;
        }

        //bug here when switch catastrophic split to expand
        long long total_mem = index.model_size();
        total_mem += index.data_size();

        int num_expansions = 0;

        dns.clear();
        index.get_all_datanodes(dns);
        for(int i = 0; i < dns.size(); i++) {
            if(dns[i]->num_resizes_ != 0) {
                num_expansions += dns[i]->num_resizes_;
            }
        }

        s.mck_num_keys += num_actual_inserts;
        s.num_mck_call++;
        s.num_data_nodes = index.stats_.num_data_nodes;
        s.total_mem = total_mem;
        s.exp_num_actions = num_actions;

        s.aut_num_actions = index.stats_.num_sideways_splits + 
            index.stats_.num_downward_splits + index.stats_.num_expand_and_retrains + index.stats_.num_expand_and_scales
            + num_expansions;
        
        return num_actual_inserts;
    }




    void print_error(std::string s) 
    {
        std::cout<<"error occured at calling " + s <<std::endl;   
    }
    void print_execution_stats(std::string s)
    {
        std::cout<< "finish " + s +"...\n";
    }

    int get_datanode_idx(int index)
    {
        return 1; 
    }

    void to_string(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> nodes, 
            std::vector<int> extra_memory, 
            std::vector<int> budget,
            int num_actions_allowed,
            int num_nodes)
    {

        for(int i = 0; i < num_nodes; i++) {
            datanode_to_string(nodes[i]);
            std::cout<<"[";
            for(int j=1; j <= num_actions_allowed; j++) {
                std::cout<<"action=" << j <<" "
                    <<"memory="<<extra_memory[j-1]<<" "
                    <<"budget="<<budget[j-1]<<" ";
            }

            std::cout << "] }"<<std::endl;
        }
    }

    void datanode_to_string(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {
        std::cout<<"------NODE INFO------\n"
            <<"data_capacity_=" << node->data_capacity_<<"\n"
            <<"expansion_threshold_="<<node->expansion_threshold_<<"\n"
            <<"num_keys_="<<node->num_keys_<<"\n"
            <<"max_key_="<<node->max_key_<<"\n"
            <<"min_key_="<<node->min_key_<<"\n"
            <<"num_resizes_="<<node->num_resizes_<<"\n"
            <<"duplication_factor_="<<node->duplication_factor_<<"\n"
            <<"level_="<<node->level_<<"\n"
            <<"------END------"
            <<std::endl;
    }


    /*  num_nodes: this veriable used to keep track of how many datanodes 
     *                  generated after performing a action(split or expandsion).
     *  current_action: this veriable used to keep track of how many actions have performed
     *  actions_chosen: a string that store the actions we have performed
     *  res: all combinations will be store in this vector
     *
     * */
    void encode_actions(int num_actions,
            int num_nodes,
            std::string actions_chosen,
            int current_action,
            std::vector<std::string> &res)
    {

        if(num_actions < current_action) {
            res.push_back(actions_chosen);
            return;
        }

        if(current_action == 1) {
            encode_actions(num_actions, num_nodes + 1, "1s,", current_action + 1, res);
            encode_actions(num_actions, num_nodes, "1e,", current_action + 1, res);
        } else {
            for(int i = 0; i < num_nodes; i++) {
                std::string temp = turn_on_node_at_index(i, num_nodes);

                encode_actions(num_actions, 
                        num_nodes + 1, actions_chosen + temp + "s,", current_action + 1, res);

                encode_actions(num_actions, 
                        num_nodes, actions_chosen + temp + "e,", current_action + 1, res);
            }
        }
    }
    void get_parents(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> &parents,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index)
    {
        //print_execution_stats("get_parents() called");
        for(int i = 0; i < dns.size(); i++) {
            parents.push_back(index.get_parent(dns[i]));
        }
    }

    void get_parents(std::vector<std::pair<KEY_TYPE, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &dns,
            std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> &parents,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index)
    {
        //print_execution_stats("get_parents() called");
        for(int i = 0; i < dns.size(); i++) {
            parents.push_back(index.get_parent(dns[i].second));
        }
    }


    std::string turn_on_node_at_index(int target_index, int total_len)
    {
        std::string res="";
        for(int i = 0; i < total_len; i++) {
            if(target_index == i) {
                res.append("1");
            } else{
                res.append("0");
            }
        }
        return res;

    }
    void print_index_stats(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index)
    {
        std::cout<<"index stats={"<<"\n"
            <<"num_keys="<<index.stats_.num_keys<<"\n"
            <<" num_data_node="<<index.stats_.num_data_nodes<<"\n"
            <<" num_expand_and_scales="<<index.stats_.num_expand_and_scales<<"\n"
            <<" num_expand_and_retrains="<<index.stats_.num_expand_and_retrains<<"\n"
            <<" num_sideways_splits="<<index.stats_.num_sideways_splits<<"\n"
            <<" num_downward_splits="<<index.stats_.num_downward_splits<<"\n"
            <<"}"
            <<std::endl;
    }

    std::pair<KEY_TYPE, KEY_TYPE> get_datanode_key_boundaries(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE> *parent)
    {
        KEY_TYPE k = node->first_key();
        int bucketID = parent->model_.predict(k);
        bucketID = std::min<int>(std::max<int>(bucketID, 0), parent->num_children_ - 1);

        int repeats = 1 << node->duplication_factor_;
        int start_bucketID =
            bucketID - (bucketID % repeats);  // first bucket with same child
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with different child
        KEY_TYPE left_boundary_value =
            (start_bucketID - parent->model_.b_) / parent->model_.a_; //Rui: x = a*m + b 
        KEY_TYPE right_boundary_value =
            (end_bucketID - parent->model_.b_) / parent->model_.a_;

        return {left_boundary_value, right_boundary_value};
    }
    int get_partition_pos(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent,
            std::pair<KEY_TYPE,KEY_TYPE> *model = nullptr)
    {
        KEY_TYPE k = node->first_key();
        int bucketID = parent->model_.predict(k);
        bucketID = std::min<int>(std::max<int>(bucketID, 0), parent->num_children_ - 1);

        int repeats = 1 << node->duplication_factor_;
        int start_bucketID =
            bucketID - (bucketID % repeats);  // first bucket with same child
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with different child
        KEY_TYPE left_boundary_value =
            (start_bucketID - parent->model_.b_) / parent->model_.a_; //Rui: x = a*m + b 
        KEY_TYPE right_boundary_value =
            (end_bucketID - parent->model_.b_) / parent->model_.a_;

        alex::LinearModel<KEY_TYPE> base_model;
        base_model.a_ = 1.0 / (right_boundary_value - left_boundary_value);
        base_model.b_ = -1.0 * base_model.a_ * left_boundary_value;


        KEY_TYPE a = base_model.a_ * 2;
        KEY_TYPE b = base_model.b_ * 2;

        if(model != nullptr){
            model->first = a;
            model->second = b;
        }

        int right_boundary = node->lower_bound((1 - b) / a);


        return right_boundary;
    }


    bool should_insert_left_range(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            std::pair<KEY_TYPE, KEY_TYPE> &boundary,
            int counts[])
    {
        KEY_TYPE k = node->first_key();
        int bucketID = parent->model_.predict(k);
        bucketID = std::min<int>(std::max<int>(bucketID, 0), parent->num_children_ - 1);

        int repeats = 1 << node->duplication_factor_;
        int start_bucketID =
            bucketID - (bucketID % repeats);  // first bucket with same child
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with different child
        KEY_TYPE left_boundary_value =
            (start_bucketID - parent->model_.b_) / parent->model_.a_; //Rui: x = a*m + b 
        KEY_TYPE right_boundary_value =
            (end_bucketID - parent->model_.b_) / parent->model_.a_;

        boundary.first = left_boundary_value;
        boundary.second = right_boundary_value;


        KEY_TYPE mid_boundary_value = (left_boundary_value + right_boundary_value) / 2;


        int left_cnt, right_cnt, left_num_gap, right_num_gap;
        left_cnt = 0;
        right_cnt = 0;
        left_num_gap = 0;
        right_num_gap = 0;

        for(int i = 0; i < node->data_capacity_; i++) {
            if(node->check_exists(i)) {
                if(node->get_key(i) < mid_boundary_value){
                    left_cnt++;
                } else if(node->get_key(i) >= mid_boundary_value){
                    right_cnt++;
                }
            }
        }
        counts[0] = left_cnt;
        counts[1] = right_cnt;

        if(left_cnt > right_cnt){
            return true;    
        }

        return false;
    }

    bool should_trigger_action(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {
        if( (node->num_inserts_ + 1) % 64 == 0 && node->catastrophic_cost()){
            return true;
        }

        if(node->num_keys_ >= node->expansion_threshold_)
            return true;

        return false;
    }

    int insert(alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            KEY_TYPE keys[],
            int size, int node_id,
            ExpeData *l)
    {
        std::mt19937_64 gen_payload(std::random_device{}());

        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node = index.get_leaf(keys[0]);
        int num_actions = node->num_resizes_ + 
            index.stats_.num_expand_and_scales + 
            index.stats_.num_sideways_splits + 
            index.stats_.num_downward_splits +
            index.stats_.num_expand_and_retrains;

        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *prev_node, *next_node, *temp;

        prev_node = node->prev_leaf_;
        next_node = node->next_leaf_;

        int current_actions = 0;

        for(int i = 0; i < size; i++) {
            //std::cout<<"key="<<keys[i]<<" at i="<<i<<"\n"; 

            if(i == 100) {
                l->empirical_cost = node->empirical_cost();
                l->actual_cost = node->cost_;
            }

            int temp_expand_and_scale = node->num_resizes_ + index.stats_.num_expand_and_scales;
            int temp_sideway = index.stats_.num_sideways_splits;
            int temp_downward = index.stats_.num_downward_splits;
            int temp_expand_and_retrain = index.stats_.num_expand_and_retrains;

            node = index.get_leaf(keys[i]);
            int datanode_size = node->data_size();
            int node_resizes = node->num_resizes_;


            index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()), l);

            if(temp_sideway != index.stats_.num_sideways_splits || 
                    temp_downward != index.stats_.num_downward_splits ||
                    temp_expand_and_retrain != index.stats_.num_expand_and_retrains || node_resizes < node->num_resizes_) {
                std::cout<<"datanode_size="<<datanode_size<<"\n";
                l->action_keys += datanode_size;
            }

            int temp_num_expand_noretrain = 0;
            temp = prev_node->next_leaf_;

            while(temp != next_node) {
                temp_num_expand_noretrain += temp->num_resizes_;
                temp = temp->next_leaf_;
            }

            current_actions = temp_num_expand_noretrain + 
                index.stats_.num_expand_and_scales +
                index.stats_.num_sideways_splits + 
                index.stats_.num_downward_splits + 
                index.stats_.num_expand_and_retrains;


            //if(current_actions - node_id >= num_actions){
            //    return i + 1;
            //}
        }

        return size;
    }




    void collect_dn_stat_before_action(KEY_TYPE key, alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index, ExpeData *l)
    {
        std::mt19937_64 gen_payload(std::random_device{}());
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node = index.get_leaf(key);
        //alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent = index.get_parent(node);

        /*
           l->left_num_keys = get_left_partion_num_keys(node, parent);
           l->right_num_keys = get_right_partion_num_keys(node, parent);
           l->left_avg = get_avg_continuous_size(node, 0, boundary);
           l->right_avg = get_avg_continuous_size(node, boundary, node->data_capacity_);
           l->whole_avg = get_avg_continuous_size(node, 0, node->data_capacity_);
           l->empirical_cost = node->empirical_cost();
           l->segments.clear();
           l->gaps.clear();
           l->segments = get_node_segments(node);
           l->gaps = get_node_gaps(node);*/

        int temp_expand_and_scale = node->num_resizes_;
        int temp_sideway = index.stats_.num_sideways_splits;
        int temp_downward = index.stats_.num_downward_splits;
        int temp_expand_and_retrain = index.stats_.num_expand_and_retrains;

        index.insert(key, static_cast<PAYLOAD_TYPE>(gen_payload()), l);


        if(node && node->num_resizes_ != temp_expand_and_scale) {
            l->action_flag = 0;
        } else if(index.stats_.num_expand_and_retrains != temp_expand_and_retrain){
            l->action_flag = 0;
        } else if(index.stats_.num_sideways_splits != temp_sideway) {
            l->action_flag = 1;
        } else if(index.stats_.num_downward_splits != temp_downward) {
            l->action_flag = 1;
        }


        /*if(index.stats_.num_expand_and_scales!= temp_expand_and_scale || 
          index.stats_.num_expand_and_retrains != temp_expand_and_retrain) {
          l->action_flag = 0;
          } else if(index.stats_.num_sideways_splits != temp_sideway || 
          index.stats_.num_downward_splits != temp_downward) {
          l->action_flag = 1;
          }*/
    }

    /*
       void fill_gaps(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
       alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index)
       {
       KEY_TYPE k = node->first_key();
       int bucketID = parent->model_.predict(k);
       bucketID = std::min<int>(std::max<int>(bucketID, 0), parent->num_children_ - 1);
       alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent = index.get_parent(node);

       int repeats = 1 << node->duplication_factor_;
       int start_bucketID =
       bucketID - (bucketID % repeats);  // first bucket with same child
       int end_bucketID =
       start_bucketID + repeats;  // first bucket with different child
       KEY_TYPE left_boundary_value =
       (start_bucketID - parent->model_.b_) / parent->model_.a_; //Rui: x = a*m + b 
       KEY_TYPE right_boundary_value =
       (end_bucketID - parent->model_.b_) / parent->model_.a_;


       int start_pos = node->predict_position(left_boundary_value);
       int size = node->data_capacity_;
       }*/

    std::pair<KEY_TYPE,KEY_TYPE> find_segment_neighbor_boundary(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node){

        int size = node->data_capacity_;
        std::pair<int,int> max_pos;
        get_max_continuous_key_size(node, 0, size, max_pos);

        KEY_TYPE left_boundary, right_boundary;
        bool found = false;
        for(int i = max_pos.first - 1; i >= 0; i--) {
            if(node->check_exists(i)) {
                found = true;
                right_boundary = node->get_key(max_pos.first);
                left_boundary = node->get_key(i); 
                break; 
            }
        }

        if(found) {
            return {left_boundary, right_boundary}; 
        }


        for(int i = max_pos.second + 1; i < size; i++) {
            if(node->check_exists(i)) {
                found = true;
                left_boundary = node->get_key(max_pos.second); 
                right_boundary = node->get_key(i);
                break; 
            }
        }

        if(found) {
            return {left_boundary, right_boundary};
        }

        return {0.0, 0.0};

    }



    std::pair<KEY_TYPE,KEY_TYPE> get_avg_gap_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent)
    {
        KEY_TYPE left_avg_gap = 0.0;
        KEY_TYPE right_avg_gap = 0.0;
        int gaps = 0;
        int counts = 0;

        int size = node->data_capacity_;
        int boundary = get_partition_pos(node, parent);


        //left
        for(int i = 0; i < boundary; i++) {
            if(!node->check_exists(i)){
                gaps++;
            }
        }

        for(int i = 1; i < boundary; i++) {
            if(node->check_exists(i) && !node->check_exists(i-1)){
                counts++;
            }
        }

        left_avg_gap = gaps / counts;
        gaps = 0;
        counts = 0;
        //right
        for(int i = boundary; i < size; i++) {
            if(!node->check_exists(i)){
                gaps++;
            }
        }

        for(int i = boundary + 1; i < size; i++) {
            if(node->check_exists(i) && !node->check_exists(i-1)){
                counts++;
            }
        }

        right_avg_gap = gaps / counts;
        return {left_avg_gap, right_avg_gap};
    }


    KEY_TYPE get_avg_continuous_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start, int end)
    {
        int keys = 0;
        int num = 0;

        for(int i = start; i < end; i++) {
            if(node->check_exists(i)){
                keys++;
            }
        }

        for(int i = start + 1; i < end; i++) {
            if(!node->check_exists(i) && node->check_exists(i-1)){
                num++;
            }
        }

        if(node->check_exists(end - 1)) {
            num++;
        }

        return static_cast<KEY_TYPE> (keys) / num;
    }

    std::pair<int, int> get_max_gap(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent, 
            std::pair<int, int> &left_pos, std::pair<int, int> &right_pos)
    {
        int left_max_gap = -1;
        int right_max_gap = -1;
        int gaps = 0;
        int size = node->data_capacity_;
        int boundary = get_partition_pos(node, parent);

        //left
        for(int i = 0; i < boundary; i++) {
            if(!node->check_exists(i)){
                gaps++;
            } else {
                if(gaps > left_max_gap){
                    left_pos.second = i - 1;
                    left_pos.first = i  - gaps;
                    left_max_gap = gaps;
                }
                gaps = 0;
            }
        }

        gaps = 0;
        //right
        for(int i = boundary; i < size; i++) {
            if(!node->check_exists(i)){
                gaps++;
            }else {
                if(gaps > right_max_gap){
                    right_pos.first = i - gaps;
                    right_pos.second = i - 1;
                    right_max_gap = gaps;
                }
                gaps = 0;
            }
        }

        return {left_max_gap - 1, right_max_gap - 1};

    }

    std::vector<int> get_node_segments(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {
        std::vector<int> res;
        int count = 0;

        int end = node->data_capacity_;
        for(int i = 0; i < end; i++) {
            if(node->check_exists(i)){
                count++;
            } else {
                if(i >= 1 && node->check_exists(i - 1)){
                    res.push_back(count);
                }
                count = 0;
            }

        }

        if(node->check_exists(end - 1)) {
            res.push_back(count);
        }

        return res;

    }

    std::vector<int> get_node_gaps(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {
        std::vector<int> res;
        int count = 0; 
        int end = node->data_capacity_;

        for(int i = 0; i < end; i++) {
            if(!node->check_exists(i)){
                count++;
            } else {
                if(i >= 1 && !node->check_exists(i - 1)){
                    res.push_back(count);
                }
                count = 0;
            }
        }

        if(!node->check_exists(end - 1)) {
            res.push_back(count);
        }

        return res;
    }





    // 1st return value: a max length of consecutive key size
    // 2nd return value: start index
    std::pair<int,int> get_max_continuous_key(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end)
    {
        int cur_max = -1;
        int start_index = -1;
        int right_max_continuous = -1;
        int keys = 0;

        for(int i = start; i < end; i++) {
            if(node->check_exists(i)){
                keys++;
            } else {
                if(keys > cur_max){
                    start_index = i - keys;
                    cur_max = keys;
                }
                keys = 0;
            }
        }

        //this if statement checks the corner
        //If the largest sagment spans the last pos of array the above code wont update
        //start_index and cur_max
        if(keys != 0) {
            start_index = end - keys;
            cur_max = keys;
        }


        return {cur_max,start_index};
    }

    int get_max_continuous_key_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end,std::pair<int,int> &pos)
    {
        int cur_max = -1;
        int right_max_continuous = -1;
        int keys = 0;

        for(int i = start; i < end; i++) {
            if(node->check_exists(i)){
                keys++;
            } else {
                if(keys > cur_max){
                    pos.first = i - keys;
                    pos.second = i - 1;
                    cur_max = keys;
                }
                keys = 0;
            }
        }

        return cur_max;
    }

    int get_max_continuous_gap_size(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end)
    {
        int runing_size = 0;
        int max = 0; 

        for(int i = start; i < end; i++) {
            if(!node->check_exists(i)){
                runing_size++;
            } else {

                if(runing_size > max){
                    max = runing_size;     
                }
                runing_size = 0;
            }
        }

        return max;
    }

    int get_mid_occupied_pos(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            int start, int end) {
        
        int mid = (start + end) / 2;
        int left = mid;
        int right = mid;

        while(left > start && !node->check_exists(left) ) {
            left--;
        }



        while(right < end && !node->check_exists(right) ) {
            right++;
        }


        if(right == mid ) {
            return right;
        }

        if(left == mid) {
            return left;
        }

        if( node->check_exists(right) && std::abs(mid - left) > std::abs(mid - right)) {
            return right;
        }


        if(node->check_exists(left) && std::abs(mid - left) < std::abs(mid - right)) {
            return left;
        }



        return -1;
        
    }



    void get_poison_key_boundaries_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start_pos, int end_pos, std::pair<KEY_TYPE, KEY_TYPE> &poison_key_boundaries)
    {
        KEY_TYPE left_boundary_value = node->get_key(start_pos);
        KEY_TYPE right_boundary_value = node->get_key(end_pos);

        poison_key_boundaries.first = left_boundary_value;
        poison_key_boundaries.second = right_boundary_value;
    }

    std::pair<KEY_TYPE, KEY_TYPE> get_node_key_range(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, const 
            alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent)
    {

        if(!node || !parent) {
            return {-1, -1};
        }

        KEY_TYPE k = node->first_key();
        int bucketID = parent->model_.predict(k);

        if(parent->children_[bucketID] != node) {
            return {-1, -1};
        }

        bucketID = std::min<int>(std::max<int>(bucketID, 0), parent->num_children_ - 1);

        int repeats = 1 << node->duplication_factor_;
        int start_bucketID =
            bucketID - (bucketID % repeats);  // first bucket with same child
        int end_bucketID =
            start_bucketID + repeats;  // first bucket with different child
        KEY_TYPE left_boundary_value =
            (start_bucketID - parent->model_.b_) / parent->model_.a_; //Rui: x = a*m + b 
        KEY_TYPE right_boundary_value =
            (end_bucketID - parent->model_.b_) / parent->model_.a_;

        return {left_boundary_value, right_boundary_value};

    }


    void get_poison_key_boundaries_gaps(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            int start_pos, int end_pos, std::pair<KEY_TYPE, KEY_TYPE> &poison_key_boundaries)
    {
        if(start_pos == 0){
            poison_key_boundaries.first = node->get_key(end_pos + 1);
            return;
        }

        if(end_pos == 0){
            poison_key_boundaries.first = node->get_key(start_pos - 1);
            return;
        }


        poison_key_boundaries.first = node->get_key(start_pos - 1);
        poison_key_boundaries.second = node->get_key(end_pos + 1);
    }

    bool validate_keys(KEY_TYPE keys[], int size,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, 
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index)
    {
        for(int i = 0; i < size; i++) {
            if(index.get_leaf(keys[i]) != node){
                std::cout<<"key="<<keys[i]<<"\n";
                return false;
            }
        }
        return true;
    }

    int get_left_partion_num_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent){

        int boundary = get_partition_pos(node, parent);

        int cnt = 0;
        for (int i = 0; i < boundary; i++) {
            if(node->check_exists(i)) {
                cnt++;
            }
        } 
        return cnt;

    }

    int get_right_partion_num_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            const alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent){

        int boundary = get_partition_pos(node, parent);

        int size = node->data_capacity_;
        int cnt = 0;
        for (int i = boundary; i < size; i++) {
            if(node->check_exists(i)) {
                cnt++;
            }
        } 
        return cnt;
    }

    void dump_dn_stats(std::ofstream &dn_stat_out_file, 
            ExpeData *l, int i, int desired_key_size, int total_insert_keys,
            int data_capacity, std::pair<KEY_TYPE, KEY_TYPE> *model)
    {
        dn_stat_out_file<< i <<","
            <<l->left_avg<<","
            <<l->right_avg<<","
            <<l->whole_avg<<","
            <<desired_key_size<<","
            <<total_insert_keys<<","
            <<l->left_num_keys<<","
            <<l->right_num_keys<<","
            <<l->empirical_cost<<","
            <<l->action_flag<<","
            <<l->expand_expected_shifts_cost<<","
            <<l->expand_expected_search_iterations_cost<<","
            <<l->split_expected_shifts_cost<<","
            <<l->split_expected_search_iterations_cost<<","
            <<l->expand_expected_shifts<<","
            <<l->expand_expected_search_iterations<<","
            <<l->split_expected_shifts<<","
            <<l->split_expected_search_iterations<<","
            <<l->boundary<<","
            <<l->capacity
            <<std::endl;
    }

    void dump_keys(std::ofstream &out, std::vector<KEY_TYPE> &keys, int dn_idx)
    {
        int size = keys.size();
        for(int i = 0; i < size; i++) {
            out<<std::setprecision(std::numeric_limits<KEY_TYPE>::digits10)<<keys[i]<<",";
        }
    }

    void dump_segments(std::ofstream &out, std::vector<int> &segments, int dn_idx)
    {
        int size = segments.size();

        if(size == 0) {
            return;
        }

        for(int i = 0; i < size - 1; i++) {
            out<< segments[i] <<",";    
        } 
        out<<segments[size - 1]<<std::endl;
    }

    void dump_gaps(std::ofstream &out, std::vector<int> &gaps, int dn_idx)
    {
        int size = gaps.size();
        if(size == 0) {
            return;
        }


        for(int i = 0; i < size - 1; i++) {
            out<< gaps[i] <<",";    
        } 
        out<<gaps[size - 1]<<std::endl;
    }

    void merge_keys(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, ExpeData *l)
    {
        int size = node->data_capacity_;

        for(int i = 0; i < size; i++) {
            if(node->check_exists(i)) {
                l->keys.push_back(node->get_key(i));
            }
        }
    }

    void print_datanode_stats(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node){
        std::cout<<"datanode stat: {" <<"\n"
            <<" data_capacity_:"<<node->data_capacity_<<"\n"
            <<" action_threshold_:"<<node->expansion_threshold_<<"\n"
            <<" num_keys_::"<<node->num_keys_<<"\n"
            <<" expected_cost:"<<node->cost_<<"\n"
            <<" empirical_cost:"<<node->empirical_cost()<<"\n"
            <<" num_resizes_:"<<node->num_resizes_<<"\n"
            <<" shifts_per_insert:"<<node->shifts_per_insert()<<"\n"
            <<"}\n";
    }
    void get_nodes_with_max_gap_size(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &nodes, 
            std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> &get_nodes_with_max_gap_pairs)
    {    
        int size = nodes.size();
        for(int i = 0; i < size; i++) {
            get_nodes_with_max_gap_pairs.push_back({get_max_continuous_gap_size(nodes[i], 0, nodes[i]->data_capacity_), nodes[i]});
        }

    }


    int is_populated(KEY_TYPE injected_key, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {

        int predicted_pos = node->model_.predict(injected_key);

        if(!node->check_exists(predicted_pos)){
            return 0;
        }

        KEY_TYPE occupied_key = node->get_key(predicted_pos); 

        KEY_TYPE left_boundary = (predicted_pos - node->model_.b_) / node->model_.a_;
        KEY_TYPE right_boundary = (predicted_pos + 1 - node->model_.b_) / node->model_.a_;

        if(left_boundary <= occupied_key && occupied_key < right_boundary) {
            return 0;
        }

        return 1;
    }

    std::vector<int> get_pupolated_slots(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {
        std::vector<int> res;
        int size = node->data_capacity_;

        for(int i = 0; i < size; i++) {
            if(node->check_exists(i)) {
                std::pair<KEY_TYPE,KEY_TYPE> boundary = get_slot_key_boundary(node, i);
                KEY_TYPE k = node->get_key(i);

                if(k >= boundary.first && k < boundary.second) {
                    res.push_back(i);
                }
            }
        }
        return res;
    }

    std::vector<int> get_unpupolated_slots(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node)
    {


        int end = node->data_capacity_;
        std::vector<int> res;

        for(int i = 0; i < end; i++) {
            if(node->check_exists(i)) {
                std::pair<KEY_TYPE,KEY_TYPE> boundary = get_slot_key_boundary(node, i);
                KEY_TYPE k = node->get_key(i);

                if(k < boundary.first || k >= boundary.second) {
                    res.push_back(i);
                }
            } else {
                res.push_back(i);
            }
        }
        return res;
    }

    std::pair<KEY_TYPE, KEY_TYPE> get_slot_key_boundary(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node, int pos)
    {
        std::pair<KEY_TYPE, KEY_TYPE> res;

        res.first = (pos - node->model_.b_) / node->model_.a_;
        res.second = (pos + 1 - node->model_.b_) / node->model_.a_;

        return res;

    }


    void gen_keys_over_slots_evenly(alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            std::vector<int> slots, 
            int keys_per_slot,
            KEY_TYPE poison_keys[], 
            int desired_key_size)
    {
        std::pair<KEY_TYPE, KEY_TYPE> slot_boundary;
        std::pair<KEY_TYPE, KEY_TYPE> node_key_boundary;
        int cnt = 0;
        node_key_boundary = get_datanode_key_boundaries(node, index.get_parent(node));

        for(int i = 0; i < slots.size(); i++) {
            slot_boundary = get_slot_key_boundary(node, slots[i]);
            if(slot_boundary.first < node_key_boundary.first || slot_boundary.second >= node_key_boundary.second) {
                continue;
            }
            cnt += keys_per_slot;
        }

        int i = slots.size() / 2;
        while(cnt < desired_key_size) {

            slot_boundary = get_slot_key_boundary(node, slots[i]);
            if(slot_boundary.first < node_key_boundary.first || slot_boundary.second >= node_key_boundary.second) {
                i++;
                continue;
            }
            break;
        } 
    }

    void dump_keys(std::vector<KEY_TYPE> &poison_keys,
            std::vector<KEY_TYPE> &substitude_poisoning_keys)
    {
        for(int i = 0 ; i < poison_keys.size(); i++) {
            substitude_poisoning_keys.push_back(poison_keys[i]);
        }
    }

    void insert_left_keys(std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> &dns,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int num_keys,
            std::vector<KEY_TYPE> &poison_keys,
            std::vector<double> &chosen_node_density,
            std::vector<KEY_TYPE> &substitude_poisoning_keys)
    {
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node;


        dns.clear();
        poison_keys.clear();
        index.get_all_datanodes(dns);

        std::sort(dns.begin(), 
                dns.begin() + dns.size(),
                [](auto const* a, auto const* b) { 
                    return a->data_size() > b->data_size();
                });


        node = dns[0];
        gen_poison_keys(poison_keys, index, node, num_keys);
        insert_poison_keys(poison_keys, index, poison_keys.size());

    }

    int launch_greedy_attack(
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int total_budget, 
            int num_actions_allowed,
            int key_distribution_flag,
            Stats &s,
            std::vector<double> &chosen_node_density)
    {

        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
        std::vector<int> selected_dn_indexes;
        std::vector<int64_t> budgets;
        std::vector<KEY_TYPE> poison_keys;

        dns.clear();
        selected_dn_indexes.clear();
        budgets.clear();

        index.get_all_datanodes(dns);

        std::sort(dns.begin(), 
                dns.begin() + dns.size(),
                [](auto const* a, auto const* b) { 
                return a->data_size() / (round_up(a->expansion_threshold_) - a->num_keys_+1) 
                > b->data_size() / (round_up(b->expansion_threshold_) - b->num_keys_+1);});

        compute_budgets(budgets, dns, num_actions_allowed);


        int expense = total_budget;

        for(int i = 0; i <budgets.size();i++) {

            if(expense <= 0) {
                break;
            }
            selected_dn_indexes.push_back(i);
            expense -= budgets[i];
            chosen_node_density.push_back(dns[i]->num_keys_ / static_cast<double> (dns[i]->data_capacity_) );
        }


        gen_poison_keys(poison_keys, selected_dn_indexes, 
                budgets, dns, index, num_actions_allowed, key_distribution_flag, total_budget);

    }



    int launch_greedy_by_data_size_attack(
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int total_budget, 
            int num_actions_allowed,
            int key_distribution_flag,
            Stats &s)
    {
        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
        std::vector<int> selected_dn_indexes;
        std::vector<int64_t> budgets;
        std::vector<KEY_TYPE> poison_keys;

        dns.clear();
        selected_dn_indexes.clear();
        budgets.clear();

        index.get_all_datanodes(dns);

        std::sort(dns.begin(), 
                dns.begin() + dns.size(),
                [](auto const* a, auto const* b) { return a->data_size() > b->data_size(); });

        compute_budgets(budgets, dns, num_actions_allowed);


        int expense = total_budget;

        for(int i = 0; i <budgets.size();i++) {

            if(expense <= 0) {
                break;
            }
            selected_dn_indexes.push_back(i);
            expense -= budgets[i];
        }

        gen_poison_keys(poison_keys, selected_dn_indexes, 
                budgets, dns, index, num_actions_allowed, key_distribution_flag, total_budget);

        int ins = insert_poison_keys(poison_keys, index, poison_keys.size());


        return ins;

    }

    void gen_keys(std::vector<KEY_TYPE> &poison_keys,
            alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node,
            alex::Alex<KEY_TYPE, PAYLOAD_TYPE> &index,
            int budget)
    {

        std::random_device rd;
        std::default_random_engine eng(rd());
        std::pair<KEY_TYPE, KEY_TYPE> node_key_boundary;

        if(std::is_floating_point<KEY_TYPE>::value){
            //double
            node_key_boundary = get_datanode_key_boundaries(node, 
                    index.get_parent(node));

            std::uniform_real_distribution<double> uniform_distr(node_key_boundary.first,
                    node_key_boundary.second);

            for(int i = 0; i< budget; i++) {
                poison_keys.push_back(uniform_distr(eng));
            }
        }else {
            //int64_t

            node_key_boundary = get_datanode_key_boundaries(node, 
                    index.get_parent(node));


            std::uniform_int_distribution<int> uniform_distr(node_key_boundary.first,
                    node_key_boundary.second);

            for(int i = 0; i< budget; i++) {
                poison_keys.push_back(uniform_distr(eng));
            }

        }
    }
}
