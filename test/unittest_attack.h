#include "doctest.h"
#include "attack.h"
#include "alex.h"
#include "../src/benchmark/utils.h"
#include "../src/core/alex_base.h"
#include "../src/core/leap_util.h"
#include <iterator>
#include <experimental/algorithm>
#include <random>


#define KEY_TYPE double
#define PAYLOAD_TYPE double


TEST_SUITE("attack") {




    /*

    TEST_CASE("Test1") {

        std::cout<<"======Test1 starts======\n";
        int total_num_keys = 20000000;
        int init_num_keys = total_num_keys / 2;

        auto keys = new KEY_TYPE[total_num_keys];
        auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
        std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
        std::vector<std::pair<double, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> node_cost_diff_pairs;
        std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> node_max_gap_pairs;

        std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

        attack::print_execution_stats("Initial");

        load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
        attack::print_execution_stats("Finishi load_binary_data()");

        std::mt19937_64 gen_payload(std::random_device{}());
        for (int i = 0; i < init_num_keys; i++) {
            values[i].first = keys[i];
            values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
        }


        std::ofstream outputFile;
        //outputFile.open("/root/data/ALEX/attack/pupolated_slots.txt");
        ExpeData l;
        alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE> *node, *prev_node, *next_node, *temp;
        std::vector<int> pup_slots, unpup_slots;
        int desired_key_size;
        std::vector<int> actions = {20};
        int action_size = 4;
        int error_node = 0;
        int size = 10;
        int num_pure_expends = 0.3 * size;
        int expected_expands = 30;
        int sample_size = 1;

        std::vector<int> sampled_node_index;

        std::random_device rd;
        std::default_random_engine eng(rd());
        std::uniform_int_distribution<int> distr(0, 8000);

        for(int i = 0; i < sample_size; i++) 
        {
            sampled_node_index.push_back(distr(eng));
        }


        for(int n = 0; n < sample_size; n++) {

            for(int i = 0; i < actions.size(); i++) {
                outputFile.open("/root/data/ALEX/attack/nodeid_"+
                        std::to_string(n) + 
                        "_actions_" + std::to_string(actions[i])+ ".txt");

                for(int k = 0; k < size; k++) {
                    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
                    std::sort(values, values + init_num_keys,
                            [](auto const& a, auto const& b) { return a.first < b.first; });

                    dns.clear();
                    node_cost_diff_pairs.clear();
                    node_max_gap_pairs.clear();
                    parents.clear();


                    index.bulk_load(values, init_num_keys);
                    attack::print_execution_stats("Finish bulk_load() insertions");
                    attack::print_index_stats(index);

                    for(int i = init_num_keys; i < total_num_keys; i++) {
                        index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
                    }

                    attack::print_execution_stats("Finish insert() insertions");
                    //attack::print_index_stats(index);

                    index.get_all_datanodes(dns);

                    //attack::get_nodes_with_max_gap_size(dns, node_max_gap_pairs);

                    //std::sort(node_max_gap_pairs.begin(), node_max_gap_pairs.begin() + node_max_gap_pairs.size(),
                    //        [](auto const& a, auto const& b) { return a.first > b.first; });

                    //attack::print_execution_stats("Finish sorting node_cost_diff_pairs vector");
                    //attack::get_parents(dns, parents, index);

                    node = dns[sampled_node_index[n]];
                    prev_node = node->prev_leaf_;
                    next_node = node->next_leaf_;

                    int memory = 0;
                    int memory_before_actions = node->data_size();
                    int num_expand_noretrain_before_insert = node->num_resizes_ + index.stats_.num_expand_and_scales;
                    int num_expand_retrain_before_insert = index.stats_.num_expand_and_retrains;
                    int num_splits_before_insert = index.stats_.num_downward_splits + index.stats_.num_sideways_splits;
                    int num_keys_before_insert = index.stats_.num_keys;


                    for(int j = 0; j < expected_expands; j++) {
                        if(j == 0) {
                            desired_key_size = attack::get_num_keys(node->expansion_threshold_);
                        } else{
                            desired_key_size = attack::get_num_keys(desired_key_size / 0.6 * 0.8);
                        }
                    }

                    desired_key_size = desired_key_size - node->num_keys_ + 1;

                    auto keys_1 = new KEY_TYPE[desired_key_size];
                    int keys_per_slot = 0;
                    std::vector<int> slots;
                    slots.clear();

                    pup_slots = attack::get_unpupolated_slots(node);
                    unpup_slots = attack::get_pupolated_slots(node);

                    int combin_sample_lenth = (unpup_slots.size() + pup_slots.size()) / 2;
                    int combin_mid = combin_sample_lenth / 2;
                    int unpup_sample_lenth = unpup_slots.size() / 2;
                    int pup_sample_lenth = pup_slots.size() / 2;


                    if(k == 0) {
                        std::vector<int> new_slots  = pup_slots;
                        keys_per_slot = desired_key_size / new_slots.size();

                        attack::gen_keys_over_slots_evenly(
                                node,
                                index,
                                new_slots,
                                keys_per_slot,
                                keys_1, 
                                desired_key_size); 
                    } else if(k == 1) {
                        std::vector<int> new_slots  = unpup_slots;
                        keys_per_slot = desired_key_size / new_slots.size();

                        attack::gen_keys_over_slots_evenly(
                                node,
                                index,
                                new_slots,
                                keys_per_slot,
                                keys_1, 
                                desired_key_size); 

                    } else if(k > 1 && k < size - num_pure_expends){
                        std::vector<int> new_slots;
                        new_slots.clear();
                        if(k % 3 == 0) {
                            std::experimental::sample(unpup_slots.begin(), unpup_slots.end(), std::back_inserter(new_slots), 
                                    unpup_sample_lenth, std::mt19937{std::random_device{}()});

                        } else if( k % 3 == 1){
                            std::experimental::sample(pup_slots.begin(), pup_slots.end(), std::back_inserter(new_slots), 
                                    pup_sample_lenth, std::mt19937{std::random_device{}()});

                        } else {
                            std::experimental::sample(pup_slots.begin(), pup_slots.end(), std::back_inserter(new_slots), 
                                    combin_mid, std::mt19937{std::random_device{}()});
                            std::experimental::sample(unpup_slots.begin(), unpup_slots.end(), std::back_inserter(new_slots), 
                                    combin_sample_lenth - combin_mid, std::mt19937{std::random_device{}()});
                        }

                        keys_per_slot = desired_key_size / new_slots.size();
                        attack::gen_keys_over_slots_evenly(
                                node,
                                index,
                                new_slots,
                                keys_per_slot,
                                keys_1, 
                                desired_key_size); 


                    } else{

                        attack::gen_keys_by_boundary(
                                node,
                                index,
                                keys_1, 
                                desired_key_size,
                                k); 
                    }

                    if(!attack::validate_keys(keys_1, desired_key_size, node, index)){
                        std::cout<<"invalidate key found at k=" << k <<"\n";
                        continue;
                    }


                    resetExpeData(&l);
                    attack::print_index_stats(index);
                    int actual_num_key_inserted = attack::insert(index,keys_1, desired_key_size, actions[i], &l);
                    attack::print_index_stats(index);

                    temp = prev_node->next_leaf_;
                    int total_mem = 0;
                    int temp_num_expand_noretrain = 0;

                    int ta = 0;

                    while(temp != next_node) {
                        total_mem += temp->data_size();
                        temp_num_expand_noretrain += temp->num_resizes_;
                        temp = temp->next_leaf_;
                        ta++;
                    }

                    std::cout<<"=========end=======\n";
                    std::cout<<"t="<<ta<<"\n"
                        <<"temp_num_expand_noretrain="<<temp_num_expand_noretrain<<"\n"
                        <<"num_expand_noretrain_before_insert="<<num_expand_noretrain_before_insert<<"\n";


                    outputFile
                        <<k<<","
                        <<index.stats_.num_keys - num_keys_before_insert<<","
                        <<memory_before_actions<<","
                        <<total_mem<<","
                        <<l.action_keys<<","
                        <<index.stats_.num_expand_and_scales + temp_num_expand_noretrain - num_expand_noretrain_before_insert<<","
                        <<index.stats_.num_expand_and_retrains - num_expand_retrain_before_insert<<","
                        <<index.stats_.num_downward_splits + index.stats_.num_sideways_splits - num_splits_before_insert
                        <<"\n";

                    delete[] keys_1;
                }

                outputFile.close();
            }
        }
        delete[] keys;
        delete[] values;
    }*/

    /*


       TEST_CASE("Test2") {

       std::cout<<"======Test2======\n";
       int total_num_keys = 20000000;
       int init_num_keys = total_num_keys / 2;

       auto keys = new KEY_TYPE[total_num_keys];
       auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
       std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
       std::vector<std::pair<double, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> node_cost_diff_pairs;
       std::vector<std::pair<int, alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*>> node_max_gap_pairs;
       std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

       attack::print_execution_stats("Initial");

       load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
       attack::print_execution_stats("Finishi load_binary_data()");

       std::mt19937_64 gen_payload(std::random_device{}());
       for (int i = 0; i < init_num_keys; i++) {
       values[i].first = keys[i];
       values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
       }

       attack::print_execution_stats("key-value pair generation");

       std::sort(values, values + init_num_keys,
       [](auto const& a, auto const& b) { return a.first < b.first; });

    //outputFile.open("/root/data/ALEX/attack/longitudes_single_point_inject_largest_num.txt");

    attack::get_nodes_with_max_gap_size(dns, node_max_gap_pairs);


    for(int i = 0; i < 3; i++) {

    //clean up   
    node_max_gap_pairs.clear();
    dns.clear();

    //init ALEX and insert 200M keys
    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    index.bulk_load(values, init_num_keys);

    for(int i = init_num_keys; i < total_num_keys; i++) {
    index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    index.get_all_datanodes(dns);
    attack::get_nodes_with_max_gap_size(dns, node_max_gap_pairs);

    std::sort(node_max_gap_pairs.begin(), node_max_gap_pairs.begin() + node_max_gap_pairs.size(),
    [](auto const& a, auto const& b) { return a.first > b.first; });

    alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>* node = node_max_gap_pairs[0].second;
    alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>* parent = index.get_parent(node);


    int desired_key_size = 
    attack::get_num_keys(node->expansion_threshold_) - 
    node->num_keys_ + 1;

    auto poison_keys =  new KEY_TYPE[desired_key_size];

    //get key boundaries 
    ExpeData l;

    if(i == 0 ){
    std::vector<int> pupolated_slots = attack::get_pupolated_slots(node);
    int selected_slot = pupolated_slots[pupolated_slots.size() / 2];
    std::pair<double,double> boundary = attack::get_slot_key_boundary(node, selected_slot);
    attack::gen_key(boundary.first, boundary.second, poison_keys, 0, desired_key_size);  
} else if( i == 1) {

    std::vector<int> pupolated_slots = attack::get_pupolated_slots(node);

    attack::gen_keys_over_pupolated_slots_uniformlly(
            node, 
            pupolated_slots, 
            poison_keys, 
            desired_key_size);
} else {
    std::vector<int> unpupolated_slots = attack::get_unpupolated_slots(node);
    attack::gen_keys_over_unpupolated_slots_uniformlly(
            node,
            unpupolated_slots, 
            poison_keys, 
            desired_key_size); 
}


if (!attack::validate_keys(poison_keys, desired_key_size, node, index)) {
    std::cout<<"invalidate key found at i=" << i <<"\n";
    continue;
}



attack::print_index_stats(index);
attack::print_datanode_stats(node);
int actual_num_key_inserted = attack::insert(index,poison_keys, desired_key_size, 0, &l);

attack::print_index_stats(index);



std::cout<<"stats={" <<"\n"
<<" i=" <<i<<"\n"
<<" actual_num_key_inserted="<<actual_num_key_inserted<<"\n"
<<" desired_key_size="<<desired_key_size<<"\n"
<<" trigger fanout_tree=" << l.action_flag<<"\n"
<<" empirical_iters="<<l.empirical_search_iters<<"\n"
<<" empirical_shifts="<<l.empirical_shifts<<"\n"
<<" combined_iters="<<l.expected_combin_search_iters<<"\n"
<<" combined_shifts="<<l.expected_combin_shifts<<"\n"
<<" left_search_iters="<<l.expected_left_search_iters<<"\n"
<<" left_shifts="<<l.expected_left_shifts<<"\n"
<<" right_search_iters="<<l.expected_right_search_iters<<"\n"
<<" right_shifts="<<l.expected_right_shifts<<"\n"
<<"}\n";



outputFile<<k<<","
<<l.action_flag<<","
<<l.empirical_search_iters<<","
<<l.empirical_shifts<<","
<<populated<<","
<<l.expected_combin_search_iters<<","
<<l.expected_combin_shifts<<","
<<l.expected_combin_pos_occup<<","
<<l.expected_left_search_iters<<","
<<l.expected_left_shifts<<","
<<l.expected_left_pos_occup<<","
<<l.expected_right_search_iters<<","
<<l.expected_right_shifts<<","
<<l.expected_right_pos_occup<<","
<<l.expected_seperated_insert_flag
<<"\n";

std::cout<<"=================\n";
//clean up index
} 

//outputFile.close();

delete[] keys;
delete[] values;
}



TEST_CASE("ML_Model"){


    alex::LinearModel<double> model_1;
    alex::LinearModelBuilder<double> builder_1(&model_1);


    int desired_key_size = 10;
    auto poison_keys = new KEY_TYPE[desired_key_size];

    attack::gen_key(1, 2, poison_keys, 0, desired_key_size);

    for(int i = 0; i < desired_key_size; i++) {
        builder_1.add(poison_keys[i], i);
        std::cout<<"add point(" << poison_keys[i] <<"," << i << ")\n";
    }

    builder_1.build();




    alex::LinearModel<double> model_2;
    alex::LinearModelBuilder<double> builder_2(&model_2);

    desired_key_size = 20;

    for(int i = 0; i < desired_key_size; i++) {
        if(i == 1) {
            for(int j = 0; j < 11; j++) {
                builder_2.add(poison_keys[i], j + i);
                std::cout<<"add point(" << poison_keys[i] <<"," << j + i<< ")\n";
            }
            i += 10;

        } else {
            builder_2.add(poison_keys[i % 10], i);
            std::cout<<"add point(" << poison_keys[i % 10] <<"," << i << ")\n";
        }

    }

    builder_2.build();

    if(model_1.a_ != model_2.a_) {
        std::cout<< "model_1.a_="<<model_1.a_<<"\n"
            <<"model_2.a_="<<model_2.a_<<"\n";
    }

    if(model_1.b_ != model_2.b_) {
        std::cout<< "model_1.b_="<<model_1.b_<<"\n"
            <<"model_2.b_="<<model_2.b_<<"\n";
    }


}





TEST_CASE("RandomTest") {

    std::cout<<"======RandomTest Starts======\n";
    int total_num_keys = 20000000;
    int init_num_keys = total_num_keys / 2;

    auto keys = new KEY_TYPE[total_num_keys];
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

    attack::print_execution_stats("Initial");

    load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
    attack::print_execution_stats("Finishi load_binary_data()");

    std::mt19937_64 gen_payload(std::random_device{}());
    for (int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }


    attack::print_execution_stats("key-value pair generation");

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index.bulk_load(values, init_num_keys);
    attack::print_execution_stats("Finish bulk_load() insertions");
    attack::print_index_stats(index);

    std::ofstream outputFile, outputFile1, outputFile2, outputFile3, outputFile4;
    outputFile.open("/root/data/ALEX/attack/longitudes_dn_stats_100m.txt");
    //outputFile1.open("/root/data/ALEX/attack/split_segments_distr_100m.txt");
    //outputFile2.open("/root/data/ALEX/attack/split_gaps_distr_100m.txt");

    //outputFile3.open("/root/data/ALEX/attack/expand_segments_distr_100m.txt");
    //outputFile4.open("/root/data/ALEX/attack/expand_gaps_distr_100m.txt");


    ExpeData l;

    for(int i = init_num_keys; i < total_num_keys; i++) {
        l.action_flag = -1;
        attack::collect_dn_stat_before_action(keys[i], index, &l);

        //std::cout<<"l.action_flag="<<l.action_flag<<"\n";
        if(l.action_flag == 1 || l.action_flag == 0) {
            //std::cout<<"l.action_flag="<<l.action_flag<<"\n";
            attack::dump_dn_stats(outputFile, &l, i, 0, 0, 0, nullptr);
        }
    }
    attack::print_index_stats(index);
    std::cout<<"======RandomTest Ends======\n";

    outputFile.close();
    outputFile1.close();
    outputFile2.close();
    outputFile3.close();
    outputFile4.close();

    delete[] keys;
    delete[] values;
}




TEST_CASE("Test2") {

    std::cout<<"======targeting consective key half-half starts======\n";
    int total_num_keys = 20000000;
    int init_num_keys = total_num_keys / 2;

    auto keys = new KEY_TYPE[total_num_keys];
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

    attack::print_execution_stats("Initial");

    load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
    attack::print_execution_stats("Finishi load_binary_data()");

    std::mt19937_64 gen_payload(std::random_device{}());
    for (int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }


    attack::print_execution_stats("key-value pair generation");

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index.bulk_load(values, init_num_keys);
    attack::print_execution_stats("Finish bulk_load() insertions");
    attack::print_index_stats(index);


    for(int i = init_num_keys; i < total_num_keys; i++) {
        index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    attack::print_execution_stats("Finish insert() insertions");
    attack::print_index_stats(index);

    index.get_all_datanodes(dns);
    attack::get_parents(dns, parents, index);

    std::pair<int, int> con_left_pos, con_right_pos, gap_left_pos, gap_right_pos;
    std::pair<int, int> max_left_right_keys, max_left_right_gaps;
    std::pair<double, double> poison_key_boundaries;


    attack::print_index_stats(index);
    std::ofstream outputFile, outputFile1;
    outputFile.open("/root/data/ALEX/attack/longitudes_stats_max_segment_half_half.txt");

    int start_pos = 0;
    int end_pos = 0;
    int error_node = 0;
    double avg_con_key = 0.0;

    double left_avg, right_avg, whole_avg;
    int left_max, right_max, desired_key_size, boundary;
    std::pair<int,int> left_max_pos, right_max_pos;


    for(int i= 0; i < dns.size(); i++) {
        boundary = attack::get_partition_pos(dns[i], parents[i]);
        left_max = attack::get_max_continuous_key_size(dns[i], 0, boundary, left_max_pos);
        right_max = attack::get_max_continuous_key_size(dns[i], boundary, dns[i]->data_capacity_, right_max_pos);

        desired_key_size = attack::get_num_keys(dns[i]->expansion_threshold_) - dns[i]->num_keys_ + 1;
        auto keys_1 = new KEY_TYPE[desired_key_size];

        int left_num = desired_key_size / 2;
        //left
        attack::get_poison_key_boundaries_keys(dns[i], left_max_pos.first, 
                left_max_pos.second, poison_key_boundaries);

        attack::gen_key(poison_key_boundaries.first, poison_key_boundaries.second, 
                keys_1, 0, left_num);
        //right
        attack::get_poison_key_boundaries_keys(dns[i], right_max_pos.first, 
                right_max_pos.second, poison_key_boundaries);
        attack::gen_key(poison_key_boundaries.first, poison_key_boundaries.second, 
                keys_1, left_num, desired_key_size);


        std::sort(keys_1, keys_1 + desired_key_size,
                [](auto const& a, auto const& b) { return a < b; });

        if(!attack::validate_keys(keys_1, desired_key_size, dns[i],index)){
            outputFile1<<i
                <<std::endl;
            error_node++;
            continue;
        }

        int total_insert_keys = dns[i]->num_keys_;
        ExpeData l;

        attack::insert(index,keys_1, desired_key_size, &l);

        outputFile<< i <<","
            <<l.left_avg<<","
            <<l.right_avg<<","
            <<l.whole_avg<<","
            <<desired_key_size<<","
            <<total_insert_keys<<","
            <<l.left_num_keys<<","
            <<l.right_num_keys<<","
            <<l.node_intra_node_cost_1<<","
            <<l.node_traversal_cost_1<<","
            <<l.node_intra_node_cost_2_left<<","
            <<l.node_intra_node_cost_2_right<<","
            <<l.node_traversal_cost_2<<","
            <<l.empirical_cost<<","
            <<l.action_flag
            <<std::endl;
        delete[] keys_1;
    }

    attack::print_index_stats(index);
    std::cout<<"total error_node="<<error_node<<"\n";
    std::cout<<"======targeting consective key half-half ends======\n";

    outputFile.close();

    delete[] keys;
    delete[] values;
}
TEST_CASE("Test3") {

    std::cout<<"======targeting connect segment&neighbor starts======\n";
    int total_num_keys = 20000000;
    int init_num_keys = total_num_keys / 2;

    auto keys = new KEY_TYPE[total_num_keys];
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

    attack::print_execution_stats("Initial");

    load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
    attack::print_execution_stats("Finishi load_binary_data()");

    std::mt19937_64 gen_payload(std::random_device{}());
    for (int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }


    attack::print_execution_stats("key-value pair generation");

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index.bulk_load(values, init_num_keys);
    attack::print_execution_stats("Finish bulk_load() insertions");
    attack::print_index_stats(index);


    for(int i = init_num_keys; i < total_num_keys; i++) {
        index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    attack::print_execution_stats("Finish insert() insertions");
    attack::print_index_stats(index);

    index.get_all_datanodes(dns);
    attack::get_parents(dns, parents, index);

    std::pair<int, int> con_left_pos, con_right_pos, gap_left_pos, gap_right_pos;
    std::pair<int, int> max_left_right_keys, max_left_right_gaps;
    std::pair<double, double> poison_key_boundaries;


    attack::print_index_stats(index);
    std::ofstream outputFile, outputFile1;
    outputFile.open("/root/data/ALEX/attack/longitudes_stats_fill_gap_segment_neighbor.txt");

    int start_pos = 0;
    int end_pos = 0;
    int error_node = 0;
    double avg_con_key = 0.0;

    double left_avg, right_avg, whole_avg;
    int left_max, right_max, desired_key_size, boundary;
    std::pair<int,int> left_max_pos, right_max_pos;


    for(int i= 0; i < dns.size(); i++) {

        poison_key_boundaries = attack::find_segment_neighbor_boundary(dns[i]);

        desired_key_size = attack::get_num_keys(dns[i]->expansion_threshold_) - dns[i]->num_keys_ + 1;
        auto keys_1 = new KEY_TYPE[desired_key_size];

        attack::gen_key(poison_key_boundaries.first, poison_key_boundaries.second, 
                keys_1, 0, desired_key_size);

        std::sort(keys_1, keys_1 + desired_key_size,
                [](auto const& a, auto const& b) { return a < b; });

        if(!attack::validate_keys(keys_1, desired_key_size, dns[i],index)){
            outputFile1<<i
                <<std::endl;
            error_node++;
            continue;
        }

        int total_insert_keys = dns[i]->num_keys_;
        ExpeData l;

        attack::insert(index,keys_1, desired_key_size, &l);

        outputFile<< i <<","
            <<l.left_avg<<","
            <<l.right_avg<<","
            <<l.whole_avg<<","
            <<desired_key_size<<","
            <<total_insert_keys<<","
            <<l.left_num_keys<<","
            <<l.right_num_keys<<","
            <<l.node_intra_node_cost_1<<","
            <<l.node_traversal_cost_1<<","
            <<l.node_intra_node_cost_2_left<<","
            <<l.node_intra_node_cost_2_right<<","
            <<l.node_traversal_cost_2<<","
            <<l.empirical_cost<<","
            <<l.action_flag
            <<std::endl;
        delete[] keys_1;
    }

    attack::print_index_stats(index);
    std::cout<<"total error_node="<<error_node<<"\n";
    std::cout<<"======targeting increase segment number ends======\n";

    outputFile.close();

    delete[] keys;
    delete[] values;
}
*/


/*
   TEST_CASE("Test2") {

   std::cout<<"======targeting relative small max consective key starts======\n";
   int total_num_keys = 20000000;
   int init_num_keys = total_num_keys / 2;

   auto keys = new KEY_TYPE[total_num_keys];
   auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
   std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
   std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

   attack::print_execution_stats("Initial");

   load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
   attack::print_execution_stats("Finishi load_binary_data()");

   std::mt19937_64 gen_payload(std::random_device{}());
   for (int i = 0; i < init_num_keys; i++) {
   values[i].first = keys[i];
   values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
   }


   attack::print_execution_stats("key-value pair generation");

   alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
   std::sort(values, values + init_num_keys,
   [](auto const& a, auto const& b) { return a.first < b.first; });

   index.bulk_load(values, init_num_keys);
   attack::print_execution_stats("Finish bulk_load() insertions");
   attack::print_index_stats(index);

   for(int i = init_num_keys; i < total_num_keys; i++) {
   index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
   }

   attack::print_execution_stats("Finish insert() insertions");
   attack::print_index_stats(index);

   index.get_all_datanodes(dns);
   attack::get_parents(dns, parents, index);

   std::pair<int, int> con_left_pos, con_right_pos, gap_left_pos, gap_right_pos;
   std::pair<int, int> max_left_right_keys, max_left_right_gaps;
   std::pair<double, double> poison_key_boundaries;


   attack::print_index_stats(index);
   std::ofstream outputFile, outputFile1;
   outputFile.open("/root/data/ALEX/attack/max_consecutive_key_stats_longitudes_small.txt");
   outputFile1.open("/root/data/ALEX/attack/longitudes_small_error_node_logs.txt");

   int start_pos = 0;
   int end_pos = 0;
   int error_node = 0;
   double avg_con_key = 0.0;

   double left_avg, right_avg, whole_avg;
   int left_max, right_max, desired_key_size, boundary;
   std::pair<int,int> left_max_pos, right_max_pos;


   for(int i= 0; i < dns.size(); i++) {
   boundary = attack::get_partition_pos(dns[i], parents[i]);
   left_avg = attack::get_avg_continuous_size(dns[i], 0, boundary);
   right_avg = attack::get_avg_continuous_size(dns[i], boundary, dns[i]->data_capacity_);
   whole_avg = attack::get_avg_continuous_size(dns[i], 0, dns[i]->data_capacity_);

   left_max = attack::get_max_continuous_key_size(dns[i], 0, boundary, left_max_pos);
right_max = attack::get_max_continuous_key_size(dns[i], boundary, dns[i]->data_capacity_, right_max_pos);

desired_key_size = attack::get_num_keys(dns[i]->expansion_threshold_) - dns[i]->num_keys_ + 1;
auto keys_1 = new KEY_TYPE[desired_key_size];

if(left_max <= right_max) {
    attack::get_poison_key_boundaries_keys(dns[i],left_max_pos.first, 
            left_max_pos.second,
            poison_key_boundaries);
    start_pos = left_max_pos.first;
    end_pos = left_max_pos.second;
    avg_con_key = left_avg;
} else {
    attack::get_poison_key_boundaries_keys(dns[i],right_max_pos.first, 
            right_max_pos.second,
            poison_key_boundaries);
    start_pos = right_max_pos.first;
    end_pos = right_max_pos.second;
    avg_con_key = right_avg;
}


attack::get_poison_key_boundaries_keys(dns[i], start_pos, end_pos, poison_key_boundaries);
attack::gen_key(poison_key_boundaries.first, poison_key_boundaries.second, keys_1, 0, desired_key_size);


std::sort(keys_1, keys_1 + desired_key_size,
        [](auto const& a, auto const& b) { return a < b; });

if(!attack::validate_keys(keys_1, desired_key_size, dns[i],index)){
    outputFile1<<i
        <<std::endl;
    error_node++;
    continue;
}

int total_keys = dns[i]->num_keys_;
int temp_sideway = index.stats_.num_sideways_splits;
int temp_downward = index.stats_.num_downward_splits;

double node_cost =dns[i]->empirical_cost();
int action_flag = 0;
Cost c;
attack::insert(index, dns[i],keys_1, i, desired_key_size, node_cost, action_flag, &c);

int num_sideway_splits = index.stats_.num_sideways_splits - temp_sideway;
int num_downward_splits = index.stats_.num_downward_splits - temp_downward;

int split_flag = 0;

if(num_sideway_splits != 0 || num_downward_splits != 0) {
    split_flag = 1;
}

outputFile<< i <<","
<<start_pos<<","
<<end_pos<<","
<<end_pos - start_pos + 1<<","
<<avg_con_key<<","
<<whole_avg<<","
<<node_cost<<","
<<desired_key_size<<","
<<total_keys<<","
<<num_sideway_splits<<","
<<num_downward_splits<<","
<<split_flag<<","
<<action_flag
<<std::endl;
delete[] keys_1;
}

attack::print_index_stats(index);
std::cout<<"total error_node="<<error_node<<"\n";
std::cout<<"======targeting relative small max consective key ends======\n";

outputFile.close();
outputFile1.close();

delete[] keys;
delete[] values;
}


TEST_CASE("Test3") {

    std::cout<<"======targeting key density starts======\n";
    int total_num_keys = 20000000;
    int init_num_keys = total_num_keys / 2;

    auto keys = new KEY_TYPE[total_num_keys];
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

    attack::print_execution_stats("Initial");

    load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
    attack::print_execution_stats("Finishi load_binary_data()");

    std::mt19937_64 gen_payload(std::random_device{}());
    for (int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }


    attack::print_execution_stats("key-value pair generation");

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index.bulk_load(values, init_num_keys);
    attack::print_execution_stats("Finish bulk_load() insertions");
    attack::print_index_stats(index);

    for(int i = init_num_keys; i < total_num_keys; i++) {
        index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    attack::print_execution_stats("Finish insert() insertions");
    attack::print_index_stats(index);

    index.get_all_datanodes(dns);
    attack::get_parents(dns, parents, index);

    std::pair<int, int> con_left_pos, con_right_pos, gap_left_pos, gap_right_pos;
    std::pair<int, int> max_left_right_keys, max_left_right_gaps;
    std::pair<double, double> poison_key_boundaries;


    attack::print_index_stats(index);
    std::ofstream outputFile;
    outputFile.open("/root/data/ALEX/attack/density_longitudes_larger.txt");
    int error_node = 0;

    for(int i= 0; i < dns.size(); i++) {
        int desired_key_size = (attack::get_num_keys(dns[i]->expansion_threshold_) - dns[i]->num_keys_) + 1;
        auto keys_1 = new KEY_TYPE[desired_key_size];
        auto counts = new int[2];
        std::pair<double, double> key_boundaries;
        bool left = attack::should_insert_left_range(dns[i], parents[i], index, key_boundaries, counts);
        if(left) {
            attack::gen_key(key_boundaries.first, (key_boundaries.first + key_boundaries.second)/2, keys_1, 0, desired_key_size);
        } else {
            attack::gen_key((key_boundaries.first + key_boundaries.second)/2, key_boundaries.second, keys_1, 0, desired_key_size);
        }

        std::sort(keys_1, keys_1 + desired_key_size,
                [](auto const& a, auto const& b) { return a < b; });

        if(!attack::validate_keys(keys_1, desired_key_size, dns[i],index)){
            //outputFile1<<i<<","
            //    <<start_pos<<","
            //    <<end_pos
            //    <<std::endl;
            //
            error_node++;
            continue;
        }

        int total_keys = dns[i]->num_keys_;
        int temp_sideway = index.stats_.num_sideways_splits;
        int temp_downward = index.stats_.num_downward_splits;
        double node_cost = dns[i]->empirical_cost();
        int action_flag = 0;
        Cost c;
        attack::insert(index, dns[i],parents[i],keys_1, i, desired_key_size, node_cost, action_flag, &c);

        int num_sideway_splits = index.stats_.num_sideways_splits - temp_sideway;
        int num_downward_splits = index.stats_.num_downward_splits - temp_downward;

        outputFile<< i <<","
            <<desired_key_size<<","
            <<node_cost<<","
            <<num_sideway_splits<<","
            <<num_downward_splits<<","
            <<action_flag
            <<std::endl;
        delete[] keys_1;
        delete[] counts;
    }

    attack::print_index_stats(index);
    outputFile.close();
    std::cout<<"======targeting key density ends======\n";

    delete[] keys;
    delete[] values;
}




TEST_CASE("Test4") {

    std::cout<<"======targeting consective gap starts======\n";
    int total_num_keys = 20000000;
    int init_num_keys = total_num_keys / 2;

    auto keys = new KEY_TYPE[total_num_keys];
    auto values = new std::pair<KEY_TYPE, PAYLOAD_TYPE>[total_num_keys];
    std::vector<alex::AlexDataNode<KEY_TYPE, PAYLOAD_TYPE>*> dns;
    std::vector<alex::AlexModelNode<KEY_TYPE, PAYLOAD_TYPE>*> parents;

    attack::print_execution_stats("Initial");

    load_binary_data(keys, total_num_keys, "/root/data/ALEX/longitudes-200M.bin.data");
    attack::print_execution_stats("Finishi load_binary_data()");

    std::mt19937_64 gen_payload(std::random_device{}());
    for (int i = 0; i < init_num_keys; i++) {
        values[i].first = keys[i];
        values[i].second = static_cast<PAYLOAD_TYPE>(gen_payload());
    }

    attack::print_execution_stats("key-value pair generation");

    alex::Alex<KEY_TYPE, PAYLOAD_TYPE> index;
    std::sort(values, values + init_num_keys,
            [](auto const& a, auto const& b) { return a.first < b.first; });

    index.bulk_load(values, init_num_keys);
    attack::print_execution_stats("Finish bulk_load() insertions");
    attack::print_index_stats(index);

    for(int i = init_num_keys; i < total_num_keys; i++) {
        index.insert(keys[i], static_cast<PAYLOAD_TYPE>(gen_payload()));
    }

    attack::print_execution_stats("Finish insert() insertions");
    attack::print_index_stats(index);

    index.get_all_datanodes(dns);
    attack::get_parents(dns, parents, index);

    std::pair<int, int> con_left_pos, con_right_pos, gap_left_pos, gap_right_pos;
    std::pair<int, int> max_left_right_keys, max_left_right_gaps;
    std::pair<double, double> poison_key_boundaries;

    std::cout<<"size of datanodes="<<dns.size()
        <<" size of parents="<<parents.size()
        <<std::endl;
    std::ofstream outputFile, outputFile1;
    int temp = 0;
    int desired_key_size = 0;

    outputFile.open("/root/data/ALEX/attack/max_consecutive_gaps_stats_longitudes.txt");
    outputFile1.open("/root/data/ALEX/attack/gaps_log.txt");
    int node_error = 0;
    for(int i= 0; i < dns.size(); i++) {

        max_left_right_gaps = attack::get_max_gap(dns[i], parents[i], gap_left_pos, gap_right_pos);
        desired_key_size = attack::get_num_keys(dns[i]->expansion_threshold_) - dns[i]->num_keys_ + 1;


        auto keys_1 = new KEY_TYPE[desired_key_size];
        int start_pos = 0;
        int end_pos = 0;
        if(max_left_right_gaps.first > max_left_right_gaps.second) {
            attack::get_poison_key_boundaries_gaps(dns[i], 
                    gap_left_pos.first, gap_left_pos.second, poison_key_boundaries);
            start_pos = gap_left_pos.first;
            end_pos = gap_left_pos.second; 
        } else {
            attack::get_poison_key_boundaries_gaps(dns[i],gap_right_pos.first, 
                    gap_right_pos.second, poison_key_boundaries);
            start_pos = gap_right_pos.first;
            end_pos = gap_right_pos.second; 
        }

        attack::gen_key(poison_key_boundaries.first, poison_key_boundaries.second, keys_1, 0, desired_key_size);
        std::sort(keys_1, keys_1 + desired_key_size,
                [](auto const& a, auto const& b) { return a < b; });

        if(!attack::validate_keys(keys_1, desired_key_size, dns[i],index)){
            //outputFile1<<i<<","
            //    <<start_pos<<","
            //    <<end_pos
            //    <<std::endl;
            node_error++;
            continue;
        }

        int total_keys = dns[i]->num_keys_;
        int temp_sideway = index.stats_.num_sideways_splits;
        int temp_downward = index.stats_.num_downward_splits;

        double node_cost = dns[i]->empirical_cost();
        int action_flag = 0;
        Cost c;

        attack::insert(index, dns[i],keys_1, i, desired_key_size, node_cost, action_flag, &c);

        int num_sideway_splits = index.stats_.num_sideways_splits - temp_sideway;
        int num_downward_splits = index.stats_.num_downward_splits - temp_downward;

        int split_flag = 0;

        if(num_sideway_splits != 0 || num_downward_splits != 0) {
            split_flag = 1;
        }

        outputFile<< i <<","
            <<start_pos<<","
            <<end_pos<<","
            <<node_cost<<","
            <<desired_key_size<<","
            <<total_keys<<","
            <<num_sideway_splits<<","
            <<num_downward_splits<<","
            <<c.node_intra_node_cost_1<<","
            <<c.node_traversal_cost_1<<","
            <<c.node_intra_node_cost_2_left<<","
            <<c.node_intra_node_cost_2_right<<","
            <<c.node_traversal_cost_2
            <<std::endl;
        delete[] keys_1;
    }

    delete[] keys;
    delete[] values;
    outputFile.close();
    outputFile1.close();
    attack::print_index_stats(index);
    std::cout<<"======targeting consective gap ends======\n";
}*/
}
