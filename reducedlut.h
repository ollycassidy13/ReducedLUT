/*
----------------------------------------------------
ReducedLUT
----------------------------------------------------
*/
#  ifndef __COMPRESSEDLUT_H
#  define __COMPRESSEDLUT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <climits>  
#include <numeric>
#include <utility>
#include <stack>
#include <cstdlib> 
#include <ctime> 
#include <unordered_map>
#include "exprtk.hpp"

using namespace std;

namespace reducedlut {
    struct struct_configs {int mdbw;  bool hbs; bool ssc; bool mlc;};

    void reducedlut(vector<long int>& table_data, const string& table_name, const string& output_path, struct struct_configs configs, long int* initial_size, vector<long int>& final_size, const bool dc, int out, std::vector<long int>& dont_care_indices, long int exiguity);
    long int hb_compression(bool ssc, const vector<long int>& t_hb, int w_s, vector<long int>& t_ust, vector<long int>& t_bias, vector<long int>& t_bias_dc, vector<long int>& t_idx, vector<long int>& t_rsh, const int wo_ust, const bool dc, const std::vector<long int>& dont_care_indices, const int w_l, int& yay, long int exiguity);
    long int max_value_ignoring_dont_cares(const std::vector<long int>& values);
    void rtl(const string& file_path, const string& table_name, const vector<int>& all_w_in, const vector<int>& all_w_out, const vector<int>& all_w_l, const vector<int>& all_w_s, const vector<vector<long int>>& all_t_lb, const vector<vector<long int>>& all_t_ust, const vector<vector<long int>>& all_t_bias, const vector<vector<long int>>& all_t_idx, const vector<vector<long int>>& all_t_rsh, int max_level, const int wo_ust);
    void plaintable_rtl(const string& file_path, const string& table_name, const vector<long int>& table_data, bool set, const int width);
    int bit_width(long int value);
    int bit_width_signed(long int min_value, long int max_value);
    std::vector<long int> get_unique_indices(std::vector<long int> sv, 
                                        std::vector<std::vector<bool>> sm, 
                                        std::vector<std::vector<int>> sm_rsh, 
                                        std::vector<long int> t_idx, 
                                        std::vector<long int> t_rsh, 
                                        std::vector<long int> t_st, 
                                        std::vector<long int> t_ust, 
                                        long int num_sub_table, 
                                        long int len_sub_table);
    long int max_value_for_bitwidth(unsigned int bitwidth);
    void help();
}

# endif
