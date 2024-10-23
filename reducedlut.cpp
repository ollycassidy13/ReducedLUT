/*
----------------------------------------------------
ReducedLUT
----------------------------------------------------
*/

#include "reducedlut.h"

int main(int argc, char* argv[]) {
    std::string filename = "logo.txt";
    std::ifstream inputFile(filename);

    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            std::cout << line << std::endl;
        }
        inputFile.close();
    }

    bool is_table = false;
    bool dc = false;
    std::string table_path;
    std::string input_path;
    std::string table_name = "reducedlut";
    std::string output_path = ".";
    reducedlut::struct_configs configs = {2, 1, 1, 0, 0};
    int rarity = 0;
    long int exiguity = 0;

    if ((argc < 2) || (argc % 2 == 0)) {
        reducedlut::help();
        return 1;
    }
   
    for (int i = 1; i < argc; i += 2) {
        std::string current_arg = argv[i];
        
        if (current_arg == "-table") {
            is_table = true;
            table_path = argv[i + 1];
        } else if (current_arg == "-input") {
            input_path = argv[i + 1];
            dc = true;
        } else if (current_arg == "-rarity") {
            rarity = std::stoi(argv[i + 1]); 
        } else if (current_arg == "-exiguity") {
            exiguity = std::stoi(argv[i + 1]); 
        } else if (current_arg == "-name") {
            table_name = argv[i + 1];
        } else if (current_arg == "-output") {
            output_path = argv[i + 1];
        } else if (current_arg == "-mdbw") {
            configs.mdbw = std::stoi(argv[i + 1]);
        } else if (current_arg == "-hbs") {
            configs.hbs = std::stoi(argv[i + 1]);
        } else if (current_arg == "-ssc") {
            configs.ssc = std::stoi(argv[i + 1]);
        } else if (current_arg == "-mlc") {
            configs.mlc = std::stoi(argv[i + 1]);
        } else if (current_arg == "-bits") {
            configs.bits = std::stoi(argv[i + 1]);
        }else {
            reducedlut::help();
            return 1;
        }
    }

    std::cout << "Results\n-------------------------------------------------------------------\n";
    std::vector<long int> table_data;
    std::vector<std::string> dont_care_values;
    bool is_signed = false;

    if (is_table) {
        std::ifstream table_file(table_path);
        if (table_file.is_open()) {
            std::string line;
            while (getline(table_file, line)) {
                table_data.push_back(std::stol(line, nullptr, 16));
            }
            table_file.close();

            if (table_data.size() != (1 << reducedlut::bit_width(table_data.size() - 1))) {
                std::cerr << "Error: The table size must be a power of 2." << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Error: Could not open the table file." << std::endl;
            return 1;
        }
    } else {
        reducedlut::help();
        return 1;
    }

    long int initial_size;
    std::vector<long int> final_size;

    int w_in = reducedlut::bit_width(table_data.size() - 1);
    int w_out = reducedlut::bit_width(*std::max_element(table_data.begin(), table_data.end()));
    std::cout << "--> Input Bit Width: " << w_in << std::endl;
    std::cout << "--> Output Bit Width: " << w_out  << std::endl;

    std::vector<long int> dont_care_indices;
    if (dc) {
        std::set<long int> unique_dont_care_indices;
        std::unordered_map<long int, int> frequency_map;  
        std::ifstream input_file(input_path);

        if (input_file.is_open()) {
            std::string line;
            while (std::getline(input_file, line)) {
                if (line.length() != static_cast<size_t>(w_in)) {  
                    std::cerr << "Error: Input value \"" << line << "\" is not " << w_in << " bits long." << std::endl;
                    return 1;
                }

                long int dont_care_index = std::stoi(line, nullptr, 2);
                if (dont_care_index < static_cast<long int>(table_data.size())) {
                    frequency_map[dont_care_index]++;
                } else {
                    std::cerr << "Error: Don't care index out of bounds." << std::endl;
                    return 1;
                }
            }
            input_file.close();
        } else {
            std::cerr << "Error: Could not open the input file." << std::endl;
            return 1;
        }

        for (long int i = 0; i < static_cast<long int>(table_data.size()); ++i) {
            if (frequency_map[i] < rarity) {
                unique_dont_care_indices.insert(i);
            }
        }

        dont_care_indices.assign(unique_dont_care_indices.begin(), unique_dont_care_indices.end());

        long int max_value_in_table = reducedlut::max_value_ignoring_dont_cares(table_data);
        int bit_width = reducedlut::bit_width(max_value_in_table);
    }

    reducedlut::reducedlut(table_data, table_name, output_path, configs, &initial_size, final_size, dc, w_out, dont_care_indices, exiguity);

    if (!final_size.empty()) {
        std::cout << "Information: The table was compressed successfully!" << std::endl;

        std::cout << "Information: The results are as follows." << std::endl;
        for (size_t i = 0; i < final_size.size(); i++) {
            std::cout << "-->  File Name: " << table_name << "_v" << i + 1 << "  Initial Size (bit): " << initial_size << "  Final Size (bit): " << final_size.at(i) << std::endl;
        }

    } else {
        std::cout << "Information: Unable to compress the table!" << std::endl;
    }

    if (configs.bits) {
        std::ofstream bits_file(output_path + "/../bits.txt", std::ios::app);
        if (bits_file.is_open()) {
            if (!final_size.empty()) {
                bits_file << "Initial Size (bit): " << initial_size << "\n";
                for (size_t i = 0; i < final_size.size(); i++) {
                    bits_file << "Final Size (bit) for " << table_name << "_v" << i + 1 << ": " << final_size.at(i) << "\n";
                }
                bits_file.close();
            }
            else {
                bits_file << "Initial Size (bit): " << initial_size << "\n";
                bits_file << "Final Size (bit) for " << table_name << "_v1" << ": " << initial_size << "\n";
                bits_file.close();
            }
        } 
        else {
            std::cerr << "Error: Could not open bits.txt for writing." << std::endl;
        }
    }

    return 0;
}

void reducedlut::reducedlut(std::vector<long int>& table_data, const std::string& table_name, const std::string& output_path, struct struct_configs configs, long int* initial_size, std::vector<long int>& final_size, bool dc, int out, std::vector<long int>& dont_care_indices, long int exiguity)
{
    std::vector<std::vector<long int>> all_cost;
    std::vector<int> all_w_in, all_w_out, all_w_l, all_w_s;
    std::vector<std::vector<long int>> all_t_lb, all_t_ust, all_t_bias, all_t_idx, all_t_rsh;
    
    std::vector<long int> t = table_data;
    std::vector<long int> bias_dc;

    int wo_ust = bit_width(*std::max_element(t.begin(), t.end()));
    int passes = 0;

    while (true) 
    {   
        passes++;
        // std::cout << passes << std::endl;
        int w_in = bit_width(t.size()-1);
        int w_out = bit_width(*std::max_element(t.begin(), t.end()));

        bool compressed = false;
        long int best_cost = (1 << w_in) * w_out;
        int best_w_l = 0, best_w_s = 0, best_rem_dc = 0;
        std::vector<long int> best_t_lb, best_t_ust, best_t_bias, best_t_idx, best_t_rsh, best_t_bias_dc;

        int max_w_l = (configs.hbs)? w_out-1 : 0;

        for (int w_l = 0; w_l <= max_w_l; w_l++) 
        {
            std::vector<long int> t_hb;
            std::vector<long int> t_lb(t.size());
            for(int i = 0; i < static_cast<int>(t.size()); i++) {
                t_hb.push_back(t.at(i) >> w_l);
                t_lb[i] = t.at(i) & (((long int)1 << w_l) - 1);
            }
            
            long int cost_t_lb = (1 << w_in) * bit_width(*std::max_element(t_lb.begin(), t_lb.end()));

            for(int w_s = configs.mdbw; w_s < w_in; w_s++) {
                std::vector<long int> t_ust, t_bias, t_idx, t_rsh, t_bias_dc;
                long int cost_t_hb;
                int yay = 0;
                if (passes == 1) {
                    cost_t_hb = hb_compression(configs.ssc, t_hb, w_s, t_ust, t_bias, t_bias_dc, t_idx, t_rsh, wo_ust, dc, dont_care_indices, w_l, yay, exiguity);
                } else {
                    cost_t_hb = hb_compression(configs.ssc, t_hb, w_s, t_ust, t_bias, t_bias_dc, t_idx, t_rsh, wo_ust, dc, bias_dc, w_l, yay, exiguity);
                }

                if ((cost_t_hb + cost_t_lb) < best_cost) {
                    compressed = true;
                    best_cost = cost_t_hb + cost_t_lb;
                    best_w_l = w_l;
                    best_w_s = w_s;
                    best_t_lb = t_lb;
                    best_t_ust = t_ust;
                    best_t_bias = t_bias;
                    best_t_idx = t_idx;
                    best_t_rsh = t_rsh;
                    best_t_bias_dc = t_bias_dc;
                    best_rem_dc = yay;
                }
            }
        }

        if (compressed) {
            std::vector<long int> cost = {(1 << w_in) * w_out, best_cost};
            all_cost.push_back(cost);
            all_w_in.push_back(w_in);
            all_w_out.push_back(w_out);
            all_w_l.push_back(best_w_l);
            all_w_s.push_back(best_w_s);
            all_t_lb.push_back(best_t_lb);
            all_t_ust.push_back(best_t_ust);
            all_t_bias.push_back(best_t_bias);
            all_t_idx.push_back(best_t_idx);
            all_t_rsh.push_back(best_t_rsh);
            std::cout << best_rem_dc << " dont cares remain out of " << dont_care_indices.size() << " origionally" << std::endl;

            if ((best_t_bias.size() >= (1 << (configs.mdbw + 1))) && configs.mlc) {
                t = best_t_bias; 

                bias_dc.clear(); 
                bias_dc = best_t_bias_dc; 
            } else {
                break;
            }
        } else {
            break;
        }
    }

    if (!all_w_in.empty()) {
        *initial_size = all_cost[0][0];
        final_size.push_back(all_cost[0][1]);
        for (size_t i = 1; i < all_w_in.size(); i++) {
            final_size.push_back(final_size[i - 1] - all_cost[i][0] + all_cost[i][1]);
        }

        for (size_t i = 0; i < all_w_in.size(); i++) {
            std::string rtl_file_path = output_path + "/" + table_name + "_v" + std::to_string(i + 1) + ".v";
            rtl(rtl_file_path, table_name, all_w_in, all_w_out, all_w_l, all_w_s, all_t_lb, all_t_ust, all_t_bias, all_t_idx, all_t_rsh, static_cast<int>(i + 1), wo_ust);
        }
    }
}

long int reducedlut::hb_compression(bool ssc, const std::vector<long int>& t_hb, int w_s, std::vector<long int>& t_ust, std::vector<long int>& t_bias, std::vector<long int>& t_bias_dc, std::vector<long int>& t_idx, std::vector<long int>& t_rsh, const int wo_ust, bool dc, const std::vector<long int>& dont_care_indices, const int w_l, int& yay, long int exiguity)
{
    const int w_in = bit_width(t_hb.size() - 1);
    const int w_out = wo_ust;
    const long int num_sub_table = (1 << (w_in - w_s));
    const int len_sub_table = (1 << w_s);

    std::vector<long int> t_st(num_sub_table * len_sub_table, 0);
    std::vector<bool> boolean_mask(num_sub_table * len_sub_table, 0);

    for (long int i = 0; i < num_sub_table; i++) {
        long int bias = *std::min_element(t_hb.begin() + i * len_sub_table, t_hb.begin() + (i + 1) * len_sub_table);
        t_bias.push_back(bias);
        for (int j = 0; j < len_sub_table; j++) {
            t_st[i * len_sub_table + j] = t_hb[i * len_sub_table + j] - bias;
        }
    }

    t_bias_dc.clear();
    if (dc) {
        for (auto idx : dont_care_indices) {
            boolean_mask[idx] = 1;
        }
    }

    // Bias dc
    
    for (long int i = 0; i < num_sub_table; i++) {
        bool all_dont_care = true;
        for (int j = 0; j < len_sub_table; j++) {
            if (boolean_mask[i * len_sub_table + j] != 1) {
                all_dont_care = false;
                break;
            }
        }
        if (all_dont_care) {
            t_bias_dc.push_back(i); 
        }
    }
    
    if (!ssc) {
        t_ust = t_st;
        int w_bias = bit_width(*std::max_element(t_bias.begin(), t_bias.end()));
        if (w_bias != 0) {
            w_bias = wo_ust - w_l;
        }
        int w_ust = wo_ust - w_l;
        return (1 << w_in) * w_ust + (1 << (w_in - w_s)) * w_bias; 
    } else {
        std::vector<std::vector<bool>> sm(num_sub_table, std::vector<bool>(num_sub_table, false));
        std::vector<std::vector<int>> sm_rsh(num_sub_table, std::vector<int>(num_sub_table, 0));
        std::vector<long int> sv(num_sub_table);

        for(long int i = 0; i < num_sub_table; i++)
        {
            for(long int j = i+1; j < num_sub_table; j++)
            {
                for(int rsh = 0; rsh < 4; rsh++)
                {
                    bool i_generates_j = true, j_generates_i = true;
                    for (int p = 0; p < len_sub_table; p++)
                    {
                        long int value_i = t_st.at(i*len_sub_table+p);
                        long int value_j = t_st.at(j*len_sub_table+p);
                        i_generates_j = i_generates_j && ((value_i >> rsh) == value_j);
                        j_generates_i = j_generates_i && ((value_j >> rsh) == value_i);
                        if(i_generates_j == false && j_generates_i == false)
                            break;
                    }
                    
                    if(i_generates_j)
                    {
                        sm.at(i).at(j) = 1;
                        sv.at(i)++;
                        sm_rsh.at(i).at(j) = rsh;
                        
                    }
                    if(j_generates_i)
                    {
                        sm.at(j).at(i) = 1;
                        sv.at(j)++;
                        sm_rsh.at(j).at(i) = rsh;
                    }
                    if(i_generates_j || j_generates_i)
                        break;
                }

            }
        }

        if (dc) {
            struct Modification {
                long int index;
                long int original_value;
            };
            std::stack<Modification> modifications;

            struct similarity {
                long int gen_j;
                long int k;
                int rsh; 
            };

            std::vector<long int> sorted_j_indices(num_sub_table);

            std::vector<long int> sorted_i_indices = get_unique_indices(sv, sm, sm_rsh, t_idx, t_rsh, t_st, t_ust, num_sub_table, len_sub_table);

            std::sort(sorted_i_indices.begin(), sorted_i_indices.end(), [&](long int a, long int b) {
                return sv[a] < sv[b];
            });

            for (long int idx_i : sorted_i_indices) {

                for (long int idx = 0; idx < num_sub_table; idx++) {
                    sorted_j_indices[idx] = idx;
                }

                std::sort(sorted_j_indices.begin(), sorted_j_indices.end(), [&](long int a, long int b) {
                    return sv[a] > sv[b];
                });
                long int i = idx_i;
                if (sv[i] == 0) {
                    bool modified = false;
                    bool found = false;  

                    for (long int sorted_idx : sorted_j_indices) {
                        long int j = sorted_idx;

                        if (i != j) {
                            for (int rsh = 0; rsh < 4; rsh++) {
                                while (!modifications.empty()) {
                                    modifications.pop();
                                }

                                bool can_generate = true;

                                for (int p = 0; p < len_sub_table; p++) {
                                    long int value_i = t_st[i * len_sub_table + p];
                                    long int value_j = t_st[j * len_sub_table + p];

                                    if (boolean_mask[i * len_sub_table + p] && (value_j << rsh) <= max_value_for_bitwidth(wo_ust - w_l)) {
                                        modifications.push({i * len_sub_table + p, value_i});
                                        t_st.at(i * len_sub_table + p) = (value_j << rsh) + (value_j & 1);
                                        modified = true;
                                    }

                                    can_generate = can_generate && ((t_st.at(i * len_sub_table + p) >> rsh) == value_j);

                                    if (!can_generate) {
                                        while (!modifications.empty()) {
                                            Modification mod = modifications.top();
                                            modifications.pop();
                                            t_st.at(mod.index) = mod.original_value;
                                        }
                                        break;
                                    }
                                }

                                if (can_generate) {
                                    sv[i] = 1;
                                    for (int b = 0; b < num_sub_table; b++) {
                                        if (b != i) {
                                            sm[b][i] = 0;
                                        }
                                    }
                                    sm[j][i] = true;
                                    sv[j]++;
                                    sm_rsh[j][i] = rsh;
                                    for (int q = 0; q < len_sub_table; q++) {
                                        boolean_mask.at(i * len_sub_table + q) = 0;
                                    }
                                    found = true;  
                                    break;
                                }
                            }
                        }

                        if (found) {  
                            break;
                        }
                    }
                }
                
                else if (sv[i] >= 1 && sv[i] <= exiguity) {
                    bool found = false;

                    for (long int idx_j : sorted_j_indices) {
                        long int j = idx_j;

                        if (i != j) {
                            std::stack<Modification> local_modifications;
                            bool can_generate = true;
                            int rsh_o = 0;
                            for (int rsh = 0; rsh < 4; rsh++) {
                                rsh_o = rsh;
                                can_generate = true;

                                for (int p = 0; p < len_sub_table; p++) {
                                    long int value_i = t_st[i * len_sub_table + p];  
                                    long int value_j = t_st[j * len_sub_table + p];

                                        if (boolean_mask[i * len_sub_table + p] && (value_j << rsh) <= max_value_for_bitwidth(wo_ust - w_l)) {
                                            local_modifications.push({i * len_sub_table + p, value_i});
                                            t_st.at(i * len_sub_table + p) = (value_j << rsh) + (value_j & 1);
                                        }

                                    can_generate = can_generate && ((t_st.at(i * len_sub_table + p) >> rsh) == value_j);

                                    if (!can_generate) {
                                        while (!local_modifications.empty()) {
                                            Modification mod = local_modifications.top();
                                            local_modifications.pop();
                                            t_st.at(mod.index) = mod.original_value;
                                        }
                                        break;
                                    }
                                }
                                if (can_generate) {
                                    break;
                                }
                            }
                            if (can_generate) {
                                bool can_move = true;
                                std::stack<Modification> total_modifications = local_modifications;
                                std::vector<similarity> matrix_updates;
                                bool moveall = true;
                                
                                for (long int gen_j: sorted_i_indices) {
                                    if (!sm[i][gen_j] || gen_j == idx_j) continue;

                                    else {
                                        for (long int k_idx : sorted_j_indices) {
                                            long int k = k_idx;
                                            if (k == gen_j) continue;
                                            else if (sm[k][gen_j]) {
                                                continue;
                                            }
                                            else {
                                                for (int rsh = 0; rsh < 4; rsh++) {
                                                    can_move = true;
                                                    for (int p = 0; p < len_sub_table; p++) {
                                                        long int value_gen_j = t_st[gen_j * len_sub_table + p];
                                                        long int value_k = t_st[k * len_sub_table + p];

                                                        if (boolean_mask[gen_j * len_sub_table + p] && (value_k << rsh) <= max_value_for_bitwidth(wo_ust - w_l)) {
                                                            total_modifications.push({gen_j * len_sub_table + p, value_gen_j});
                                                            t_st.at(gen_j * len_sub_table + p) = (value_k << rsh) + (value_k & 1);
                                                        }
                                                        
                                                        can_move = can_move && ((value_k >> rsh) == value_gen_j);

                                                        if (!can_move) {
                                                            break;
                                                        }
                                                    }
                                                    if (can_move) {
                                                        matrix_updates.push_back({gen_j, k, rsh});
                                                        break;
                                                    }
                                                }
                                                if (can_move) {
                                                    break;
                                                }
                                                else {
                                                    moveall = false;
                                                }
                                            } 
                                        }
                                    }
                                }
                                if (moveall) {
                                    sv[i] = 0;
                                    for (int a = 0; a < num_sub_table; a++) {
                                        if (a != i) {
                                            sm[a][i] = 0;
                                        }
                                    }
                                    sv[j]++;
                                    sm[j][i] = 1;
                                    sm_rsh[j][i] = rsh_o;
                                    for (int p = 0; p < len_sub_table; p++) {
                                        boolean_mask.at(i * len_sub_table + p) = 0;
                                    }
                                    for (int a = 0; a < matrix_updates.size(); a++) {
                                        similarity update = matrix_updates[a];
                                        sv[update.gen_j] = 0;
                                        for (int b = 0; b < num_sub_table; b++) {
                                            if (b != update.gen_j) {
                                                sm[b][update.gen_j] = 0;
                                            }
                                        }
                                        sv[update.k]++;
                                        sm[update.k][update.gen_j] = 1;
                                        sm_rsh[update.k][update.gen_j] = update.rsh;
                                        for (int p = 0; p < len_sub_table; p++) {
                                            boolean_mask.at(update.gen_j * len_sub_table + p) = 0;
                                        }
                                    }
                                    while(!total_modifications.empty()) {
                                        total_modifications.pop();
                                    }
                                    matrix_updates.clear();
                                    found = true;
                                    // recalc sv
                                    for (int i = 0; i < num_sub_table; ++i) {
                                        long int row_sum = 0;
                                        for (int j = 0; j < num_sub_table; ++j) {
                                            if (sm[j][i]) { // order i and j
                                                row_sum += 1;
                                            }
                                        }
                                        sv[i] = row_sum - 1;
                                    }
                                    break;
                                }
                                else {
                                    while (!total_modifications.empty()) {
                                        Modification mod = total_modifications.top();
                                        total_modifications.pop();
                                        t_st.at(mod.index) = mod.original_value;
                                    }
                                    matrix_updates.clear();
                                }
                            }
                        }
                        if (found) {
                            break;
                        }
                    }
                }
            }
        }

        yay = 0;
        for (int i = 0; i < boolean_mask.size(); i++) {
            if (boolean_mask[i]) {
                yay++;
            }
        }

        t_idx.resize(num_sub_table);
        t_rsh.resize(num_sub_table);
        std::vector<bool> sub_table_derived(num_sub_table, false);
        long int unique_idx = 0;

        while (std::any_of(sv.begin(), sv.end(), [](int v) { return v > 0; })) {
            long int idx_unique = std::distance(sv.begin(), std::max_element(sv.begin(), sv.end()));
            t_idx[idx_unique] = unique_idx;
            t_rsh[idx_unique] = 0;
            sub_table_derived[idx_unique] = true;
            for (long int idx = 0; idx < num_sub_table; idx++) {
                if (sm[idx_unique][idx]) {
                    t_idx[idx] = unique_idx;
                    t_rsh[idx] = sm_rsh[idx_unique][idx];
                    sub_table_derived[idx] = true;
                    std::fill(sm[idx].begin(), sm[idx].end(), false);
                    sv[idx] = 0;
                    for (long int i = 0; i < num_sub_table; i++) {
                        if (sm[i][idx]) {
                            sm[i][idx] = false;
                            sv[i]--;
                        }
                    }
                }
            }

            for (long int i = 0; i < num_sub_table; i++) {
                if (sm[i][idx_unique]) {
                    sm[i][idx_unique] = false;
                    sv[i]--;
                }
            }
            sv[idx_unique] = 0;

            t_ust.insert(t_ust.end(), t_st.begin() + idx_unique * len_sub_table, t_st.begin() + (idx_unique + 1) * len_sub_table);

            unique_idx++;
        }

        for (long int i = 0; i < num_sub_table; i++) {
            if (!sub_table_derived[i]) {
                t_idx[i] = unique_idx;
                t_rsh[i] = 0;

                t_ust.insert(t_ust.end(), t_st.begin() + i * len_sub_table, t_st.begin() + (i + 1) * len_sub_table);

                unique_idx++;
            }
        }
    }

    int w_ust = wo_ust - w_l;
    int w_bias = bit_width(*std::max_element(t_bias.begin(), t_bias.end()));
    if (w_bias != 0) {
        w_bias = wo_ust - w_l;
    }
    int w_rsh = bit_width(*std::max_element(t_rsh.begin(), t_rsh.end()));
    int w_idx = bit_width(*std::max_element(t_idx.begin(), t_idx.end()));
    return t_ust.size() * w_ust + (1 << (w_in - w_s)) * (w_bias + w_rsh + w_idx);
}

long int reducedlut::max_value_for_bitwidth(unsigned int bitwidth) {
    if (bitwidth >= sizeof(long int) * 8) {
        return std::numeric_limits<long int>::max();
    }
    if (std::is_signed<long int>::value) {
        return (1L << (bitwidth - 1)) - 1; 
    } else {
        return (1UL << bitwidth) - 1; 
    }
}

std::vector<long int> reducedlut::get_unique_indices(std::vector<long int> sv, 
                                         std::vector<std::vector<bool>> sm, 
                                         std::vector<std::vector<int>> sm_rsh, 
                                         std::vector<long int> t_idx, 
                                         std::vector<long int> t_rsh, 
                                         std::vector<long int> t_st, 
                                         std::vector<long int> t_ust, 
                                         long int num_sub_table, 
                                         long int len_sub_table) {
    std::vector<long int> unique_indices;
    long int unique_idx = 0;
    t_idx.resize(num_sub_table);
    t_rsh.resize(num_sub_table);
    std::vector<bool> sub_table_derived(num_sub_table, false);
    
    while (std::any_of(sv.begin(), sv.end(), [](int v) { return v > 0; })) {
        long int idx_unique = std::distance(sv.begin(), std::max_element(sv.begin(), sv.end()));
        unique_indices.push_back(idx_unique);  
        
        t_idx[idx_unique] = unique_idx;
        t_rsh[idx_unique] = 0;
        sub_table_derived[idx_unique] = true;
        for (long int idx = 0; idx < num_sub_table; idx++) {
            if (sm[idx_unique][idx]) {
                t_idx[idx] = unique_idx;
                t_rsh[idx] = sm_rsh[idx_unique][idx];
                sub_table_derived[idx] = true;
                std::fill(sm[idx].begin(), sm[idx].end(), false);
                sv[idx] = 0;

                for (long int i = 0; i < num_sub_table; i++) {
                    if (sm[i][idx]) {
                        sm[i][idx] = false;
                        sv[i]--;
                    }
                }
            }
        }

        for (long int i = 0; i < num_sub_table; i++) {
            if (sm[i][idx_unique]) {
                sm[i][idx_unique] = false;
                sv[i]--;
            }
        }

        sv[idx_unique] = 0;

        t_ust.insert(t_ust.end(), t_st.begin() + idx_unique * len_sub_table, t_st.begin() + (idx_unique + 1) * len_sub_table);

        unique_idx++;
    }
    return unique_indices;
}

void reducedlut::rtl(const string& file_path, const string& table_name, const vector<int>& all_w_in, const vector<int>& all_w_out, const vector<int>& all_w_l, const vector<int>& all_w_s, const vector<vector<long int>>& all_t_lb, const vector<vector<long int>>& all_t_ust, const vector<vector<long int>>& all_t_bias, const vector<vector<long int>>& all_t_idx, const vector<vector<long int>>& all_t_rsh, int max_level, const int wo_ust)
{
    ofstream file_init(file_path);
    file_init.close();

    for (int level = max_level; level >= 1; level--) 
    {
        const vector<long int> t_lb = all_t_lb.at(level-1);
        const vector<long int> t_ust = all_t_ust.at(level-1);
        const vector<long int> t_bias = all_t_bias.at(level-1);
        const vector<long int> t_idx = all_t_idx.at(level-1);
        const vector<long int> t_rsh = all_t_rsh.at(level-1);

        const int w_in = all_w_in.at(level-1);
        const int w_out = wo_ust;
        const int w_l = all_w_l.at(level-1);
        const int w_s = all_w_s.at(level-1);

        const int w_lb = (t_lb.size() == 0)? 0 : bit_width(*max_element(t_lb.begin(), t_lb.end()));
        int w_ust = wo_ust;
        if (w_l != 0) {
            w_ust = w_ust - w_l;
        }
        int w_bias = (t_bias.size() == 0)? 0 : bit_width(*max_element(t_bias.begin(), t_bias.end()));
        if (w_bias != 0) {
            w_bias = w_ust;
        }
        const int w_idx = (t_idx.size() == 0)? 0 : bit_width(*max_element(t_idx.begin(), t_idx.end()));
        const int w_rsh = (t_rsh.size() == 0)? 0 : bit_width(*max_element(t_rsh.begin(), t_rsh.end()));

        if (w_ust != 0) {
            string name = table_name + "_ust_" + to_string(level);
            plaintable_rtl(file_path, name, t_ust, true, w_ust);
        }
        if (level == max_level && w_bias != 0) {
            string name = table_name + "_bias_" + to_string(level);
            plaintable_rtl(file_path, name, t_bias, true, w_bias);
        }
        if (w_idx != 0) 
        {
            string name = table_name + "_idx_" + to_string(level);
            plaintable_rtl(file_path, name, t_idx, false, 0);
        }
        if (w_rsh != 0) 
        {
            string name = table_name + "_rsh_" + to_string(level);
            plaintable_rtl(file_path, name, t_rsh, false, 0);
        }
        if (w_lb != 0) 
        {
            string name = table_name + "_lb_" + to_string(level);
            plaintable_rtl(file_path, name, t_lb, false, 0);
        }

        ofstream file(file_path, ios::app);

        if(level == 1)
            file << "\nmodule " << table_name << "(address, data);\n";
        else 
            file << "\nmodule " << table_name << "_" << level << "(address, data);\n";
        file << "input wire [" << w_in - 1 << ":0] address;\n";
        file << "output reg [" << w_out - 1 << ":0] data;\n\n";


        if(w_idx != 0)
            file << "wire [" << w_idx - 1 << ":0] i; " << table_name << "_idx_" << level << " idx_" << level << "_inst(address[" << w_in - 1 << ":" << w_s << "], i);\n";
      
        if (w_rsh != 0)
            file << "wire [" << w_rsh - 1 << ":0] t; " << table_name << "_rsh_" << level << " rsh_" << level << "_inst(address[" << w_in - 1 << ":" << w_s << "], t);\n";

        if (w_bias != 0) 
        {
            if (level == max_level) 
                file << "wire [" << w_bias - 1 << ":0] b; " << table_name << "_bias_" << level << " bias_" << level << "_inst(address[" << w_in - 1 << ":" << w_s << "], b);\n";
            else
                file << "wire [" << w_bias - 1 << ":0] b; " << table_name << "_" << level + 1 << " " << table_name << "_" << level + 1  << "_inst(address[" << w_in - 1 << ":" << w_s << "], b);\n";
        }

        if (w_lb != 0)
            file << "wire [" << w_lb - 1 << ":0] lb; " << table_name << "_lb_" << level << " lb_" << level << "_inst(address, lb);\n";

        if (w_ust != 0) 
        {
            if (w_idx != 0) 
                file << "wire [" << w_ust - 1 << ":0] u; " << table_name << "_ust_" << level << " ust_" << level << "_inst({i, address[" << w_s - 1 << ":0]}, u);";
            else 
            {
                if(t_idx.size() == 0)
                    file << "wire [" << w_ust - 1 << ":0] u; " << table_name << "_ust_" << level << " ust_" << level << "_inst(address, u);";
                else
                    file << "wire [" << w_ust - 1 << ":0] u; " << table_name << "_ust_" << level << " ust_" << level << "_inst(address[" << w_s - 1 << ":0], u);";
            }
        }

        file << "\n\nalways @(*) begin\n";

        if (w_l != 0)
            file << "\tdata = {";
        else 
            file << "\tdata = ";
        
        if (w_ust != 0 && w_rsh != 0 && w_bias != 0)
            file << "((u >> t) + b)";
        else if (w_ust != 0 && w_rsh != 0 && w_bias == 0) 
            file << "(u >> t)";
        else if (w_ust != 0 && w_rsh == 0 && w_bias != 0)
            file << "(u + b)";
        else if (w_ust != 0 && w_rsh == 0 && w_bias == 0)
            file << "u";
        else if (w_ust == 0 && w_bias != 0)
            file << "b";
        else 
            file << "'" << (w_out - w_l) << "'d0";

        if (w_l != 0) 
            if(w_l == w_lb)
                file << ", lb};\n";
            else if(w_lb != 0)
                file << ", " << (w_l - w_lb) << "'d0, lb};\n";
            else
                file << ", " << w_l << "'d0};\n";
        else
            file << ";\n";

        file << "end\nendmodule\n";

        file.close();
    }
}

void reducedlut::plaintable_rtl(const string& file_path, const string& table_name, const vector<long int>& table_data, bool set, const int width) 
{
    ofstream file(file_path, ios::app);
    int w_out = 0;
    int w_in = bit_width(table_data.size() - 1);
    if (set) {
        w_out = width;
    }
    else {
        w_out = bit_width(*max_element(table_data.begin(), table_data.end()));
    }

    file << "\nmodule " << table_name << "(address, data);\n";
    file << "input wire [" << w_in - 1 << ":0] address;\n";
    file << "output reg [" << w_out - 1 << ":0] data;\n\n";
    file << "always @(*) begin\n";
    file << "\tcase(address)\n";
    for (long int i = 0; i < table_data.size(); i++) 
    {
        file << "\t\t" << w_in << "'d" << i << ": ";
        file << "data = " << w_out << "'d" << table_data.at(i) << ";\n";
    }
    file << "\t\tdefault: data = " << w_out << "'d0" << ";\n\tendcase\nend\n";
    file << "endmodule\n";

    file.close();
}

long int reducedlut::max_value_ignoring_dont_cares(const std::vector<long int>& values) {
    long int max_val = LONG_MIN;
    for (const auto& val : values) {
        if (val != LONG_MAX && val > max_val) {
            max_val = val;
        }
    }
    return max_val;
}

int reducedlut::bit_width(long int value) 
{
    if (value == 0)
        return 0;
    return floor(log2(abs(value)))+1;
}

int reducedlut::bit_width_signed(long int min_value, long int max_value) 
{
    int w = 1;
    while(1)
    {
        if(min_value >= -((long int)1 << (w-1)) && max_value <= (((long int)1 << (w-1))-1))
            return w;
        else
            w++;
    }
}

void reducedlut::help() 
{
    string filename = "help.txt";
    ifstream inputFile(filename);
    if (inputFile.is_open()) 
    {
        string line;
        while (std::getline(inputFile, line)) 
        {
            std::cout << line << std::endl;
        }
        inputFile.close();
    }
}