#ifndef PEAK_EXAMINING
#define PEAK_EXAMINING
#include <vector>
#include <math.h>

#include "optimize_selection.h"

vector<double> selection_opt::examine_peaks(){
    vector<double> ans;

    double interaction_distance = 0.05;
    double cut_off = 0;
    int local_maxima_radius = 4;

    vector<int> peaks;
    vector<bool> isolated;

    

    for(int i = local_maxima_radius; i < options.gss_out_pos.size() - local_maxima_radius; i++){
        bool at_least_as_high = true;
        for(int j = -local_maxima_radius; j <= local_maxima_radius; j++){
            if(options.gss_out_llr[i] < options.gss_out_llr[i + j])
                at_least_as_high = false;
        }
        if(at_least_as_high && options.gss_out_llr[i] >= cut_off){
            peaks.push_back(i);
            if(peaks.size() > 1){
                if( options.gss_out_pos[peaks[peaks.size() - 1]] - 
                    options.gss_out_pos[peaks[peaks.size() - 2]] < interaction_distance){
                        isolated[isolated.size() - 1] = false;
                        isolated.push_back(false);
                    }
            }else{
                isolated.push_back(true);
            }
        }
    }

    
    for(int i = 0; i < peaks.size(); i++){
        cout << "peak: " << peaks[i] << " " << isolated[i] << " " << options.gss_out_pos[peaks[i]] << "\n";
    }
    
    if(peaks.size() == 0){
        cerr << "\nNo peaks found in gss output.\n";
        return ans;
    }

    
    


    double nelder_mead_reflection  = 1;
    double nelder_mead_contraction = 0.5;
    double nelder_mead_expansion   = 2;
    double nelder_mead_shrinkage   = 0.5;

    nelder_mead optimizer(
        nelder_mead_reflection,
        nelder_mead_contraction,
        nelder_mead_expansion,
        nelder_mead_shrinkage
    );

    context = *this;

    vector<double> empty(0);                    
    double neutral_lnl = to_be_optimized(empty);
    
    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;

    double chrom_size = 0;
    for(uint i = 0; i < n_recombs.size(); i++){
        chrom_size += n_recombs[i];
    }

    if(false){
        double s = options.gss_out_sel[peaks[0]];
        //double s = 0.0119829190899;

        vector<double> poss(6);
        poss[0] = 0.143197;
        poss[1] = 0.973179;
        poss[2] = 1.00108;
        poss[3] = 0.199923;
        poss[4] = 0.995748;
        poss[5] = 1.01841;
        double poss_lnl = to_be_optimized(poss);

        vector<double> poss2(3);
        poss2[0] = 0.21021;
        poss2[1] = (1-s)/(1-s/2);
        poss2[2] = 1/(1-s/2);
        
        double poss2_lnl = to_be_optimized(poss2);

        vector<double> poss3(6);
        poss3[0] = 0.2;
        poss3[1] = 0.98989898989;
        poss3[2] = 1.0101010;
        poss3[3] = 0.21;
        poss3[4] = 0.98989898989;
        poss3[5] = 1.0101010;
        double poss3_lnl = to_be_optimized_only_near_sites(poss3);

        vector<double> real(6);
        real[0] = 0.2;
        real[1] = 0.98989898989;
        real[2] = 1.0101010;
        real[3] = 0.21;
        real[4] = 0.98989898989;
        real[5] = 1.0101010;
        double real_lnl = to_be_optimized(real);
        cerr << "\n";

        cerr << "Neut point  lnl: \t" << neutral_lnl << "\n";
        cerr << "Poss1 point lnl: \t" << poss_lnl << "\n";
        cerr << "Poss2 point lnl: \t" << poss2_lnl - neutral_lnl<< "\n";
        cerr << "Poss3 point lnl: \t" << poss3_lnl << "\n";
        cerr << "Real point  lnl: \t" << real_lnl  - neutral_lnl<< "\n";
    }

    vector<vector<vector<double>>> bottle_necks;

    vector<vector<double>> shallow;

    vector<double> shallow_short(4);
    shallow_short[0] = 5;
    shallow_short[1] = 0.03;
    shallow_short[2] = 0.01;
    shallow_short[3] = 20;

    vector<double> shallow_tall(4);
    shallow_tall[0] = 5;
    shallow_tall[1] = 0.03;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 20;

    shallow.push_back(shallow_short);
    shallow.push_back(shallow_tall);

    vector<vector<double>> deep;

    vector<double> deep_short(4);
    deep_short[0] = 5;
    deep_short[1] = 0.01;
    deep_short[2] = 0.005;
    deep_short[3] = 5;

    vector<double> deep_tall(4);
    deep_tall[0] = 5;
    deep_tall[1] = 0.01;
    deep_tall[2] = 0.01;
    deep_tall[3] = 5;

    deep.push_back(deep_short);
    deep.push_back(deep_tall);

    bottle_necks.push_back(shallow);
    bottle_necks.push_back(deep);



    vector<vector<vector<double>>> two_site_bottle_necks;

    //vector<vector<double>> shallow;

    //vector<double> shallow_short(4);
    shallow_short[0] = 3;
    shallow_short[1] = 0.01;
    shallow_short[2] = 0.01;
    shallow_short[3] = 20;

    //vector<double> shallow_tall(4);
    shallow_tall[0] = 3;
    shallow_tall[1] = 0.01;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 20;

    shallow[0] = shallow_short;
    shallow[1] = shallow_tall;

    //vector<vector<double>> deep;

    //vector<double> deep_short(4);
    deep_short[0] = 3;
    deep_short[1] = 0.005;
    deep_short[2] = 0.005;
    deep_short[3] = 5;

    //vector<double> deep_tall(4);
    deep_tall[0] = 3;
    deep_tall[1] = 0.005;
    deep_tall[2] = 0.01;
    deep_tall[3] = 5;

    deep[0] = deep_short;
    deep[1] = deep_tall;

    two_site_bottle_necks.push_back(shallow);
    two_site_bottle_necks.push_back(deep);


    //FIX replace 1 with peaks.size()
    for(int i = 0; i < 1; i++){
        //FIX replace true with isolated[i]
        if(true){
            //if isolated, we're gonna do two things
            // 1 see if the selection is non-additive
            // 2 see if there might be two sites here instead of 1

            //////
                //Seeing if site is non-additive

            //FIX replace false with !options.restrict_to_dominance
            if(false){

                double s = options.gss_out_sel[peaks[i]];

                vector<double> starting_parameters(3);
                starting_parameters[0] = options.gss_out_pos[peaks[i]];
                starting_parameters[1] = 1/(1-s/2);
                starting_parameters[2] = (1-s)/(1-s/2);

                vector<vector<double>> sites = parameters_to_sites(starting_parameters);

                

                vector<double> best_parameters = multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast
                );


                
                double new_lnl = to_be_optimized(best_parameters);

                cout << "\n\nGss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio: " << new_lnl - neutral_lnl << "\n";
                cout << "site:\t" << sites[0][0] << "\t" << sites[sites.size() - 1][1] << ",1," << sites[sites.size() - 1][2] << " or " << sites[sites.size() - 1][1]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2]) << ","<< 1/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2]) <<"," << sites[sites.size() - 1][2]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2]) << "\n";
            
            }
            
            //FIX replace false with else and connect to last statement
            if(true){
                
                vector<double> dom_starting_parameters(2);
                //FIX replace 0.2 with options.gss_out_pos[peaks[i]]
                dom_starting_parameters[0] = 0.2;
                dom_starting_parameters[1] = 1;

                vector<double> add_starting_parameters(2);
                add_starting_parameters[0] = 0.2;
                add_starting_parameters[1] = 0;

                //FIX uncomment
                vector<vector<double>> sites = parameters_to_sites(dom_starting_parameters, 2);
                
                vector<double> best_parameters = multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast_dom0,
                    2
                );

                double new_lnl = to_be_optimized_pop0_dominant(best_parameters);

                cout << "\n\nGss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio (restricted to dom 0): " << new_lnl - neutral_lnl << "\n";
                cerr << "\n\nGss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio (restricted to dom 0): " << new_lnl - neutral_lnl << "\n";
                cout << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";
                cerr << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";
                



                
		sites = parameters_to_sites(add_starting_parameters, 2);

                best_parameters = multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast_additive,
                    2
                );

                new_lnl = to_be_optimized_additive(best_parameters);

                cout << "\n\nlnl ratio (restricted to additive): " << new_lnl - neutral_lnl << "\n";
                cerr << "\n\nlnl ratio (restricted to additive): " << new_lnl - neutral_lnl << "\n";
                cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
                cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
               



                //FIX uncomment this
                /*
                sites = parameters_to_sites(dom_starting_parameters, 2);

                best_parameters = multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast_dom1,
                    2
                );

                new_lnl = to_be_optimized_pop1_dominant(best_parameters);

                cout << "\n\nGss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio (restricted to dom 1): " << new_lnl - neutral_lnl << "\n";
                cerr << "\n\nGss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio (restricted to dom 1): " << new_lnl - neutral_lnl << "\n";
                cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1,1 or " << sites[0][1]/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << 1/max(1.0, sites[0][1]) << "\n";
                cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1,1 or " << sites[0][1]/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << 1/max(1.0, sites[0][1]) << "\n";
                */
            }
                //
            //////

            //FIX replace true with something idk
            if(false){
            //////
                //Seeing if this site is actually 2 seperate sites
            

            //FIX replace both 0.2 with options.gss_out_pos[peaks[i]]
            vector<double> starting_params(6);
            starting_params[0] = options.gss_out_pos[peaks[i]] - 0.01;
            starting_params[1] = 1;
            starting_params[2] = 1;
            starting_params[3] = options.gss_out_pos[peaks[i]] + 0.01;
            starting_params[4] = 1;
            starting_params[5] = 1;


            double best_ratio = -DBL_MAX;

            vector<vector<double>> sites = parameters_to_sites(starting_params);
            vector<double> found_parameters;

            

            vector<double> best_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                two_site_bottle_necks,
                &search_sites_fast
            );
            
            
            

            sites = parameters_to_sites(best_parameters);

            double new_lnl = to_be_optimized(best_parameters);

            cout << "Gss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio: " << new_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            cerr << "site:\t" << sites[1][0] << "\t" << sites[1][1] << ",1," << sites[1][2] << " or " << sites[1][1]/max(sites[1][1], sites[1][2]) << ","<< 1/max(sites[1][1], sites[1][2]) <<"," << sites[1][2]/max(sites[1][1], sites[1][2]) << "\n";

            cout << "\n\nTesting Possibility of multiple sites\n";
            cout << "Gss lnl ratio: " << options.gss_out_llr[peaks[i]] << " new lnl ratio: " << new_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            cout << "site:\t" << sites[1][0] << "\t" << sites[1][1] << ",1," << sites[1][2] << " or " << sites[1][1]/max(sites[1][1], sites[1][2]) << ","<< 1/max(sites[1][1], sites[1][2]) <<"," << sites[1][2]/max(sites[1][1], sites[1][2]) << "\n";
            }
                //
            //////
        }
    }
    

    
    

    return ans;
}

#endif
