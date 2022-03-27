#ifndef OPTIMIZE_SELECTION
#define OPTIMIZE_SELECTION
#include <vector>
#include <math.h>



#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


selection_opt context;

void selection_opt::set_context(){
    context = *this;
}


double timer = 0;



vector<double> get_local_ancestry(vector<mat> neutral_model){
    
    map<int,vector<mat> > transition_matrix ;

    alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[0].ploidy_switch[0]], context.n_recombs, context.position, context.markov_chain_information[0].ploidy_switch[0], neutral_model ) ;

    vector<mat> interploidy_transitions;


    int sample_count = context.markov_chain_information.size();
    int site_count = context.markov_chain_information[0].alphas.size();


    //populating alphas
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        context.markov_chain_information[m].compute_forward_probabilities(  transition_matrix, interploidy_transitions) ;
    }

    //populating betas
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        context.markov_chain_information[m].compute_backward_probabilities( transition_matrix, interploidy_transitions ) ;
    }

    
    vector<double> expected_ancestry(context.markov_chain_information[0].alphas.size());

    //looping through samples
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {

        //looping through sites
        for( int i = 0; i < expected_ancestry.size(); i++){
            vec smoothed_probs = context.markov_chain_information[m].alphas[i] % context.markov_chain_information[m].betas[i] ;
            normalize( smoothed_probs ) ;

            //ploidy 1
            if(smoothed_probs.size() == 2){
                expected_ancestry[i] += smoothed_probs[0];
            }

            //ploidy 2
            if(smoothed_probs.size() == 3){
                expected_ancestry[i] += smoothed_probs[0];
                expected_ancestry[i] += smoothed_probs[1] * 0.5;
            }

            //TODO generalize ploidy here
        }

    }

    for( int i = 0; i < expected_ancestry.size(); i++){
        expected_ancestry[i] /= sample_count;
    }

    return expected_ancestry;
}





vector<mat> last_calculated_transition_matricies;

double to_be_optimized(vector<double> parameters){

    int selected_sites_count = parameters.size()/3;
    vector<double> selection_recomb_rates(selected_sites_count + 1);
    vector<vector<double>> fitnesses(selected_sites_count);

    if(selected_sites_count > 0){
        
        cerr << "Calculating likelihood for\n";
        for (uint i = 0; i < selected_sites_count; i++){
            cerr << "selection site: " << parameters[3*i + 0] << " with fitness: " << parameters[3*i + 1] << ",1," << parameters[3*i + 2] << "\n";
        }
        //Sort selected sites and fitnesses///////////////////////////
        struct Selected_pair{                                       //
            double site;
            vector<double> fitness;

            bool operator<(const Selected_pair x) const
                { return site < x.site;}
        };
        
        //Fill vector of structs
        vector<Selected_pair> selected_pairs(selected_sites_count);
        for(int i = 0; i < selected_sites_count; i++){
            Selected_pair ss;
            ss.site = parameters[i*3];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[i*3 + 2];

            selected_pairs[i] = ss;
        }
              
        //sort
        sort(selected_pairs.begin(), selected_pairs.end());         //
        //////////////////////////////////////////////////////////////


        double last = 0;
        for (uint i = 0; i < selected_pairs.size(); i++){
            selection_recomb_rates[i] = selected_pairs[i].site - last;
            last = selected_pairs[i].site;
        }

        selection_recomb_rates[selected_pairs.size()] = 1 - last;

        for(uint i = 0; i < selected_sites_count; i++){
            fitnesses[i] = selected_pairs[i].fitness;
        }
        
    }else{
        selection_recomb_rates[0] = 1;
    }


    vector<mat> transition_matrices = calculate_transition_rates(
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        context.options.m,
        context.options.generations
    );
    last_calculated_transition_matricies = transition_matrices;
    
    //sometimes i need this? sometimes i dont?
    //No i definitely need this
    cerr << '\n';
    


    map<int,vector<mat> > transition_matrix ;

    alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[0].ploidy_switch[0]], context.n_recombs, context.position, context.markov_chain_information[0].ploidy_switch[0], transition_matrices ) ;


    vector<mat> interploidy_transitions;

    
    double lnl = 0 ;
    
    
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        lnl += context.markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions) ;
    }


    cerr << "lnl = " << setprecision(15) << lnl << "\n";

    cerr << "TIME PASSED IN ONE ITERATION: " << (get_wall_time() - timer) << "\n";
    timer = get_wall_time();
    return lnl;
}


vector<mat> neutral_transition_matrices;

double to_be_optimized_only_near_sites(vector<double> parameters) {

    int selected_sites_count = parameters.size()/3;
    vector<double> selection_recomb_rates(selected_sites_count + 1);
    vector<vector<double>> fitnesses(selected_sites_count);

    if(selected_sites_count > 0){

        cerr << "Calculating likelihood for\n";
        for (uint i = 0; i < selected_sites_count; i++){
            cerr << "selection site: " << parameters[3*i + 0] << " with fitness: " << parameters[i*3 + 1] << ",1," << parameters[i*3 + 2] << "\n";
        }

        //Sort selected sites and fitnesses///////////////////////////
        struct Selected_pair{                                       //
            double site;
            vector<double> fitness;

            bool operator<(const Selected_pair x) const
                { return site < x.site;}
        };
        
        //Fill vector of structs
        vector<Selected_pair> selected_pairs(selected_sites_count);
        for(int i = 0; i < selected_sites_count; i++){
            Selected_pair ss;
            ss.site = parameters[i*3];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[3*i + 2];

            selected_pairs[i] = ss;
        }
              
        //sort
        sort(selected_pairs.begin(), selected_pairs.end());         //
        //////////////////////////////////////////////////////////////


        double last = 0;
        for (uint i = 0; i < selected_pairs.size(); i++){
            selection_recomb_rates[i] = selected_pairs[i].site - last;
            last = selected_pairs[i].site;
        }

        selection_recomb_rates[selected_pairs.size()] = 1 - last;

        for(uint i = 0; i < selected_sites_count; i++){
            fitnesses[i] = selected_pairs[i].fitness;
        }
        
    }else{
        selection_recomb_rates[0] = 1;
    }

    

    vector<mat> transition_matrices = fast_transition_rates(
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        context.options.m,
        context.options.generations
    );

    
    //sometimes i need this? sometimes i dont?
    //No i definitely need this
    cerr << '\n';
    


    map<int,vector<mat> > transition_matrix ;

    alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[0].ploidy_switch[0]], context.n_recombs, context.position, context.markov_chain_information[0].ploidy_switch[0], transition_matrices ) ;

    vector<mat> interploidy_transitions;

    double lnl = 0 ;
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        lnl += context.markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions) ;
    }



    
    vector<mat> these_neutral_transition_rates(neutral_transition_matrices.size());
    
    
    

    mat filler_matrix(2,2,fill::zeros);

    filler_matrix(0,0) = 0.5;
    filler_matrix(0,1) = 0.5;
    filler_matrix(1,0) = 0.5;
    filler_matrix(1,1) = 0.5;
    
    double sum = 0;
    for(uint i = 0; i < neutral_transition_matrices.size(); i++){
        sum += context.n_recombs[i];
        bool near_a_selected_site = false;
        for(uint j = 0; j < selected_sites_count; j++){
            double distance = sum - parameters[j*3];
            if (distance <= fast_transitions_radius_in_morgans && distance >= -fast_transitions_radius_in_morgans){
                near_a_selected_site = true;
            }
        }
        if(near_a_selected_site){
            these_neutral_transition_rates[i] = neutral_transition_matrices[i];
        }else{
            these_neutral_transition_rates[i] = filler_matrix;
        }
    }
    

    transition_matrix.clear();

    alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[0].ploidy_switch[0]], context.n_recombs, context.position, context.markov_chain_information[0].ploidy_switch[0], these_neutral_transition_rates ) ;
    
    interploidy_transitions.clear();
    
    double neutral_lnl = 0 ;
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        neutral_lnl += context.markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions) ;
    }
    


    cerr << "lnl ratio = " << setprecision(15) << lnl - neutral_lnl << "\n";

    cerr << "TIME PASSED IN ONE ITERATION: " << (get_wall_time() - timer) << "\n";
    timer = get_wall_time();
    return lnl - neutral_lnl;
}


//This function updates local_ancestries as a side effect
double to_be_optimized_ignore_far_sites(vector<double> parameters){

    int selected_sites_count = parameters.size() / 3;
    vector<double> selection_recomb_rates(selected_sites_count + 1);
    vector<vector<double>> fitnesses(selected_sites_count);

    if(selected_sites_count > 0){

        cerr << "Calculating likelihood for\n";
        for (uint i = 0; i < selected_sites_count; i++){
            cerr << "selection site: " << parameters[3*i + 0] << " with fitness: " << parameters[i*3 + 1] << ",1," << parameters[i*3 + 2] << "\n";
        }

        //Sort selected sites and fitnesses///////////////////////////
        struct Selected_pair{                                       //
            double site;
            vector<double> fitness;

            bool operator<(const Selected_pair x) const
                { return site < x.site;}
        };
        
        //Fill vector of structs
        vector<Selected_pair> selected_pairs(selected_sites_count);
        for(int i = 0; i < selected_sites_count; i++){
            Selected_pair ss;
            ss.site = parameters[i*3];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[3*i + 2];

            selected_pairs[i] = ss;
        }
              
        //sort
        sort(selected_pairs.begin(), selected_pairs.end());         //
        //////////////////////////////////////////////////////////////


        double last = 0;
        for (uint i = 0; i < selected_pairs.size(); i++){
            selection_recomb_rates[i] = selected_pairs[i].site - last;
            last = selected_pairs[i].site;
        }

        selection_recomb_rates[selected_pairs.size()] = 1 - last;

        for(uint i = 0; i < selected_sites_count; i++){
            fitnesses[i] = selected_pairs[i].fitness;
        }
        
    }else{
        selection_recomb_rates[0] = 1;
    }

    
    fast_transitions_radius_in_morgans = 0.1;
    vector<mat> transition_matrices = fast_transition_rates(
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        context.options.m,
        context.options.generations
    );
    
    
    //sometimes i need this? sometimes i dont?
    //No i definitely need this
    cerr << '\n';
    

    double sum = 0;
    for(uint i = 0; i < neutral_transition_matrices.size(); i++){
        
        sum += context.n_recombs[i];

        bool near_a_selected_site = false;

        for(uint j = 0; j < selected_sites_count; j++){

            double distance = sum - parameters[j*3];

            if ( distance <=  fast_transitions_radius_in_morgans 
              && distance >= -fast_transitions_radius_in_morgans){

                near_a_selected_site = true;
            }
        }

        if(transition_matrices[i](0,0) == 0.5){
            
            transition_matrices[i] = neutral_transition_matrices[i];
            local_ancestries[i] = context.options.m;
        }
    }
    fast_transitions_radius_in_morgans = 0.01;


    last_calculated_transition_matricies = transition_matrices;



    map<int,vector<mat> > transition_matrix ;

    alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[0].ploidy_switch[0]], context.n_recombs, context.position, context.markov_chain_information[0].ploidy_switch[0], transition_matrices ) ;


    vector<mat> interploidy_transitions;

    double lnl = 0 ;
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        lnl += context.markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions) ;
    }

    
    return lnl;
}



//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        if(i % 2 == 0){
            new_parameters.push_back(parameters[i]);
        }else{
            new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
            new_parameters.push_back(1/(1-parameters[i]/2));
        }
        
    }

    return to_be_optimized(new_parameters);
}




//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}


//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive_fast(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        if(i % 2 == 0){
            new_parameters.push_back(parameters[i]);
        }else{
            new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
            new_parameters.push_back(1/(1-parameters[i]/2));
        }
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}



vector<double> sites_to_parameters(vector<vector<double>> sites, int parameters_per_site = 3){
    if(parameters_per_site == 0){
        parameters_per_site = 3;
    }

    vector<double> parameters(sites.size()*parameters_per_site);

    for(uint i = 0; i < parameters.size(); i++){
        parameters[i] = sites[i/parameters_per_site][i%parameters_per_site];
    }

    return parameters;
}

vector<vector<double>> parameters_to_sites(vector<double> parameters, int parameters_per_site = 3){

    vector<vector<double>> sites(parameters.size()/parameters_per_site);

    for(int j = 0; j < sites.size(); j++){
        vector<double> site(parameters_per_site);

        for(int i = 0; i < parameters_per_site; i++){
            site[i] = parameters[j*parameters_per_site + i];
        }

        sites[j] = site;
    }

    return sites;
}





vector<double> search_site                       (double chrom_size, nelder_mead &opt, vector<double> site, double width, double height, double depth);


vector<double> search_sites                      (double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth);
vector<double> search_sites_fast                 (double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth);
vector<double> search_sites_fast_fix_all_but_last(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth);



vector<double> multi_level_optimization(
    double chrom_size,
    nelder_mead &opt,
    vector<vector<double>> &sites,
    vector<vector<vector<double>>> &bottle_necks,
    vector<double> (*search) (double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth),
    int parameters_per_site = 3
){
    vector<double> best_parameters;
    double best_ratio = -DBL_MAX;

    vector<double> found_parameters;
    

    for(int j = 0; j < bottle_necks.size(); j++){
        for(int k = 0; k < bottle_necks[j].size(); k++){
            for(int l = 0; l < bottle_necks[j][k][0]; l++){
                cerr << "\n SEARCH " << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                found_parameters = search(chrom_size, opt, sites, bottle_necks[j][k][1], bottle_necks[j][k][2], bottle_necks[j][k][3]);
                        
                cout << "\n Result of search" << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                cout << opt.max_value << "\n";
                for(uint j = 0; j < found_parameters.size(); j++){
                    cout << found_parameters[j] << "\n";
                }
                cout << "\n";

                if(opt.max_value > best_ratio) {
                    best_ratio = opt.max_value;
                    best_parameters = found_parameters;
                }
            }
        }
        sites = parameters_to_sites(best_parameters, parameters_per_site);
        best_ratio = -DBL_MAX;
    }

    return best_parameters;
}





vector<double> selection_opt::enact_optimization(){

    // Hyper Parameters of optimization ///////////
double nelder_mead_reflection  = 1;
double nelder_mead_contraction = 0.5;
double nelder_mead_expansion   = 2;
double nelder_mead_shrinkage   = 0.5;


    //////////////////////////////////////////////

vector<vector<vector<double>>> bottle_necks;

vector<vector<double>> shallow;

vector<double> shallow_short(4);
shallow_short[0] = 3;
shallow_short[1] = 0.03;
shallow_short[2] = 0.01;
shallow_short[3] = 20;

vector<double> shallow_medium(4);
shallow_medium[0] = 3;
shallow_medium[1] = 0.03;
shallow_medium[2] = 0.03;
shallow_medium[3] = 20;

vector<double> shallow_tall(4);
shallow_tall[0] = 3;
shallow_tall[1] = 0.03;
shallow_tall[2] = 0.05;
shallow_tall[3] = 20;

shallow.push_back(shallow_short);
shallow.push_back(shallow_medium);
shallow.push_back(shallow_tall);

vector<vector<double>> deep;

vector<double> deep_short(4);
deep_short[0] = 3;
deep_short[1] = 0.01;
deep_short[2] = 0.005;
deep_short[3] = 5;

vector<double> deep_tall(4);
deep_tall[0] = 3;
deep_tall[1] = 0.01;
deep_tall[2] = 0.01;
deep_tall[3] = 5;

deep.push_back(deep_short);
deep.push_back(deep_tall);

bottle_necks.push_back(shallow);
bottle_necks.push_back(deep);





vector<vector<vector<double>>> shallow_bottle_necks;

vector<vector<double>> shallow1;

//vector<double> shallow_short(4);
shallow_short[0] = 3;
shallow_short[1] = 0.03;
shallow_short[2] = 0.01;
shallow_short[3] = 20;

//vector<double> shallow_tall(4);
shallow_tall[0] = 3;
shallow_tall[1] = 0.03;
shallow_tall[2] = 0.05;
shallow_tall[3] = 20;

shallow1.push_back(shallow_short);
shallow1.push_back(shallow_tall);

vector<vector<double>> shallow2;

//vector<double> deep_short(4);
deep_short[0] = 3;
deep_short[1] = 0.01;
deep_short[2] = 0.005;
deep_short[3] = 5;

shallow2.push_back(deep_short);

shallow_bottle_necks.push_back(shallow1);
shallow_bottle_necks.push_back(shallow2);

    // Create optimizer
    nelder_mead optimizer(
        nelder_mead_reflection,
        nelder_mead_contraction,
        nelder_mead_expansion,
        nelder_mead_shrinkage
    );

    //set_context();
    context = *this;
    
    
    double chrom_size = 0;
    for(uint i = 0; i < n_recombs.size(); i++){
        chrom_size += n_recombs[i];
    }

    //TESTER CODE
    if(false){

        vector<double> poss(3);
        poss[0] = 0.1;
        poss[1] = 1;
        poss[2] = 1.0204081632;
        double poss_lnl = to_be_optimized(poss);


        vector<double> real(3);
        real[0] = 0.1;
        real[1] = 1;
        real[2] = 1.0204081632;
        double real_lnl = to_be_optimized(real);
        cerr << "\n";

        cerr << "Poss point lnl: \t" << poss_lnl << "\n";
        cerr << "Real point lnl: \t" << real_lnl << "\n";
    }

    // Calculating lnl for neutral model
    vector<double> empty(0);
    double neutral_lnl = to_be_optimized(empty);
    

    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;




    



    
    vector<double> real(3);
    real[0] = 0.2;
    real[1] = 0.995;
    real[2] = 0.99;
    double real_lnl = to_be_optimized(real);
    



    //iterative selected site adding
    vector<vector<double>> sites;
    vector<bool> site_has_been_deep_searched; 

    double last_lnl = neutral_lnl;
    cout << "Neutral lnl\t" << setprecision(15) << neutral_lnl << "\n";
    vector<double> data_ancestry;
    vector<double> expected_ancestry;

    //STANDARD: for now, only adds 5 sites
    while(sites.size() < 5){

        data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
        expected_ancestry = local_ancestries;

        vector<double> smoothed_data_ancestry(data_ancestry.size());

        double largest_deviation = 0;
        int largest_deviator = 0;
        for(int i = 1; i < n_recombs.size() - 1; i++){

            double total = 0;
            double count = 0;
            for(int j = 0; morgan_position[i] - morgan_position[j] < 0.001 && j >= 0; j--){
                total += data_ancestry[j];
                count ++;
            }
            for(int j = 1; morgan_position[j] - morgan_position[i] < 0.001 && j < n_recombs.size(); j++){
                total += data_ancestry[j];
                count ++;
            }
            smoothed_data_ancestry[i] = total/count;


            if (abs(smoothed_data_ancestry[i] - expected_ancestry[i]) > largest_deviation){
                largest_deviation = abs(smoothed_data_ancestry[i] - expected_ancestry[i]);
                largest_deviator = i;
                //cerr << "largestdeviator" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] <<"\n";
            }
            if(i%100 == 0){
                cerr << setprecision(5) << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] << "\n";
            }
        }

        vector<double> new_site(3);
        new_site[0] = morgan_position[largest_deviator];
        new_site[1] = 1;
        new_site[2] = 1;
        

        cerr << "placing at :\n";
        cerr << new_site[0] << "\n";
        cout << "Placing new site at: " << new_site[0] << "\n";

        vector<double> best_parameters;
        double best_ratio = -DBL_MAX;

        vector<double> found_parameters;

        
        //generalized bottlenecking
        if(true){

            //Checking if new site is near an old site
            bool new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = false;
            int site_its_close_to = 0;
            for(int s = 0; s < sites.size(); s++){
                if( abs(sites[s][0] - new_site[0]) < 0.05 && !site_has_been_deep_searched[s]){
                    new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = true;
                    site_its_close_to = s;
                }
            }

            cout << "close: " << new_site_is_close_to_existing_site_that_hasnt_been_deep_searched << "\n";

            if(new_site_is_close_to_existing_site_that_hasnt_been_deep_searched){
                //pop close one from sites and place it at the end
                // then do a bottle_necks search with
                // search sites fast fix all but last

                //Remove site its close to and place it at the end
                cout << "its close\n";

                vector<double> near_site = sites[site_its_close_to];
                sites.erase(sites.begin() + site_its_close_to);
                sites.push_back(near_site);

                multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast_fix_all_but_last
                );

                site_has_been_deep_searched[site_has_been_deep_searched.size() - 1] = true;

            }else{
                
                cout << "its not close\n";

                sites.push_back(new_site);
                site_has_been_deep_searched.push_back(false);
                // do a bottle_necks search with
                // search sites fast fix all but last

                multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    shallow_bottle_necks,
                    &search_sites_fast_fix_all_but_last
                );
            }
            
            /*
            multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );
            */


            // i need a fine tune afterward when theres more than one site

            /*
            for(int j = 0; j < bottle_necks.size(); j++){
                for(int k = 0; k < bottle_necks[j].size(); k++){
                    for(int l = 0; l < bottle_necks[j][k][0]; l++){
                        cerr << "\n SEARCH " << j << "/" << bottle_necks.size() << " " << k << "/" << bottle_necks[j].size() << " " << l << "/" << bottle_necks[j][k][0] << "\n";
                        found_parameters = search_sites_fast_fix_all_but_last(chrom_size, optimizer, sites, bottle_necks[j][k][1], bottle_necks[j][k][2], bottle_necks[j][k][3]);
                        
                        cout << "\n Result of search " << j << "/" << bottle_necks.size() << " " << k << "/" << bottle_necks[j].size() << " " << l << "/" << bottle_necks[j][k][0] << "\n";
                        cout << optimizer.max_value << "\n";
                        for(uint j = 0; j < found_parameters.size(); j++){
                            cout << found_parameters[j] << "\n";
                        }
                        cout << "\n";

                        if(optimizer.max_value > best_ratio) {
                            best_ratio = optimizer.max_value;
                            best_parameters = found_parameters;
                        }
                    }
                }
                sites = parameters_to_sites(best_parameters);
                best_ratio = -DBL_MAX;
            }
            */
        }else{

        }
        //end of else


        cerr << "\n\nBest sites so far:\n";
        for(int k = 0; k < sites.size(); k++){

            cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
        }
        cerr << "\n";

        cout << "\nnew site added:\t" << sites[sites.size() - 1][0]
        << "\n" << sites[sites.size() - 1][1]
        << ",1," << sites[sites.size() - 1][2]
        << " or\n" << sites[sites.size() - 1][1]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2])
        << ","<< 1/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2])
        << "," << sites[sites.size() - 1][2]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2]) << "\n";

        best_parameters = sites_to_parameters(sites);
        double new_lnl = to_be_optimized(best_parameters);

        cout <<"lnl after new site\t" << setprecision(15) << new_lnl << "\n";
        

        cerr << "lnl: " << new_lnl << "\n";

        if(new_lnl - last_lnl < 1){
            sites.pop_back();
            break;
        }
        last_lnl = new_lnl;

        //expected_ancestry = local_ancestries;

        //exit(0);
    }

    cerr << "\n\nTERMINATED SEARCH\n\n\n";
    cout << "TERMINATED SEARCH\n";

    // I have to re-check sites by themselves to see if they really add to the lnl
        

    cerr << "Found sites:\n\n";
    for(int k = 0; k < sites.size(); k++){
        cout << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
        cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
    }
    
    
    
    return sites_to_parameters(sites);
}





vector<double> search_sites(
    double chrom_size, 
    nelder_mead &opt,
    vector<vector<double>> sites,
    double width,
    double height,
    double depth
)
{
    vector<double> center_point(sites.size()*3);
    vector<double> scales(sites.size()*3);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*3+0] = sites[i][0];
        center_point[i*3+1] = sites[i][1];
        center_point[i*3+2] = sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }

    opt.populate_points(sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }
    
    opt.calculate_points(&to_be_optimized);

    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized);
        
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    return opt.points[opt.max_index];
}



vector<double> search_sites_fast(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*3);
    vector<double> scales(sites.size()*3);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*3+0] = sites[i][0];
        center_point[i*3+1] = sites[i][1];
        center_point[i*3+2] = sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }

    opt.populate_points(sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.calculate_points(&to_be_optimized_only_near_sites);

    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_only_near_sites);
        
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    return opt.points[opt.max_index];
}


vector<double> search_sites_fast_fix_all_but_last(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<vector<double>> new_sites;

    vector<bool> in_new_sites(sites.size());

    for(uint i = 0; i < sites.size(); i++){

        in_new_sites[i] = false;

        double distance_from_last = abs(sites[i][0] - sites[sites.size() - 1][0]);

        if (distance_from_last <= fast_transitions_radius_in_morgans) {
            new_sites.push_back(sites[i]);
            in_new_sites[i] = true;
        }

    }

    vector<double> center_point(new_sites.size()*3);
    vector<double> scales(new_sites.size()*3);

    for(uint i = 0; i < new_sites.size(); i++){
        center_point[i*3+0] = new_sites[i][0];
        center_point[i*3+1] = new_sites[i][1];
        center_point[i*3+2] = new_sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }


    opt.populate_points(new_sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.calculate_points(&to_be_optimized_only_near_sites);

    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_only_near_sites);
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    vector<double> ans(sites.size()*3);
    
    int counter = 0;

    for(uint i = 0; i < sites.size(); i++){
        if(in_new_sites[i]){
            ans[i*3 + 0] = opt.points[opt.max_index][counter*3 + 0];
            ans[i*3 + 1] = opt.points[opt.max_index][counter*3 + 1];
            ans[i*3 + 2] = opt.points[opt.max_index][counter*3 + 2];

            counter++;
        }else{ 
            ans[i*3 + 0] = sites[i][0];
            ans[i*3 + 1] = sites[i][1];
            ans[i*3 + 2] = sites[i][2];
        }
    }

    return ans;
}






//BAD
vector<double> search_sites_fast_dom0(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.calculate_points(&to_be_optimized_pop0_dominant_fast);
    int iterations = 0;
    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4 && iterations < 16){
        opt.iterate(&to_be_optimized_pop0_dominant_fast);
        iterations++;
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    return opt.points[opt.max_index];
}

vector<double> search_sites_fast_dom1(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.calculate_points(&to_be_optimized_pop1_dominant_fast);
    int iterations = 0;
    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4 && iterations < 16){
        opt.iterate(&to_be_optimized_pop1_dominant_fast);
        iterations++;
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    return opt.points[opt.max_index];
}


vector<double> search_sites_fast_additive(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.calculate_points(&to_be_optimized_additive_fast);
    int iterations = 0;
    opt.repeated_shrinkages = 0;
    while(opt.max_value - opt.min_value > depth && opt.repeated_shrinkages < 4 && iterations < 16){
        opt.iterate(&to_be_optimized_additive_fast);
        iterations++;
        cerr << "INFO simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
    }

    return opt.points[opt.max_index];
}
#endif