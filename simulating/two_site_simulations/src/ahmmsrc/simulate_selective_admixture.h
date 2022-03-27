#ifndef SIMULATE_SELECTIVE_ADMIXTURE_H
#define SIMULATE_SELECTIVE_ADMIXTURE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

#include <pthread.h>

#define EPSILON 0.0000001


using namespace arma;
using namespace std;

double gwall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int gamete_to_index(vector<bool> gamete){
    int ans = 0;
    for(unsigned int i = 0; i < gamete.size(); i++)
        ans += gamete[i] * (1 << i);
    
    return ans;
}

vector<bool> index_to_gamete(int index, int sites){
    vector<bool> ans(sites);
    for(int s = 0; s < sites; s++)
        ans[s] = (index >> s) % 2;
    return ans;
}





void HtoD(Row<double> *H, Row<double> *D){
    for (uint i = 0 ; i < H->n_elem; i++){
        for (uint j = 0; j < H->n_elem; j++){
            (*D)(i*H->n_elem + j) = (*H)(i) * (*H)(j);
        }
    }
}


void normH(Row<double> *H){
    double total = 0;
    for (uint i = 0; i < H->n_elem; i++){
        total += (*H)(i);
    }
    for (uint i = 0; i < H->n_elem; i++){
        (*H)(i) = (*H)(i)/total;
    }
}



vector<double> adjacent_transition_rate(vector<double> recomb_rates, vector<vector<double>> fitnesses, double m,
    double n_site1, double n_site2, int generations){
    
    
    // inserting netural sites //////////////////////////////////////////////////////////////
    int sites = fitnesses.size() + 2;                                                      //
    int n_sites[2];

    vector<double> neutral_fit(3);
    neutral_fit[0] = 1;
    neutral_fit[1] = 1;
    neutral_fit[2] = 1;
    
    vector<double> new_recomb_rates = recomb_rates;
    vector<vector<double>> new_fitnesses = fitnesses;

    
    //insert neutral sites into vectors
    bool inserted_first = false;
    double recomb_sum = 0;
    
    for(unsigned int i = 0; i < recomb_rates.size(); i++){
        recomb_sum += recomb_rates[i];

        if(!inserted_first){
            if(n_site1 - EPSILON <= recomb_sum){
                auto fit = new_fitnesses.begin();
                new_fitnesses.insert(fit + i, neutral_fit);
                n_sites[0] = i;

                auto rit = new_recomb_rates.begin();
                new_recomb_rates.insert(rit + i + 1, recomb_sum - n_site1);
                new_recomb_rates[i] = n_site1 - (recomb_sum - recomb_rates[i]);
                

                inserted_first = true;
            }
        }
        
        if(inserted_first){
            if(n_site2 - EPSILON <= recomb_sum){
                auto fit = new_fitnesses.begin();
                new_fitnesses.insert(fit + i + 1, neutral_fit);
                n_sites[1] = i + 1;

                auto rit = new_recomb_rates.begin();
                new_recomb_rates.insert(rit + i + 2, recomb_sum - n_site2);
                new_recomb_rates[i + 1] =  n_site2 - (recomb_sum - new_recomb_rates[i+1]);

                break;
            }
        }
    }

    fitnesses = new_fitnesses;
    recomb_rates = new_recomb_rates;                                                       //
    /////////////////////////////////////////////////////////////////////////////////////////
    
    

    int haploids = pow(2,sites);
    
    Row<double> H(haploids, fill::zeros);
    H(0) = m;
    H(haploids - 1) = 1 - m;

    mat M(haploids*haploids, haploids, fill::zeros);

    Row<double> D(haploids*haploids, fill::zeros);

    // populate M matrix //////////////////////////////////////////////////////////////
    //loop through all diploids                                                      //
    
    for (int i = 0; i < haploids; i++){
        for (int j = 0; j < haploids; j++){

            int di = i*haploids + j;
            vector<bool> chr1 = index_to_gamete(i, sites);
            vector<bool> chr2 = index_to_gamete(j, sites);


            double fitness = 1;

            for(int s = 0; s < sites; s++){
                if     (chr1[s] == 0 && chr2[s] == 0)
                    fitness *= fitnesses[s][0];
                else if(chr1[s] == 1 && chr2[s] == 1)
                    fitness *= fitnesses[s][2];
                else
                    fitness *= fitnesses[s][1];
            }

            //calculating possible recombinants
            vector<bool> recomb_chr1 = chr1;
            vector<bool> recomb_chr2 = chr2;

            // split before first site
            M(di, gamete_to_index(recomb_chr1)) += recomb_rates[0]/2 * fitness;
            M(di, gamete_to_index(recomb_chr2)) += recomb_rates[0]/2 * fitness;

            // split directly after site s
            for(int s = 0; s < sites; s++){
                recomb_chr1[s] = chr2[s];
                recomb_chr2[s] = chr1[s];

                M(di, gamete_to_index(recomb_chr1)) += recomb_rates[s+1]/2 * fitness;
                M(di, gamete_to_index(recomb_chr2)) += recomb_rates[s+1]/2 * fitness;

            }

        }
    }                                                                                //
    ///////////////////////////////////////////////////////////////////////////////////
    
    
    // Run through generations /////////////
    for(int g = 0; g < generations; g++) {//
        HtoD(&H, &D);
        H = D * M;
        normH(&H);
    }                                     //
    ////////////////////////////////////////

    
    // Calculating transition rates ///////////////////////////////
    vector<double> transitions(4);                               //
    
    for(int i = 0; i < haploids; i++){
        int anc_of_site_1 = index_to_gamete(i,sites)[n_sites[0]];
        int anc_of_site_2 = index_to_gamete(i,sites)[n_sites[1]];

        if(anc_of_site_1 == 0 && anc_of_site_2 == 0)
            transitions[0] += H(i);
        
        if(anc_of_site_1 == 0 && anc_of_site_2 == 1)
            transitions[1] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 0)
            transitions[2] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 1)
            transitions[3] += H(i);

    }                                                            //
    ///////////////////////////////////////////////////////////////
    

    return transitions;
}




vector<double> local_ancestries;

int complete_transition_rates_number_of_threads = 1;



vector<double> *global_neutral_sites;
vector<double> *global_recomb_rates_around_selected_sites;
vector<vector<double>> *global_fitnesses;
double global_m;
int global_generations;
vector<mat> *global_transition_matrices;

void *single_window_process(void *threadid){

    long t = (long)threadid;
    
    for(uint i = t; i < (*global_transition_matrices).size(); i+= complete_transition_rates_number_of_threads){

        vector<double> trans = adjacent_transition_rate(*global_recomb_rates_around_selected_sites, *global_fitnesses, global_m, (*global_neutral_sites)[i], (*global_neutral_sites)[i + 1], global_generations);

        mat transition_matrix(2,2,fill::zeros);

        transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
        transition_matrix(0,1) = 1 - transition_matrix(0,0);
        transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
        transition_matrix(1,1) = 1 - transition_matrix(1,0);

        (*global_transition_matrices)[i] = transition_matrix;
        local_ancestries[i] = trans[0] + trans[1];
    }

    pthread_exit(NULL);
    return NULL;
}




vector<mat> calculate_transition_rates(
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations
){
    vector<mat> transition_matrices(recomb_rates_around_neutral_sites.size());
    local_ancestries.clear();
    local_ancestries.resize(recomb_rates_around_neutral_sites.size());

    //NOTE im changing the relationship between neutral_sites and recomb_rates_around_neutral_sites
    // im now treating recomb_rates_around_neutral_sites as if it only contains regions delimited by two neutral sites
    // thus neutral_sites.size = recomb_rates_around_neutral_sites.size + 1

    // Creating site vectors //////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size()  + 1);//
    vector<double> selected_sites(recomb_rates_around_selected_sites.size() - 1);

    double sum = 0;
    neutral_sites[0] = 0; //extra
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++){
        sum += recomb_rates_around_selected_sites[i];
        selected_sites[i] = sum;
    }                                                                            //
    ///////////////////////////////////////////////////////////////////////////////


    // Populating global variables with local variables //////////////////////////////
    global_neutral_sites = &neutral_sites;                                          //
    global_recomb_rates_around_selected_sites = &recomb_rates_around_selected_sites;
    global_fitnesses = &fitnesses;
    global_m = m;
    global_generations = generations;
    global_transition_matrices = &transition_matrices;                              //
    //////////////////////////////////////////////////////////////////////////////////


    // Creating threads to deal with independent adjacent neutral regions //////////////////
    vector<pthread_t> threads(complete_transition_rates_number_of_threads);               //

    

    for(long t = 0; t < complete_transition_rates_number_of_threads; t++){
        int rc = pthread_create(&threads[t], NULL, single_window_process, (void *)t);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    
    //wait for all to finish by joining them
    for (int t = 0; t < complete_transition_rates_number_of_threads; t++) {
        pthread_join(threads[t], NULL);
    }                                                                                     //
    ////////////////////////////////////////////////////////////////////////////////////////

    
    
    return transition_matrices;
}






void alt_create_transition_matrix ( map<int,vector<mat> > &transition_matrix , vector<vector< map< vector<transition_information>, double > > > &transition_info, vector<double> &recombination_rate, vector<int> &positions, double &number_chromosomes, vector<mat> &transition_matrices);


//indicie_bounds[i] corrosponds to a thread i, which will find transitions for
// all grid sites in the range [indicie_bounds[i][0], indicie_bounds[i][1])
vector<vector<int>> indicie_bounds;

selection_opt grid_context;
vector<vector<double>> *global_lnl_ratios;
vector<double> *global_grid_sites;
double grid_transitions_radius_in_morgans = 0.05;
vector<vector<double>> grid_fitnesses;
vector<vector<markov_chain>> thread_specific_markov_chain_information;



void *single_grid_process(void *threadid){

    long t = (long)threadid;
    for(uint i = indicie_bounds[t][0]; i < indicie_bounds[t][1]; i++){
        
        //neutral sites index bounds
        uint bottom_trans_matrix_index = 0;
        uint top_trans_matrix_index = 0;
        
        
        uint starting_point = i * (*global_neutral_sites).size()/(*global_grid_sites).size();
        
        //setting bottom_trans_matrix_index to be the lowest  neutral site which is within grid_transitions_radius_in_morgans away from grid site
        //setting    top_trans_matrix_index to be the highest neutral site which is within grid_transitions_radius_in_morgans away from grid site //////////////////////////////////////////////
        if((*global_grid_sites)[i] - grid_transitions_radius_in_morgans < (*global_neutral_sites)[starting_point]){                                                                           //
            for(
                bottom_trans_matrix_index = starting_point;
                (*global_grid_sites)[i] - grid_transitions_radius_in_morgans < (*global_neutral_sites)[bottom_trans_matrix_index] && bottom_trans_matrix_index > 0;
                bottom_trans_matrix_index--
            ){}
            for(
                top_trans_matrix_index = bottom_trans_matrix_index;
                (*global_grid_sites)[i] + grid_transitions_radius_in_morgans > (*global_neutral_sites)[top_trans_matrix_index] && top_trans_matrix_index < (*global_neutral_sites).size() - 1;
                top_trans_matrix_index++
            ){}
        }
        else if((*global_grid_sites)[i] + grid_transitions_radius_in_morgans >= (*global_neutral_sites)[starting_point]){

            for(
                top_trans_matrix_index = starting_point;
                (*global_grid_sites)[i] + grid_transitions_radius_in_morgans > (*global_neutral_sites)[top_trans_matrix_index] && top_trans_matrix_index < (*global_neutral_sites).size() - 1;
                top_trans_matrix_index++
            ){}
            for(
                bottom_trans_matrix_index = top_trans_matrix_index;
                (*global_grid_sites)[i] - grid_transitions_radius_in_morgans < (*global_neutral_sites)[bottom_trans_matrix_index] && bottom_trans_matrix_index > 0;
                bottom_trans_matrix_index--
            ){}
        }                                                                                                                                                                                     //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        


        vector<double> recomb_around_selected_site(2);
        recomb_around_selected_site[0] = (*global_grid_sites)[i];
        recomb_around_selected_site[1] = (1 - (*global_grid_sites)[i]);

        
        for(uint j = 0; j < grid_fitnesses.size(); j++){
            
            vector<vector<double>> lone_fitness(1);
            lone_fitness[0] = grid_fitnesses[j];

            
            vector<mat> truncated_transition_matricies((*global_neutral_sites).size() - 1);
            
            for(uint k = 0; k < truncated_transition_matricies.size(); k++){
                
                mat transition_matrix(2,2,fill::zeros);
                
                if(k < bottom_trans_matrix_index || k > top_trans_matrix_index){
                    transition_matrix(0,0) = 0.5;
                    transition_matrix(0,1) = 0.5;
                    transition_matrix(1,0) = 0.5;
                    transition_matrix(1,1) = 0.5;
                }else{
                    
                    vector<double> trans = adjacent_transition_rate(recomb_around_selected_site, lone_fitness, global_m, (*global_neutral_sites)[k], (*global_neutral_sites)[k + 1], global_generations);
                    
                    transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
                    transition_matrix(0,1) = 1 - transition_matrix(0,0);
                    transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
                    transition_matrix(1,1) = 1 - transition_matrix(1,0);
                }
                
                truncated_transition_matricies[k] = transition_matrix;
                
            }
            
            
            map<int,vector<mat> > transition_matrix ;
            
            alt_create_transition_matrix( transition_matrix, grid_context.transition_matrix_information[grid_context.markov_chain_information[0].ploidy_switch[0]], grid_context.n_recombs, grid_context.position, grid_context.markov_chain_information[0].ploidy_switch[0], truncated_transition_matricies ) ;
            

            
            vector<mat> interploidy_transitions;
            
    
            double lnl = 0 ;

            
            
            for ( int m = 0 ; m < grid_context.markov_chain_information.size() ; m ++ ) {
                lnl += grid_context.markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions, bottom_trans_matrix_index, top_trans_matrix_index - 1) ;
            }
            
            if(j){
                (*global_lnl_ratios)[i][j] = lnl - (*global_lnl_ratios)[i][0];
            }else{
                (*global_lnl_ratios)[i][0] = lnl;
            }
            
        }
        
        (*global_lnl_ratios)[i][0] = 0;

        cerr << "\n lnl ratios from site " << (*global_grid_sites)[i] << "\n";
        for(uint j = 0; j < (*global_lnl_ratios)[i].size(); j++){
            cerr << (*global_lnl_ratios)[i][j] << " ";
        }
        cerr << "\n";

    }
    

    pthread_exit(NULL);
    return NULL;
}

vector<vector<double>> grid_search(
  vector<double> &recomb_rates_around_neutral_sites,
  double m,
  int generations,
  selection_opt &context
){

    grid_fitnesses.resize(0);

    vector<double> point(3);
    point[0] = 1;
    point[1] = 1;
    point[2] = 1;

    grid_fitnesses.push_back(point);

    for(int i = -1; i < 2; i+=2){
        for(int j = -1; j < 2; j+=2){
            point[0] = 1 + 0.02*i;
            point[1] = 1;
            point[2] = 1 + 0.02*j;
            grid_fitnesses.push_back(point);
        }
    }
    
    // Creating site vectors /////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size() + 1);//

    double sum = 0;
    neutral_sites[0] = 0;
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    int grid_divisions = floor(sum * 200);                                      //
    //////////////////////////////////////////////////////////////////////////////

    
    vector<double> grid_sites(grid_divisions);
    for(int i = 0; i < grid_divisions; i++){
        grid_sites[i] = ((i + 0.5)/grid_divisions * sum);
    }

    vector<vector<double>> lnl_ratios(grid_divisions);
    for(uint i = 0; i < grid_divisions; i++){
        lnl_ratios[i].resize(grid_fitnesses.size());
    }



    // Populating global variables with local variables //
    global_neutral_sites = &neutral_sites;              //
    global_m = m;
    global_generations = generations;
    global_grid_sites = &grid_sites;
    grid_context = context;
    global_lnl_ratios = &lnl_ratios;                    //
    //////////////////////////////////////////////////////


    // Assigning different grid points to different threads /////////////////////////
    int number_of_threads = 5;                                                     //
    

    vector<pthread_t> threads(number_of_threads);
    for(uint i = 0; i < indicie_bounds.size();i++){
        indicie_bounds[i].clear();
    }

    indicie_bounds.resize(number_of_threads);

    for(int i = 0; i < number_of_threads; i++){
        indicie_bounds[i].push_back((grid_sites.size()*i)/number_of_threads);
        indicie_bounds[i].push_back((grid_sites.size()*(i+1))/number_of_threads);
    }
    

    for(long t = 0; t < number_of_threads; t++){
        int rc = pthread_create(&threads[t], NULL, single_grid_process, (void *)t);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }
    
    //wait for all to finish by joining them
    for (int t = 0; t < number_of_threads; t++) {
        pthread_join(threads[t], NULL);
    }                                                                              //
    /////////////////////////////////////////////////////////////////////////////////
    
    return lnl_ratios;
}













//This uses the same global variables as calculate transition rates
int fast_transition_rates_number_of_threads = 1;

double fast_transitions_radius_in_morgans = 0.05;

vector<double> *global_selected_sites;


void *single_fast_window_process(void *threadid){

    long t = (long)threadid;
    
    for(uint i = t; i < (*global_transition_matrices).size(); i+= fast_transition_rates_number_of_threads){
        
        vector<double> pertinent_selected_sites;
        vector<vector<double>> pertinent_fitnesses;

        for(uint j = 0; j < (*global_selected_sites).size(); j++){
            double distance = (*global_neutral_sites)[i] - (*global_selected_sites)[j];
            if(distance <= fast_transitions_radius_in_morgans && distance >= -fast_transitions_radius_in_morgans){
                pertinent_selected_sites.push_back((*global_selected_sites)[j]);
                pertinent_fitnesses.push_back((*global_fitnesses)[j]);
            }
        }

        if(pertinent_selected_sites.size() >= 1){

            vector<double> pertinent_recomb_rates_around_selected_sites(pertinent_selected_sites.size() + 1);
            pertinent_recomb_rates_around_selected_sites[0] = pertinent_selected_sites[0];

            for(uint j = 1; j < pertinent_selected_sites.size(); j++){
                pertinent_recomb_rates_around_selected_sites[j] = pertinent_selected_sites[j] - pertinent_selected_sites[j - 1];
            }
            pertinent_recomb_rates_around_selected_sites[pertinent_selected_sites.size()] = 1 - pertinent_selected_sites[pertinent_selected_sites.size() - 1];


            vector<double> trans = adjacent_transition_rate(pertinent_recomb_rates_around_selected_sites, pertinent_fitnesses, global_m, (*global_neutral_sites)[i], (*global_neutral_sites)[i + 1], global_generations);

            mat transition_matrix(2,2,fill::zeros);

            transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
            transition_matrix(0,1) = 1 - transition_matrix(0,0);
            transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
            transition_matrix(1,1) = 1 - transition_matrix(1,0);

            (*global_transition_matrices)[i] = transition_matrix;
            local_ancestries[i] = trans[0] + trans[1];
        }else{
            mat transition_matrix(2,2,fill::zeros);

            transition_matrix(0,0) = 0.5;
            transition_matrix(0,1) = 0.5;
            transition_matrix(1,0) = 0.5;
            transition_matrix(1,1) = 0.5;

            (*global_transition_matrices)[i] = transition_matrix;
        }
    }
    
    pthread_exit(NULL);
    return NULL;
}



vector<mat> fast_transition_rates(
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations
){
    vector<mat> transition_matrices(recomb_rates_around_neutral_sites.size());
    local_ancestries.clear();
    local_ancestries.resize(recomb_rates_around_neutral_sites.size());

    //NOTE im changing the relationship between neutral_sites and recomb_rates_around_neutral_sites
    // im now treating recomb_rates_around_neutral_sites as if it only contains regions delimited by two neutral sites
    // thus neutral_sites.size = recomb_rates_around_neutral_sites.size + 1

    // Creating site vectors //////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size()  + 1);//
    vector<double> selected_sites(recomb_rates_around_selected_sites.size() - 1);

    double sum = 0;
    neutral_sites[0] = 0; //extra
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++){
        sum += recomb_rates_around_selected_sites[i];
        selected_sites[i] = sum;
    }                                                                            //
    ///////////////////////////////////////////////////////////////////////////////


    // Populating global variables with local variables //////////////////////////////
    global_neutral_sites = &neutral_sites;                                          //
    global_selected_sites = &selected_sites;
    global_recomb_rates_around_selected_sites = &recomb_rates_around_selected_sites;
    global_fitnesses = &fitnesses;
    global_m = m;
    global_generations = generations;
    global_transition_matrices = &transition_matrices;                              //
    //////////////////////////////////////////////////////////////////////////////////



    // Creating threads to deal with independent adjacent neutral regions //////////////////
    
    vector<pthread_t> threads(fast_transition_rates_number_of_threads);

    for(long t = 0; t < fast_transition_rates_number_of_threads; t++){
        int rc = pthread_create(&threads[t], NULL, single_fast_window_process, (void *)t);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    
    //wait for all to finish by joining them
    for (int t = 0; t < fast_transition_rates_number_of_threads; t++) {
        pthread_join(threads[t], NULL);
    }                                                                                     //
    ////////////////////////////////////////////////////////////////////////////////////////

    
    
    return transition_matrices;
}















//// create transition matrix for a given admixture model
void alt_create_transition_matrix ( map<int,vector<mat> > &transition_matrix , vector<vector< map< vector<transition_information>, double > > > &transition_info, vector<double> &recombination_rate, vector<int> &positions, double &number_chromosomes, vector<mat> &transition_matrices) {
    


    /// check if we already computed this for this sample ploidy
    if ( transition_matrix.find( number_chromosomes ) != transition_matrix.end() ) {
        return ;
    }
    
    /// else, have to create entire matrix
    /// first create data object of approporate size
    transition_matrix[number_chromosomes].resize(recombination_rate.size()) ;



    //// iterate across all positions and compute transition matrixes

    //NOTE changed p starting point from 1 to 0

    for ( int p = 0 ; p < recombination_rate.size(); p ++ ) {
        
        mat segment_transitions = transition_matrices[p];
        
        

        /// population transitions by summing across all routes

        transition_matrix[number_chromosomes][p] = mat(transition_info.size(),transition_info.size(), fill::zeros);

        
        for ( int i = 0 ; i < transition_info.size() ; i ++ ) {
            for ( int j = 0 ; j < transition_info[i].size() ; j ++ ) {
                for ( std::map<vector<transition_information>,double>::iterator t = transition_info[i][j].begin() ; t != transition_info[i][j].end() ; ++ t ) {
                    
                    double prob_t = 1 ;
                    for ( int r = 0 ; r < t->first.size() ; r ++ ) {
                        prob_t *= pow( segment_transitions(t->first[r].start_state,t->first[r].end_state), t->first[r].transition_count ) ;
                    }

                    
                    transition_matrix[number_chromosomes][p](j,i) += prob_t * t->second ;
                }
            }
        }
    }
}



















//// create all transition rates between ancestry types for a single chromosome
mat alt_create_transition_rates ( vector<pulse> admixture_pulses, double n, vector<double> ancestry_proportion ) {
    
    /// determine ancestry proportions based on the fraction of remainder
    for ( int p = 0 ; p < admixture_pulses.size() ; p ++ ) {
        admixture_pulses[p].proportion = admixture_pulses[p].fraction_of_remainder * ancestry_proportion[admixture_pulses[p].type] ;
        ancestry_proportion[admixture_pulses[p].type] -= admixture_pulses[p].proportion ;
    }
    
    //// sort by time so oldest events are at the top
    sort( admixture_pulses.begin(), admixture_pulses.end() ) ;
    
    /// all ancestry proprionts to this point are in units of final ancestry
    /// however before the pulses that came before, there was more
    double ancestry_accounted = 0 ;
    for ( int p = admixture_pulses.size() - 1 ; p > 0 ; p -- ) {
        double accounted_here = admixture_pulses[p].proportion ;
        admixture_pulses[p].proportion = admixture_pulses[p].proportion/( 1 - ancestry_accounted ) ;
        ancestry_accounted += accounted_here ;
    }

    /// matrix to hold transition rates
    mat transition_rates( admixture_pulses.size(), admixture_pulses.size(), fill::zeros ) ;
    
    /// iterate through all states
    for ( int s1 = 0 ; s1 < admixture_pulses.size() ; s1 ++ ) {
        for ( int s2 = 0 ; s2 < admixture_pulses.size() ; s2 ++ ) {
            
            //// self transition rates are just going to be 1-all others
            if ( s2 == s1 ) continue ;
            
            /// rates calculated one way if greater than
            if ( s1 > s2 ) {
                for ( int t = s1 ; t < admixture_pulses.size() ; t ++ ) {
                    
                    /// basic recombination rates in each epoch
                    double rate ;
                    if ( t != admixture_pulses.size() - 1 ) {
                        rate = n * ( 1 - exp( (admixture_pulses[t+1].time-admixture_pulses[t].time)/n ) ) ;
                    }
                    else {
                        rate = n * ( 1 - exp( (-1*admixture_pulses[t].time)/n ) ) ;
                    }
                    
                    /// probability of no recombination in prior epochs
                    /// will skip epoch of s1 since that's in the above equation by definition
                    for ( int t2 = t - 1 ; t2 > s1 - 1 ; t2 -- ) {
                        rate *= exp((admixture_pulses[t2+1].time-admixture_pulses[t2].time)/n) ;
                    }
                    
                    /// now probability of not selecting other acnestry types
                    for ( int a = t ; a > s2 ; a -- ) {
                        rate *= ( 1 - admixture_pulses[a].proportion ) ;
                    }
                    
                    /// now select the correct ancestry type
                    if ( s2 != 0 ) {
                        rate *= admixture_pulses[s2].proportion ;
                    }
                    
                    transition_rates(admixture_pulses[s2].entry_order,admixture_pulses[s1].entry_order) += rate ;
                }
            }
            else {
                for ( int t = s2 ; t < admixture_pulses.size() ; t ++ ) {
                    
                    /// basic recombination rates in each epoch
                    double rate ;
                    if ( t != admixture_pulses.size() - 1 ) {
                        rate = n * ( 1 - exp( (admixture_pulses[t+1].time-admixture_pulses[t].time)/n ) ) ;
                    }
                    else {
                        rate = n * ( 1 - exp( (-1*admixture_pulses[t].time)/n ) ) ;
                    }
                    
                    /// probability of no recombination in prior epochs
                    /// will skip epoch of s1 since that's in the above equation by definition
                    for ( int t2 = t - 1 ; t2 > s2 - 1 ; t2 -- ) {
                        rate *= exp((admixture_pulses[t2+1].time-admixture_pulses[t2].time)/n) ;
                    }
                    
                    /// now deal with selecting lineage of the correct ancestry type
                    for ( int a = t ; a > s2 ; a -- ) {
                        rate *= ( 1 - admixture_pulses[a].proportion ) ;
                    }
                    
                    /// now select the correct ancestry type
                    rate *= admixture_pulses[s2].proportion ;
                    
                    /// and augment by this rate for this epoch
                    transition_rates(admixture_pulses[s2].entry_order,admixture_pulses[s1].entry_order) += rate ;
                }
            }
        }
    }
    
    return transition_rates.t() ;
}

#endif
