/*

 copyright: Russ Corbett-Detig
            rucorbet@ucsc.edu
 
 This is software distributed under the gnu public license version 3.
 
 */

/// headers
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std ;

/// linear algebra library is armadillo
#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma ;

/// our header files in /src directory
#include "print_usage.h"
#include "factorial.h"
#include "nchoosek.h" 
#include "subsample.h" 
#include "multichoose.h"
#include "multipermute.h"
#include "normalize.h"
#include "ancestry_pulse.h"
#include "ploidy_path.h" 
#include "markov_chain.h"
#include "read_samples.h" 
#include "pulses_to_ancestry.h" 
#include "compute_forward.h"
#include "compute_backward.h"
#include "forward_backward.h"
#include "viterbi.h" 
#include "transition_information.h"
#include "exponentiate_matrix.h"
#include "cmd_line.h"
#include "create_transition_rates.h"
#include "read_cmd_line.h"
#include "evaluate_vertex.h"
#include "check_vertex.h"
#include "sort_vertices.h"
#include "create_pulses.h" 
#include "create_states.h"
#include "input_line.h"
#include "distribute_alleles.h" 
#include "binomial.h"
#include "read_emissions.h"
#include "genotype_emissions.h"
#include "read_input.h"
#include "nelder_mead.h"
#include "golden_search.h"
#include "bootstrap.h" 
#include <iomanip>

#include "selection_optimization_context.h"
#include "simulate_selective_admixture.h"
#include "new_read_input.h"
#include "nelder_mead_general.h"
#include "optimize_selection.h"

#include "peak_examining.h"
#include "read_gss.h"


int main ( int argc, char *argv[] ) {
    
    /// time tracking
    clock_t t = clock() ;
    clock_t total = clock() ;
    
    /// seed prng
    srand (t) ;
    
	// read cmd line 
	cmd_line options ;
    cerr << "reading command line" ; t = clock();
	options.read_cmd_line( argc, argv ) ;

    /// chain objects for each sample
    vector<markov_chain> markov_chain_information ;
    
    /// get sample ids and ploidy from input file
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading sample ids and ploidy" ; t = clock();
    read_samples( markov_chain_information, options.sample_file, options.viterbi ) ;

    /// create states matrix
    cerr << "\t\t\t" << (double) (clock() - t) << " ms\n" << "creating states matrix" ; t = clock();
    /// store all possible state space arranged by ploidy and then vector of state counts
    map<int,vector<vector<int> > > state_list ;
    /// now create initial state list
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_initial_states( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, options.ancestry_pulses, state_list ) ;
        }
    }
    
	/// read in panels and update matrices
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading data and creating emissions matrices\t" ; t = clock() ;
    /// store recombination rates and positions
    vector<int> position ;
    vector<double> recombination_rate ;
    vector<string> chromosomes ;


    
    if(options.selection_mode){
        selection_read_file( options, markov_chain_information, state_list, position, recombination_rate, chromosomes ) ;
    }else{
        read_file( options, markov_chain_information, state_list, position, recombination_rate, chromosomes ) ;
    }
    
    //removing first element, beware
    //recombination_rate.erase(recombination_rate.begin());
    recombination_rate[0] = 0;

    //cout << "recombination_rate size " << recombination_rate.size() << "\n";
    //cout << "position size " << position.size() << "\n";
    //cout << "last position" << position[36172];


    //checking sum of recombination rates
    double recomb_sum = 0;
    for(int i = 0; i < recombination_rate.size(); i++){
        recomb_sum += recombination_rate[i];
        if( recomb_sum >= 0.1999 && recomb_sum <= 0.2001){
            //cout << "LOOKERE " << recomb_sum << " " << i << "  " << position[i] << "\n";
        }
    }
    //cout << "Sum of recombination_rate " << recomb_sum << "\n";
    
    
        
    /// create basic transition information
    cerr << (double) (clock() - t) << " ms" << endl << "computing transition routes\t\t\t" ; t = clock() ;
    /// 3d map to look up by ploidy, start state, end state, and then relevant transition information
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_transition_information( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, transition_matrix_information, state_list[markov_chain_information.at(m).sample_ploidy_path[p].ploidy] ) ;
        }
    }
    
    /// create admixture model(s)
    cerr << (double) (clock() - t) << " ms" << endl << "creating initial admixture model(s)\t\t" ; t = clock();
    vector<vector<pulse> > vertices ;
    int nparams = create_pulses( vertices, options ) ;
    
    /// set number of restarts if unspecified default is factorial * 2
    cerr << (double) (clock() - t) << " ms" << endl << "estimating " << nparams << " parameters\n" ;
    if ( options.n_restarts < 0 ) {
        options.n_restarts = factorial[nparams] * 2 ;
    }
    
    /// vector of models to be evaluated and optimized
    vector<pulse> optimum ;
    /*
    /// if there are params to estimate, do amoeba search
    if ( nparams > 1 ) {
        cerr << "starting nelder-mead search\t\t" << endl ;
        optimum = nelder_mead_search( vertices, options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
        cerr << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
    }
    
    /// or do golden section line search for single parameter optimization
    else if ( nparams == 1 ) {
        cerr << "starting golden section search\t\t" << endl ;
        optimum = golden_search( options, markov_chain_information, transition_matrix_information, recombination_rate, position, state_list ) ;
        cerr << "\n\t\t\t\t\tSEARCH TIME: " << (double) (clock() - t) << " ms" << endl << endl << "optimal model found:\n\n" ;
    }
    
    /// otherwise just evaluate the supplied model
    else {
        optimum = options.ancestry_pulses ;
        cerr << endl << endl << "evaluating supplied model:\n\n" ;
    }
    */
    
    /// print model
    cerr << "\ttype\ttime\tproportion\n" ;
    //cout << "optimum: \n" ;
    //cout << "\ttype\ttime\tproportion\n" ;
    /*
    vector<double> a = options.ancestry_proportion ;
    for ( int p = 0 ; p < optimum.size() ; p ++ ) {
        optimum[p].proportion = a[optimum[p].type] * optimum[p].fraction_of_remainder ;
        a[optimum[p].type] -= optimum[p].proportion ;
        cerr << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
        cout << "\t" << optimum[p].type << "\t" << optimum[p].time << "\t" << optimum[p].proportion << endl ;
    }
    
    /// bootstrap models as necessary
    if ( options.n_bootstraps > 0 ) {
        cerr << "computing " << options.n_bootstraps << " bootstrap models" << endl ;
        
        vector<vector<pulse> > bootstrap = bootstraps( vertices, markov_chain_information, transition_matrix_information, recombination_rate, position, options, state_list, chromosomes ) ;
        
        //// print out bootstrapped admixture models 
        for ( int b = 0 ; b < options.n_bootstraps ; b ++ ) {
         
            cout << "bootstrap: " << b << endl ;
            cout << "\ttype\ttime\tproportion\n" ;
            
            cerr << endl << "bootstrap: " << b << endl ;
            cerr << "\ttype\ttime\tproportion\n" ;
            
            vector<double> a = options.ancestry_proportion ;
            for ( int p = 0 ; p < optimum.size() ; p ++ ) {
                bootstrap[b][p].proportion = a[optimum[p].type] * bootstrap[b][p].fraction_of_remainder ;
                a[bootstrap[b][p].type] -= bootstrap[b][p].proportion ;
                cerr << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
                cout << "\t" << bootstrap[b][p].type << "\t" << bootstrap[b][p].time << "\t" << bootstrap[b][p].proportion << endl ;
            }
        }
    }
    */

    /// create transition rates for the optimal or supplied set of pulses
    cerr << endl << "creating per morgan transition rates\t\t" ; t = clock();
    //call function 
    //mat transition_rates = create_transition_rates( optimum, options.ne, options.ancestry_proportion ) ;
    
    /// create transition information
    cerr << (double) (clock() - t) << " ms" << endl << "creating transition matrices\t\t\t" ; t = clock();
    map<int,vector<mat> > transition_matrix ;


    /*
    selection_opt testest(recombination_rate, options, markov_chain_information, transition_matrix_information, position);
    testest.set_context();
    vector<double> params(3);
    params[0] = 0.2;
    params[1] = 1;
    params[2] = 0.96;
    double out = to_be_optimized(params);
    cerr << "TESTEST output: " << out << "\n";
    //testest.enact_optimization();
    */
    

    vector<mat> transition_matrices;

    if (options.selection_mode && !options.optimize_selection){
        
        vector<double> selection_recomb_rates(options.selected_sites.size() + 1);
        
        double last = 0;

        for (int i = 0; i < options.selected_sites.size(); i++){
            selection_recomb_rates[i] = options.selected_sites[i] - last;
            last = options.selected_sites[i];
        }

        selection_recomb_rates[options.selected_sites.size()] = 1 - last;


        //trying out recombination_rate instead of 
        transition_matrices = calculate_transition_rates(
            recombination_rate,
            selection_recomb_rates,
            options.fitnesses,
            options.m,
            options.generations
        );

        //cout << "Transition rates calculated " << transition_matrices.size() << "\n";
        
    }

    if(options.optimize_selection){
        selection_opt selection_optimizer(recombination_rate, options, markov_chain_information, transition_matrix_information, position);
        selection_optimizer.enact_optimization();
        return 0;
    }
    if(options.use_ahmm_s_peaks){
        read_gss(options, recombination_rate, position);
        selection_opt selection_optimizer(recombination_rate, options, markov_chain_information, transition_matrix_information, position);
	selection_optimizer.examine_peaks();
    }

    //cerr << transition_matrices[0] << " trans mat\n";

    //THIS NEWLINE IS HIGHLY CRITICAL, DONT ASK WHY
    cerr << '\n';



    /* Prints out expected ancestry
    double window_size = 0.01;
    
    double last_window = -window_size;

    double forward_sum = 0;

    double forward_ancestry_proportion = options.m;

    vector<double> locations;
    vector<double> forward_props;

    

    for(
        int i = 0;
        i < transition_matrices.size();
        i++
    ){
        forward_sum += recombination_rate[i];

        if(forward_sum >= last_window+window_size){
            locations.push_back(forward_sum);
            forward_props.push_back(forward_ancestry_proportion);
            last_window += window_size;
        }

        double from_ancestry_0 = 0;
        double from_ancestry_1 = 0;


        from_ancestry_0  =      forward_ancestry_proportion  * transition_matrices[i](0,0);
        from_ancestry_0 += (1 - forward_ancestry_proportion) * transition_matrices[i](1,0);

        from_ancestry_1  =      forward_ancestry_proportion  * transition_matrices[i](0,1);
        from_ancestry_1 += (1 - forward_ancestry_proportion) * transition_matrices[i](1,1);

        forward_ancestry_proportion = from_ancestry_0/(from_ancestry_0 + from_ancestry_1);
    }

    for(int i = 0; i < locations.size(); i ++){
        cout << locations[i] << "\t" << forward_props[i] << "\n";
    }
    */


    //cout << "transition_matrices size" << transition_matrices.size() << "\n";
    //cout << transition_matrices.back() << "\n";

    //transition_matrix.insert(pair<int,vector<mat>>(2, transition_matrices));

    
    alt_create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[0].ploidy_switch[0]], recombination_rate, position, markov_chain_information[0].ploidy_switch[0], transition_matrices ) ;

    //cout << "transition_matrix[2] size: " << transition_matrix[2].size() << "\n";
    //cout << transition_matrix[2].back() << "\n";
    /*
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information.at(m).number_chromosomes], recombination_rate, position, markov_chain_information.at(m).number_chromosomes, transition_rates ) ;
        for ( int p = 0 ; p < markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            create_transition_matrix( transition_matrix, transition_matrix_information[markov_chain_information[m].ploidy_switch[p]], recombination_rate, position, markov_chain_information[m].ploidy_switch[p], transition_rates ) ;
        }
    }./
    */
    
    //// create interploidy transition matrix
    vector<mat> interploidy_transitions; // = create_interploidy_transitions ( state_list, optimum, options.ancestry_proportion ) ;
    
    /// output viterbi path for optimized model
    /*if ( options.viterbi == true ) {
        cerr << (double) (clock() - t) << " ms" << endl << "viterbi posterior decoding and printing\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].viterbi( position, recombination_rate, state_list, chromosomes, transition_matrix, interploidy_transitions, options.output_pulses, optimum ) ;
        }
    }*/

    /// output forward-backward full probability distribution by default
    
        cerr << (double) (clock() - t) << " ms" << endl << "computing forward probabilities\t" ; t = clock() ;
        double lnl = 0 ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            lnl += markov_chain_information[m].compute_forward_probabilities( transition_matrix, interploidy_transitions ) ;
        }
        cout << setprecision(15) << lnl << "\n";
        cerr << "lnl: " << setprecision(15) << lnl << "\t\t" << (double) (clock() - t) << " ms" << endl << "computing backward probabilities\t\t\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].compute_backward_probabilities( transition_matrix, interploidy_transitions ) ;
        }
        
        //Nico has commented out the posterior printing
        /*
        cerr << (double) (clock() - t) << " ms" << endl << "forward-backward posterior decoding and printing\t\t\t" ; t = clock() ;
        for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
            markov_chain_information[m].combine_prob( position, state_list, chromosomes, options.output_pulses, optimum ) ;
        }
        */
    
        
    cerr << (double) (clock() - t) << " ms" << endl ;
    cerr << "total run time:\t\t\t" << (double) (clock() - total) << " ms" << endl ;

	return 0 ; 
}
	

