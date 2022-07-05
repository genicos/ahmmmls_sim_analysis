#ifndef __CMD_LINE_H
#define __CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// terms to bound the optimization
    double t_max ;
    double t_min ;
    
    /// to bound proportion search
    double p_max ;
    double p_min ;
    
    /// to create intial simplex points
    double t_length ;
    double p_length ;
    
    /// number of restarts
    int n_restarts ;
    
    /// proportion ancestry for 0-n must sum to 1
    /// these are therefore the final ancestry proportion, not necessary the proportion that fluxed if additional pulses occured closer to the present
    vector<double> ancestry_proportion ;
    
    /// store relevant ancestry information
    vector<pulse> ancestry_pulses ;
    
    /// diploid effective population size ( i.e. 2n )
    double ne ;
    
    /// tolerance for parameter search
    double tolerance ;
    
    /// minimum recombinational distance between markers
    double minimum_distance ;
    
    /// error rates for reads (if read based) or genotypes (if genotype based)
    double error_rate ;
    
    /// bool sample is expressed as genotypes, not read counts
    bool genotype ;
    
    /// ancestral genotype frequencies are fixed
    bool ancestral_fixed ;
    
    /// viterbi output
    /// caution: not recommended for more samples of ploidy > 1
    bool viterbi ;
    
    /// number of digits of precision to include
    int precision ;
    
    /// output actual pulses rather than ancestry states
    bool output_pulses ;
    
    /// error rates specifed
    bool error_rates ; 
    
    /// input file name
    string input_file ;
    
    /// sample file
    string sample_file ;
    
    /// bootstrap
    int n_bootstraps ;
    int block_size ; 


    

    bool admixture_parameters_set = false;

    // given ancestry proportion
    double m;

    // given generations
    int generations;

    

    // Place sites one at a time, with no info on selected sites
    bool uninformed_inference = false;
    


    // examine sites from site file
    bool use_site_file = false;

    //sites file
    string site_file;
    vector<int>    site_file_positions;
    vector<double> site_file_morgan_positions;
    vector<string> site_file_options;

    
    // calculate likelihoods for models
    bool use_model_file = false;

    string model_file;
    vector<vector<double>> models;

    int cores = 1;


    bool verbose_stdout = false;
    bool verbose_stderr = false;


    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;

} ;

#endif

