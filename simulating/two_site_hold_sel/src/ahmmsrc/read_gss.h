#ifndef READ_GSS
#define READ_GSS
#include <string>
#include <vector>

void read_gss(cmd_line &options, vector<double> &recomb_rates , vector<int> &positions){
    
    // stream in file
    ifstream in ( options.gss_output.c_str() );
    
    options.gss_out_pos.resize(0);
    options.gss_out_sel.resize(0);
    options.gss_out_llr.resize(0);

    double chrom_pos = 0;
    int i = 0;

    while ( !in.eof() ) {
        double pos, sel;
	string llr_str;
	double llr;
        in >> pos >> sel >> llr_str;
	
	try{
	    llr = stod(llr_str);
	}catch(exception &e){
	    llr = 0;
	}

	
        for(;positions[i] < pos; i++){
            chrom_pos += recomb_rates[i];
        }
        

        options.gss_out_pos.push_back(chrom_pos);
        options.gss_out_sel.push_back(sel);
        options.gss_out_llr.push_back(llr);
	
    }
}

#endif
