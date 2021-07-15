//
//  02_1.cpp
//  Estimation of integral between 0 and 1 of pi/2 cos(x pi/2)
//  Created by Andrea Caprotti on 25/03/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "../../libraries/block_stat.h"
#include "../../libraries/random.h"

int main(int argc,char *argv[]){
    
    Random rnd;
    
    int M = 100000;      // total number of steps
    int N = 100;         // number of blocks
    int L = int(M/N);
    double pih=M_PI/2.;  // just for laziness
    int i,j;             // cycle indexes
    
    // Estimation using uniform distribution
    std::vector<double> block_avg;      // Block averages

    // use 2(1-x) as probability distribution to sample
    std::vector<double> block_sample;
    double distr_sample;                //placeholder variable

    std::ofstream integral_data ("std_integral.txt");
    std::ofstream integral_samp ("sample_integral.txt");
    
    for (i = 0; i < N; i++) {
        for (j=0; j<L; j++){
            block_avg.push_back(pih*cos(pih*rnd.Rannyu()));
            
            distr_sample = 1 - sqrt(1 - rnd.Rannyu());
            block_sample.push_back(pih*cos(pih*distr_sample)/(2-2*distr_sample));
        }
        block_average(block_avg, integral_data);
        block_average(block_sample, integral_samp);
    }
    integral_data.close();
    integral_samp.close();

    return 0;
}
