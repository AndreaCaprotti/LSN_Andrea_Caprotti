#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../libraries/block_stat.h"
#include "../../libraries/random.h"

int main (int argc, char *argv[]){

    Random rnd;
   
    int M = 1000000;         // Total number of random throws
    int N = 100;            // Number of blocks
    int L = int(M/N);       // Number of throws per block
    
    std::vector<double> block_avg;
    std::vector<double> block_var_avg;
    std::ofstream avg_data ("progressive_averages_1.txt");
    std::ofstream var_data ("progressive_variance_1.txt");
    
    // defines histogram vector for chi square estimation
    int bins = 100;
    std::vector<int> histogram(bins,0);
    std::ofstream chi_sq("chi_sq.txt");
    double chi_square;
    
    double foo_bar;         // placeholder variable
    double sum;             // placeholder for average in a block
    double var_sum;         // same for estimation of variance
    int hist_index,i,j,l;   // cycle indexes
    
    for (i=0; i < N; i++){
        sum = 0;
        var_sum = 0;
        fill(histogram.begin(), histogram.end(), 0);
        chi_square = 0;
        
        // generate L = 10000 (pseudo-)random variables
        // for each block
        for(j=0;j<L;j++){
            foo_bar = rnd.Rannyu();
            sum += foo_bar;
            var_sum += pow(rnd.Rannyu()-0.5,2);
            hist_index = floor(foo_bar*bins);
            histogram[hist_index]++;
        }
        // progressive average on blocks
        block_avg.push_back(sum/L);
        block_average(block_avg,avg_data);
        
        block_var_avg.push_back(var_sum/L);
        block_average(block_var_avg,var_data);
        
        //chi square evaluation on each block
        for (l=0; l<bins; l++)
            chi_square += pow(histogram[l]-L/bins, 2)/(L/bins);
        chi_sq << chi_square << std::endl;
    }
    
    rnd.SaveSeed();
    chi_sq.close();
    avg_data.close();
    var_data.close();
    

    return 0;
}
