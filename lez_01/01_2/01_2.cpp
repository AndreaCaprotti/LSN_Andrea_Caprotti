//
//  01_2.cpp
//  Creates histograms to check Central Limit Theorem
//  using different probability distributions
//
//  Created by Andrea Caprotti on 23/03/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "../../libraries/random.h"

using namespace std;

int main (int argc, char *argv[]){

    Random rnd;

    vector<std::string>file_name={"standard_", "exponential_", "lorentzian_"};
    std::string suffix;
    int index_avg[5]={1,2,5,10,100};
    
    int M = 10000;          // number of indipendent elements
    int avg_N;
    
    double lambda = 1;      //exponential parameter
    double Gamma = 1;       //lorentzian paramenters
    double mean = 0;
    double sum;
    
    int name, i, j, n;      // indexes defined here to releive
                            // computational stress on cycles
    for (name = 0; name < file_name.size(); name++){
        for (i = 0; i<5; i++){
            avg_N = index_avg[i];
            suffix = to_string(avg_N);
            ofstream histogram_file (file_name[name]+suffix+".txt");
            //generates new file for every combination
            for (j=0; j<M; j++){
                sum = 0;
                for (n=0; n<avg_N; n++){
                    if (name == 0)        //std distribiution
                        sum += rnd.Rannyu();
            
                    else if ( name == 1)   //exponential
                        sum += rnd.Exponential(lambda);
                   
                    else if (name == 2)     //lorentzian
                        sum += rnd.Lorentzian(Gamma, mean);
                }
                sum = sum / avg_N;
                //average made differently
                histogram_file << sum << endl;
            }
            //histogram will actually be built using Python
            histogram_file.close();
        }
    }
    rnd.SaveSeed();
    return 0;
}
