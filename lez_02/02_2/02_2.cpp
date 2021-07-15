//
//  02_2.cpp
//  Simulation of random walks and estimation of
//  average square displacement
//
//  Created by Andrea Caprotti on 27/03/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include "../../libraries/block_stat.h"
#include "../../libraries/random.h"
#include "../../libraries/random_walk.h"

int main(int argc,char *argv[]){

    Random rnd;
    
    int M = 10000;                  // random walks started
    int N = 100;                    // number of blocks
    int n_steps = 100;              // number of steps of walk
    int i,k;                        // cycle and walk index
    
    std::ofstream discrete_walk("discrete_walk.txt");
    std::ofstream cont_walk ("continuous_walk.txt");
    std::ofstream walk_example_c("walk_example_c.txt");
    std::ofstream walk_example_d("walk_example_d.txt");
    
    std::ofstream d_walk_nb("discrete_walk_nb.txt");
    std::ofstream c_walk_nb("continuous_walk_nb.txt");
    
    std::ofstream s_w ("bunch_of_walks.txt");
    
    Random_Walk walk_d[M];  // initialises vector of walks
    Random_Walk walk_c[M];
    
    for (i=0; i<M; ++i){
        walk_d[i].start(3);
        walk_c[i].start(3);
    }
    
    std::vector<double> block_c; // always set up when staring
    block_c.clear();
    std::vector<double> block_d;
    block_d.clear();
    
    double a = 1;                // step increment lenght
    Random ran;                  // determines random variables
                                 // (not very elegant, otherwise
                                 // every walk risks to be same)
    
    // first position always 0
    discrete_walk << 0 << std::setw(12) << 0 << std::endl;
    cont_walk << 0 << std::setw(12) << 0 << std::endl;
    d_walk_nb << 0 << std::setw(12) << 0 << std::endl;
    c_walk_nb << 0 << std::setw(12) << 0 << std::endl;
    
    for (i=0; i < n_steps; i++){       // repeat for each
       // saves progress of random walk in vector
        walk_c[5].print_position(walk_example_c);
        walk_d[39].print_position(walk_example_d);
        
        block_d.clear();
        block_c.clear();

        for (k = 0; k < M; ++k){       // block division
            // discrete step for each walk
            walk_d[k].discrete_step_increment(a,ran.Rannyu(), ran.Rannyu(),100);
            block_d.push_back(sqrt(walk_d[k].vect_length()));
            
            // continuous step for each walk
            walk_c[k].continuous_step_increment(a, ran.Rannyu(), ran.Rannyu());
            block_c.push_back(sqrt(walk_c[k].vect_length()));
            
            if (M%(k+1)==1000)     // random representative walks
                s_w << block_d[k] << std::setw(12);
        }
        s_w << std::endl;
        
        full_block_average (block_d, discrete_walk, N, false);
        full_block_average (block_c, cont_walk, N, false);
        
        // attempt to estimate error without block average
        d_walk_nb << average_vec(block_d) << std::setw(12) << sqrt(std_dev_vec(block_d)) << std::endl;
        c_walk_nb << average_vec(block_c) << std::setw(12) << sqrt(std_dev_vec(block_c)) << std::endl;
        
    }
    
    discrete_walk.close();
    cont_walk.close();

    walk_example_c.close();
    walk_example_d.close();
    
    d_walk_nb.close();
    c_walk_nb.close();
    
    s_w.close();
    
    return 0;
}
