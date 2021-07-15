//
//  05_1.cpp
//  Metropolis algorithm for estimation of atomic radius
//  of hydrogen atom ground state and 2p state
//
//  Created by Andrea Caprotti on 22/04/21.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../../libraries/random_walk.h"
#include "../../libraries/block_stat.h"

// Bohr radius fixed as 1 so it won't appear in calculations
// Probability ditribution specific to problem, so defined here
double gs_probability(double radius){
    return exp(-2*radius)/M_PI;
}

double e2p_probability (double radius, double z){
    return z*z*exp(-radius)/(32*M_PI);
}

int main(int argc, char *argv[]){
    
    int i;                      // cycle index
    int M = 100000;             // Number of total steps
    int N = 100;                // Number of blocks
    
    std::vector<double> start_position ={0.,0.,10.}; // nomen omen
    // doesn't start froim origin since wavefunction diverges there
    double save_z;              // used to obtain coord z for 2p

    // Metropolis_Walk class, which inherits moves from random walk
    // class but adds ability to accept or refuse a step
    Metropolis_Walk gs_variable(start_position);
    Metropolis_Walk e2p_variable(start_position);
    
    std::vector<double> gs_steps;  // vectors to store radius (for
    std::vector<double> e2p_steps; // block average afterwards…)
    double p_gs, p_e2p, radius;    // storage var. for probability
    // and for radius
    
    //double fs_gs  = 2.755;       // gs scale factor
    //double fs_e2p = 6.85;        // e2p scale factor
                                 // to have acceptance at 50%
    double fs_gs = 0.755;       // optimal values for gaussian
    double fs_e2p = 1.755;      // steps
    
    // assigns probability for first configuration
    radius = gs_variable.vect_length(); // only one vector
    // element (index=0)
    p_gs = gs_probability(radius);
    gs_variable.assign_prob(p_gs);
    
    // same but for excited states
    radius = sqrt(gs_variable.vect_length());
    save_z = e2p_variable.current_position(2); // gets coord in
    // 3rd direction (z)
    p_e2p = e2p_probability(radius, save_z);
    e2p_variable.assign_prob(p_e2p);
    
    std::ofstream gs_position("gs_position.txt");
    std::ofstream e2p_position("e2p_position.txt");
    
    for(i=0; i<M; i++){
        
        //gs_variable.new_move(fs_gs); // uniform prob step for each direction
        //e2p_variable.new_move(fs_e2p);
        
        gs_variable.new_gauss_move(fs_gs); // gaussian prob step for each direction
        e2p_variable.new_gauss_move(fs_e2p);
        
        // maybe not the most condensed writing but clean
        radius = sqrt(gs_variable.vect_length());
        p_gs = gs_probability(radius);
        gs_variable.accept_move(p_gs);
        gs_steps.push_back(sqrt(gs_variable.vect_length()));
        gs_variable.print_position(gs_position);
        
        
        radius = sqrt(e2p_variable.vect_length());
        save_z = e2p_variable.current_position(2);
        p_e2p = e2p_probability(radius, save_z);
        e2p_variable.accept_move(p_e2p);
        e2p_steps.push_back(sqrt(e2p_variable.vect_length()));
        e2p_variable.print_position(e2p_position);
    }
    
    std::cout << gs_variable.no_acc/(double)M << " " << e2p_variable.no_acc/(double)M << std::endl;
    // control of acceptance rate
    
    gs_position.close();
    e2p_position.close();
    
    std::ofstream gs_avg ("gs_radius_average.txt");
    full_block_average(gs_steps, gs_avg, N, true);
    gs_avg.close();
    
    std::ofstream e2p_avg ("e2p_radius_average.txt");
    full_block_average(e2p_steps, e2p_avg, N, true);
    e2p_avg.close();
    
    return 0;
}
