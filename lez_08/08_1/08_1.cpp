//
//  08_1.cpp
//  Programme estimates the energy of a generic wavefunction in a
//  double well potential, in order to determine the best
//  parameters for the ground state.
//  Average value of energy estimated as a Monte Carlo integral,
//  using as PDF the square of the wavefunction, sampled using
//  a Metropolis algorithm.
//
//  Created by Andrea Caprotti on 20/05/21.
//

#define PRINT               // prints histogram and progressive energy average
//#define GET_PARAMS          // gets parameter mu and sigma from command line
#define EQUIL               // inlcudes "equilibration phase" for acceptance

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../../libraries/random_walk.h"
#include "../../libraries/block_stat.h"

// density probability for wavefunction
double square_wf (double x, double mu, double s_2){
    double pos_exp, neg_exp;
    pos_exp = pow((x + mu),2)/(2 * s_2);
    neg_exp = pow((x - mu),2)/(2 * s_2);
    
    return exp(-2 * pos_exp) + exp(-2 * neg_exp) + 2 * exp( -pos_exp - neg_exp);
}

// integrand to be estimated (expected value of H\psi / \psi)
double kin_energy (double x, double m, double s_2){
    return 1./(2*s_2) - ((x + m)*(x + m)*exp(- (x + m)*(x + m)/(2*s_2)) + (x - m)*(x - m)*exp(- (x - m)*(x - m)/(2*s_2)))/(2*s_2*s_2*(exp(- (x + m)*(x + m)/(2*s_2)) + exp(- (x - m)*(x - m)/(2*s_2))));
}

double pot_energy (double x){
    return pow(x,4) - 2.5 * pow(x,2);
}

double energy (double x, double mu, double s_2){
    return pot_energy(x) + kin_energy(x,mu,s_2);
}

int main(int argc, char *argv[]){
    
    int M = 100000;             // total number of steps for walk
    int N = 100;                // number of blocks
    int i;                      // index variable
    double mu, s;               // wf parameters
    bool fba;                   // for full block average
    std::string file_name = "data_8_1/min_energy_search.txt";
    double fs=0.5;
    
    fba = false;
    
#ifdef PRINT                // saves progressive energy avg
    fba = true;
    file_name = "data_8_1/energy_avg.txt";
    // wf histogram
    std::ofstream hist ("data_8_1/wavefunction.txt");
#endif
    
#ifndef GET_PARAMS
    mu = 0.811299;              // optimal parameters for wavefunction
    s = 0.62025;             // used for progressive energy avg
#endif
    
#ifdef GET_PARAMS
    // uses launcher.sh to estimate many params
    mu = atof(argv[1]) + (atof(argv[3]) - atof(argv[1]))* atof(argv[2])/atof(argv[7]);
    s = atof(argv[4]) + (atof(argv[6]) - atof(argv[4]))*  atof(argv[5])/atof(argv[7]);
#endif
    
    // wavefunction sampling
    std::vector<double> start_position = {0};  // single dof
    // clever position to start due to symmetry
    
    Metropolis_Walk wf(start_position);
    double save_prob;                   // storage variable
    
    // initiates random walk
    save_prob = square_wf(wf.current_position(0), mu, s*s);
    wf.assign_prob(save_prob);
    
    // integral sampling
    std::vector<double> en_average;
    en_average.clear();                 // never too sure…
    std::ofstream out_energy;
    out_energy.open(file_name, std::ios::app);    
    
#ifdef EQUIL
    // "equilibration": a series of attempts to determine the
    // optimal step lenght factor fs for given parameters
    // (not exact, but at least about 50% acceptance)
    double acc_rate;
    int n = 1000;       // equilibration time
    int k=0;
    do{
        Metropolis_Walk trial(start_position);
        trial.assign_prob(save_prob); // same initial condition
        for (i=0; i<n; ++i){
            fs = 0.5 + 0.25*k;
            trial.new_move(fs);
            save_prob = square_wf(trial.current_position(0), mu, s*s);
            trial.accept_move(save_prob);
        }
        k++;
        acc_rate = trial.no_acc / (double)n;
    } while(acc_rate > 0.55);
#endif
    
    //                   //
    // actual simulation //
    //                   //
    
    for (i=0; i<M; ++i){
        // Metropolis step to evaluate wf evolution
        wf.new_move(fs);
        save_prob = square_wf(wf.current_position(0), mu, s*s);
        wf.accept_move(save_prob);
        
        // Energy evaluation to be saved for average
        en_average.push_back(energy(wf.current_position(0), mu, s*s));
        
#ifdef PRINT
        wf.print_position(hist);
#endif
    }
    
#ifdef GET_PARAMS
    out_energy << mu << std::setw(15) << s << std::setw(15);
#endif
    
    full_block_average(en_average, out_energy, N, fba);
    
   // std::cout << "Acceptance rate:" << wf.no_acc/(double)M << std::endl;

    out_energy.close();
#ifdef PRINT
    hist.close();
#endif
    return 0;
}
