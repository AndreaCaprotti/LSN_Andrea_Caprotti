//
//  03_1.cpp
//  Simulation of GBM in order to establish optimal put and call
//  prices after a given time T
//
//  Created by Andrea Caprotti on 08/04/21.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../libraries/block_stat.h"
#include "../../libraries/random.h"

double price_evolution (double old_price, double r, double s, double t, double gauss_m){
    return old_price*exp((r - (s*s)/2)*t + gauss_m*s*sqrt(t));
}
                         
double call (double price, double strike, double r, double t){
    if (price - strike > 0){
        return exp(-r*t)*(price - strike);
    }
    else
        return 0;
}

double put (double price, double strike, double r, double t){
    if (strike - price > 0)
        return exp(-r*t)*(strike - price);
    else
        return 0;
}

int main(int argc,char *argv[]){
    
    Random rnd;
    int M = 10000;                  // random walks started
    int N = 100;                    // number of blocks
    int n_steps = 100;
    int i,k;                        // cycle index

    // different containers to distinguish the various
    // possible combinations
    // different file for each combination for better order
    std::ofstream call_direct ("direct_call.txt");
    std::vector<double> call_dir_vec;
    
    std::ofstream call_discrete ("discrete_call.txt");
    std::vector<double> call_disc_vec;

    std::ofstream put_direct ("direct_put.txt");
    std::vector<double> put_dir_vec;

    std::ofstream put_discrete ("discrete_put.txt");
    std::vector<double> put_disc_vec;

    std::ofstream distr_dir ("direct_distribution.txt");
    std::ofstream distr_disc ("discrete_distribution.txt");
    
    double start_price = 100 ;
    double strike_price = 100;
    double T = 1;
    double rate = 0.1;
    double volatility = 0.25;
    
    double price_disc;         // palceholder variables for saving
    double price_dir;          // temporarely price value
    double time_incr = 1/double(n_steps);
    
    // MonteCarlo simulation to fill vectors w/ estimates of optimal price for call or put
    for (i=0; i<M; i++){
        price_disc = start_price;
        // price_disc to be updated for each step
        for (k = 0; k<n_steps; k++)
            price_disc = price_evolution(price_disc,rate,volatility,time_incr, rnd.Gauss(0,1) );
        
        price_dir = price_evolution(start_price,rate,volatility,T, rnd.Gauss(0,1));
        // price_dir and price_disc contain the price at time T
        // computed respectively directly or w/ step evolution
        
        // values saved directly to check total distribution
        // (particularly that they are the same…)
        distr_dir << price_dir << std::endl;
        distr_disc << price_disc << std::endl;
        // now they are simply used to evaluate calls and puts
        call_dir_vec.push_back(call(price_dir, strike_price, rate,T));
        call_disc_vec.push_back(call(price_disc, strike_price, rate, T));
        put_dir_vec.push_back(put(price_dir, strike_price, rate,T));
        put_disc_vec.push_back(put(price_disc, strike_price, rate,T));
    }

    // this same operation is made for each possible combination (call/put & direct/discrete)
    // direct call
    full_block_average(call_dir_vec, call_direct, N, true);

    // call after discrete steps
    full_block_average(call_disc_vec, call_discrete, N,true);

    // direct put
    full_block_average(put_dir_vec, put_direct, N,true);

    // call after discrete steps
    full_block_average(put_disc_vec, put_discrete, N,true);

    call_direct.close();
    call_discrete.close();
    put_direct.close();
    put_discrete.close();
    
    return 0;
}
