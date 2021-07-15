//
//  10_1.cpp
//  Travelling Salesman Problem solved using
//  simulated annealing. Uses same structures of
//  Chromosomes used for the genetic algorithm
//
//  Created by Andrea Caprotti on 30/06/21.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include "10_1.hpp"

int main (int argc, char* argv[]){
    // first portion is exaclty what also happens
    // in the genetic algorithm
    
    Random rnd;
    srand(0);
    City foo_bar_ville;
    Chromosome route(n_city);
    
    double foo_bar;
    int foo, bar;
    int i, j;                       // cycle index
    bool goodness;                  // check for requirements
    
    //                                                  //
    // random cities generated in a square of side l=1  //
    // assings first chromosome                         //
    //                                                  //
    std::ofstream cities ("city_coordinates.txt");
    
    for (int i=0; i < n_city; ++i){    // city initialisation
        // generates cities inside unit square
        foo_bar_ville.x = (rand()/(double)RAND_MAX);
        foo_bar_ville.y = (rand()/(double)RAND_MAX);
        
        foo_bar_ville.index = i;
        country.push_back(foo_bar_ville);
        
        cities << std::setw(3) << foo_bar_ville.index << std::setw(12) << foo_bar_ville.x << std::setw(12) << foo_bar_ville.y << std::endl;
        
        // initialises length matrix (to then set weight)
        for(j = 0; j < i; ++j){
            len_matrix[i][j] = pow(country[i].x - country[j].x,2) + pow(country[i].y - country[j].y,2);
            // distance between i-th and j-th city
            len_matrix[j][i] = len_matrix[i][j]; //symmetric matrix…
            // for j=i distance is alreay set as 0
        }
        
        route.gene_sequence.push_back(i);
        // first chromosome is just ordered sequence of random cities
    }
    cities.close();
    
    route.weight = SetWeight(route);    // first weight set
    
    // repetition of code sets temperature from outside
    if (argc > 1){
        beta = atof(argv[1]);
        n_steps = pow(beta,1.5)*n_steps;
    }
    
    // evolves route with random mutations
    // acceptance of step (when allowed) w/ Boltzmann weight
    // (route length acts as energy)
    
    double weight_comp, new_weight;    // placeholder variables
    int choice;
    std::vector<int> placeholder;

    for (i=0; i < n_steps; ++i){
    
        placeholder = route.gene_sequence; //failsafe
        weight_comp = route.weight;
        
        // choice of random mutation
        foo = (int)((n_city - 1)*rnd.Rannyu()) + 1;
        bar = (int)((n_city - 1)*rnd.Rannyu()) + 1;
        // indexes of mutation (avoid first element)
        // second index always greater than first
        
        choice = (int)rnd.Rannyu(0,3);  // choice of mutation
        
        switch (choice) {
            case 0:
                route.SinglePermutation (foo, bar);
                break;
                
            case 1:
                bar = (int)(n_city/2 * rnd.Rannyu());
                //max genes that can be swapped
                route.BlockPermutation (bar, foo);
                break;
                
            default:
                if (foo < bar)
                        route.PartialReverse (foo, bar);
                else
                        route.PartialReverse (bar, foo);
                break;
        }
        
        goodness = route.GoodChromo(); // goodness check
        new_weight = SetWeight(route);
        
        if (goodness){            // conditions to accept move
            if ( new_weight < weight_comp){
                accepted++;
                route.weight = new_weight;
                continue;
            }
            
            foo_bar = exp(-beta*(new_weight-weight_comp));
            if (rnd.Rannyu() < foo_bar){
                accepted++;
                route.weight = new_weight;
                continue;
            }
        }
        
        // if move isn't accepted (or it doesn't fit
        // requirements) move is reverted
        
        route.FillChromo(placeholder);
    }
    
    std::ofstream Out("SA_STD.txt", std::ios::app);
    Out << std::setw(12) << beta << std::setw(12) << route.weight << std::endl;
    Out.close();
    
    // after each repetition updates best sequence
    std::ofstream BestRoute("final_route.txt");
    route.PrintSequence(BestRoute);
    BestRoute.close();
    
    std::cout << accepted/(double)n_steps << " acceptance rate" << std::endl;
    
    return 0;
}
