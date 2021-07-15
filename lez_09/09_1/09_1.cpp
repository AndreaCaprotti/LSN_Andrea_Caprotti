//
//  09_1.cpp
//  
//
//  Created by Andrea Caprotti on 27/05/21.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>                    // M_PI
#include <string>
#include "09_1.hpp"

//#define CIRCLE
//#define Get_Params                  // to find optimal parameters

int main (int argc, char* argv[]){
    
    Random rnd;
    
    int N = 900;                        // number of possible routes
    City foo_bar_ville;
    Chromosome route(n_city);
        
    double foo_bar;
    int foo, bar;
    int i, j, good_check;               // cycle index
    bool goodness;                      // check for requirements
    std::string filename;
    
#ifdef CIRCLE
    filename = "_circle.txt";
#endif
    
#ifndef CIRCLE
    filename = "_square.txt";
#endif
    
#ifdef Get_Params
    s_swap_p = atof(argv[1])/100;
    m_swap_p = atof(argv[2])/100;
    invert_p = atof(argv[3])/100;
    
    std::ofstream params ("params_compare_s.txt", std::ios::app);
#endif
    
    //                                                   //
    // random cities generated on circle of radius r = 1 //
    // (or in a square of side l = 1)                    //
    // assings first chromosome                          //
    //                                                   //
    std::ofstream cities ("city_coordinates"+filename);
    
    for (int i=0; i < n_city; ++i){     // city initialisation

#ifdef CIRCLE                           // generates cities on unitary circle
        foo_bar = 2*M_PI*rnd.Rannyu();  // angle for position on circle
        
        foo_bar_ville.x = cos(foo_bar); // city added to country
        foo_bar_ville.y = sin(foo_bar);
#endif
        
#ifndef CIRCLE                          // generates cities inside unit square
        foo_bar_ville.x = (rand()/(double)RAND_MAX);
        foo_bar_ville.y = (rand()/(double)RAND_MAX);
#endif
        foo_bar_ville.index = i;
        country.push_back(foo_bar_ville);
        
        cities << std::setw(3) << foo_bar_ville.index << std::setw(12) << foo_bar_ville.x << std::setw(12) << foo_bar_ville.y << std::endl;
        
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
    
    //                                                                  //
    // genrates N different random routes as permutations of the first  //
    //                                                                  //
    if (route.GoodChromo())               // checks if first chromosome
        route_atlas.AddChromosome(route); // has right requirements
    else{                                 // otherwise program stops
        std::cout << "ERROR: First route does not follow requirements, terminating program" << std::endl;
        return -1;
    }
    SetWeight(route);
    
    goodness = false;
    good_check = 0;
    for (i=1; i<N; ++i){        // one chromosome already added
        while (!goodness){
            for (j=0; j<n_city; ++j){
                foo = (int)((n_city-1)*rnd.Rannyu()) + 1; // first city fixed
                bar = (int)((n_city-1)*rnd.Rannyu()) + 1;

                route.SinglePermutation(foo, bar); // shuffles elements
            }
            goodness = route.GoodChromo();
            good_check++;
            if (good_check == attempt_limit){
                std::cout << "ERROR: Max attempt limit in generating different routes has been reached. Terminating program." << std::endl;
                return -1;
            }
        }
        goodness = false;
        SetWeight(route);                     // updates chromosome's weight
        route_atlas.AddChromosome(route);
    }
    route_atlas.Sort();
    
    //                                                                  //
    // Starts evolution of population by new consequent generations     //
    //                                                                  //

    std::ofstream Out ("fitness_comparison"+filename);
    
    for (i = 0; i < generation_no; ++i){
        route_atlas.RandomMutation(N, s_swap_p, m_swap_p, invert_p);
        for (j=0; j < N; ++j)                       // mutations need to be
            SetWeight (route_atlas.population[j]);  // weighted
            
        route_atlas.NextGen(gen_p);
        
        route_atlas.ResetPop();                     // replaces and sorts
        for (j=0; j < N; ++j){                      // new generation
            SetWeight (route_atlas.store_pop[j]);
            route_atlas.AddChromosome(route_atlas.store_pop[j]);
        }
        route_atlas.Sort();
        
#ifndef Get_Params
        route = route_atlas.GetBest();  // route used as placeholder
        Out << std::setw(12) << i       // records fitness evolution
        << std::setw(12) << route.weight
        << std::setw(12) << route_atlas.AvgWeight()
        << std::endl;
#endif
    }
    Out.close();
    
    // prints best final route
    std::ofstream BestRoute ("best_route"+filename);
    route_atlas.Sort();
    route = route_atlas.GetBest();
    route.PrintSequence(BestRoute);
    BestRoute.close();

#ifdef Get_Params
    params << std::setw(12) << s_swap_p
           << std::setw(12) << m_swap_p
           << std::setw(12) << invert_p
           << std::setw(12) << route.weight
           << std::endl;
    params.close();
#endif
    
    return 0;
}
