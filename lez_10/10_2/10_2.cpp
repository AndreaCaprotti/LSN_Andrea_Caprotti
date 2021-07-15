//
//  10_2.cpp
//  Variation of 09_1.cpp which includes MPI exchange
//
//  Created by Andrea Caprotti on 03/07/21.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>                    // M_PI (as in \pi…)
#include <string>
#include "mpi.h"                    // as in parallel program…
#include "10_2.hpp"

#define MIGRATE

using namespace std;

int main (int argc, char* argv[]){
    
    int foo, bar;
    int i, j, k, good_check;            // cycle index
    bool goodness;                      // check for requirements
    
    int numtasks, rank, tag=1;     // MPI initialisation
    int rank_send, next, prev;
    int buffer[n_city], receive[n_city];

    MPI_Request reqs[2];   // required variable for non-blocking calls

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for (i=0; i < rank +1; ++i)     // starting p1 p2 change w/ rank
            Primes >> p1 >> p2 ;
    }
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    srand(0);                       // this should be the same for all!

    City foo_bar_ville;
    Chromosome route(n_city);
    std::string filename;
    
    filename = to_string(rank);         // to distinguish the various
#ifdef MIGRATE
    filename += "_mig_"+to_string(migrations);
#endif
    //                                                   //
    // random cities generated on circle of radius r = 1 //
    // (or in a square of side l = 1)                    //
    // assings first chromosome                          //
    //                                                   //
    if (rank == 0){
        std::ofstream cities ("city_coordinates.txt");
        
        for (int i=0; i < n_city; ++i){     // city initialisation
            // generates cities inside unit square
            foo_bar_ville.x = (rand()/(double)RAND_MAX);
            foo_bar_ville.y = (rand()/(double)RAND_MAX);
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
    }
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
    //   Starts evolution of population by new consequent generations   //
    //                                                                  //

    std::ofstream Out ("fitness_comparison_"+filename+".txt");

    for (i = 0; i < generation_no; ++i){
        route_atlas.RandomMutation(N, s_swap_p, m_swap_p, invert_p);
#ifdef MIGRATE
        if ((i+10) % n_migr == 0){
            for (k = 0; k < migrations; ++k){
                route = route_atlas.GetBest();
                
                for (j = 0; j < n_city; ++j){    // copies best route to share
                    buffer[j]=route.gene_sequence[j];
                    receive[j] = 0;
                }
                
                if (rank == 0)  // determines where best is sent
                    rank_send = (int)(3*rnd.Rannyu())+1; // 1, 2 or 3 shift
                
                MPI_Bcast(&rank_send, 1, MPI_INT, 0, MPI_COMM_WORLD);
                // info on ranks shared w/ everyone
                prev = (4 + rank - rank_send) % 4;
                next = (rank + rank_send) % 4;
                
                // post non-blocking receives and sends for neighbors
                MPI_Irecv(&receive, n_city, MPI_INT, prev, tag, MPI_COMM_WORLD, &reqs[0]);
                
                MPI_Isend(&buffer, n_city, MPI_INT, next, tag, MPI_COMM_WORLD, &reqs[1]);
                // wait for all non-blocking operations to complete
                MPI_Barrier(MPI_COMM_WORLD);
                
                // adds gained routes from tail
                route.gene_sequence.clear();
                for (int j=0; j<n_city;++j)
                    route.gene_sequence.push_back(receive[j]);
                if (route.GoodChromo())
                    route_atlas.ReplaceChromo(route, N-k-1);
            }
        }
#endif
        for (j=0; j < N; ++j)                       // mutations need to be
            SetWeight (route_atlas.population[j]);  // weighted
            
        route_atlas.NextGen(gen_p);
        
        route_atlas.ResetPop();                     // replaces and sorts
        for (j=0; j < N; ++j){                      // new generation
            SetWeight (route_atlas.store_pop[j]);
            route_atlas.AddChromosome(route_atlas.store_pop[j]);
        }
        route_atlas.Sort();
        
        route = route_atlas.GetBest();  // route used as placeholder
        Out << std::setw(12) << i       // records fitness evolution
        << std::setw(12) << route.weight
        << std::setw(12) << route_atlas.AvgWeight()
        << std::endl;

    }
    Out.close();
    
    // prints best final route
    std::ofstream BestRoute ("best_route"+filename+".txt");
    route_atlas.Sort();
    route = route_atlas.GetBest();
    route.PrintSequence(BestRoute);
    BestRoute.close();
    
    MPI_Finalize();
    
    return 0;
}
