//
//  10_2.hpp
//
//
//  Created by Andrea Caprotti on 03/07/21.
//

#ifndef _10_2_hpp
#define _10_2_hpp

#include "random.h"
#include "chromosome.h"

struct City{
    double x,y;     // coordinates
    int index;      // to recognise particular city
};

// city set
const int n_city = 32 ;
double len_matrix[n_city][n_city] = {{0}};
std::vector <City> country;

// mutation probabilities
double s_swap_p = 0.05;
double m_swap_p = 0.05;
double invert_p = 0.05;
double gen_p = 0.5;
const double prob_exp = 3;

// chromosomes
PopulationChrome route_atlas(n_city, prob_exp);
int attempt_limit = 1000;
int generation_no = 200;
int N = 900;                        // number of possible routes
int n_migr = 20;                    // determines migration frequency of best
int migrations = 2;                     // element to another continent

void SetWeight (Chromosome &chromo){
    double w = 0;
    int index_1, index_2;
    
    for (int i=1; i < n_city; ++i){
        index_1 = chromo.gene_sequence[i];
        index_2 = chromo.gene_sequence[i-1];
        
        w+= len_matrix[index_1][index_2];
    }
    w+= len_matrix[index_1][0];         // closes the route
    chromo.weight = w;
    
    return;
}

#endif /* _10_2_hpp */
