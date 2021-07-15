//
//  10_1.hpp
//  
//
//  Created by Massimo Caprotti on 30/06/21.
//

#ifndef _0_1_hpp
#define _0_1_hpp

#include <vector>
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

double SetWeight (Chromosome &chromo){
    double w = 0;
    int index_1, index_2;
    
    for (int i=1; i < n_city; ++i){
        index_1 = chromo.gene_sequence[i];
        index_2 = chromo.gene_sequence[i-1];
        
        w+= len_matrix[index_1][index_2];
    }
    w+= len_matrix[index_1][0];         // closes the route
    
    return w;
}

// simulation
int n_steps = 900;  // number of steps in simulation
double beta=1;
int accepted = 0;

#endif /* _0_1_hpp */
