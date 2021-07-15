//
//  random_walk.cpp
//  Class to describe position evolution of a D-dimension random walk
//  with possible discrete or continuous steps
//
//  Expanded to describe a Metropolis walk (adds the possibility of
//  rejecting the move). Continuous step only works in 3D, though,
//  sinnce the spherical coordinates conversion cannot be easily
//  generalised…
//
//  Created by Andrea Caprotti on 31/03/21.
//

#ifndef __R_W__
#define __R_W__

#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

class Random_Walk {
// MD random walk
// this class actually simply implements a generci walk, the
// randomness needs to arrive from the outside because of
// problems reguarding vectors of these walks…
// not the most elegant solution, but the fastest which works
private:
    std::vector<double> pos;            // position vector (dim
                                        // not fixed a priori)
    int dim;                            // number of dimensions
    
public:
    Random_Walk(){};
    ~Random_Walk(){};                   // destructor
    
    void start(int N) {                 // initialises a walk
        dim = N;
        pos.clear();
        for (int i=0; i < dim; i++)     // default position is origin
            pos.push_back(0);
        return;
    }
    
    double vect_length(){               // return square vector
                                        // distance from origin
        double r_2 = 0;
        int i;                          // cycle index
        for (i = 0; i<dim; i++)
            r_2 += pow(pos[i],2);
        return r_2;
    }
    
    void discrete_step_increment (double sl, double random_1, double random_2, int direction){
    // receives step length, random variables in [0,1) from
    // outside and verse of step (+/- 1)
        double verse;
        if (direction > dim)          // eventually dir can be fixed
            direction = (int) (dim*random_1);
        
        verse = random_2;
        if (verse < 0.5)
            pos[direction] -= sl;
        else if (verse > 0.5)
            pos[direction] += sl;
        
        return;
    }
    
    void continuous_step_increment (double sl, double ran_theta, double ran_phi){
        // spherical coords to cartesian (for 3D walk)
        // params are scaled so that phi\in[0,2pi] and theta\in[0,pi]
        pos[0] += sl*cos(2*M_PI*ran_phi)*sin(M_PI*ran_theta);
        pos[1] += sl*sin(2*M_PI*ran_phi)*sin(M_PI*ran_theta);
        pos[2] += sl*cos(2*M_PI*ran_phi);
        
        return;
    }
    
    double single_dim_step(double fs, double random){
        double step_l;
        step_l = fs*2*(random - 0.5);    // generates step in either dir
        pos[0]+= step_l;
        return vect_length();
    }
    
    // direct access to position
    void print_position (std::ofstream& data_file){
        for (int i = 0; i < dim; i++ )
            data_file << pos[i] << std::setw(15);
        data_file << std::endl;
        return;
    }
    
    double current_position(int k){
        // returns single coordinate in k-th direction
        // since position is obstinately a private member
        return pos[k];
    }
    
    void assign_position (std::vector<double> coords){
        // simply substitutes entirely vector of coordinates
        // (with additional check that size is the same expected)
        if (coords.size() == dim)
            pos = coords;
        else
            std::cout << "Error in position assignment: size does not correspond" << std::endl;
        return;
    }
};


// Metropolis class
class Metropolis_Walk : public Random , public Random_Walk{

private:
    int N;                                  // degrees of freedom
    std::vector<double> save_position;      // nomen omen
    double conf_probability;                // saves conf probability
    
    double accept_bound, attempt_probabilty;
    // defined here so that it's not defined each time is called
    
public:
    Metropolis_Walk(std::vector<double> start_pos) { // constructor
        N = start_pos.size();           // assings degrees of fr.
        start(N);
        
        assign_position(start_pos);
        save_current_position();
        no_acc = 0;
        return;
    };
    ~Metropolis_Walk(){};               // destructor
    
    void new_move (double fs){          // step lenght given
        save_current_position();

        if (N==1)
            single_dim_step(fs, Rannyu());
        else
            continuous_step_increment(fs*Rannyu(), Rannyu(), Rannyu());
            // evolves position by passing two random values
        
        return;
    }
    
    void new_gauss_move(double fs){
        save_current_position();
        
        double step_l;                // gaussian step corresponds to
                                      // movement in all three dirs
        for (int i=0; i<N; i++){
            step_l = fs*Gauss(0,1);
            discrete_step_increment(step_l,Rannyu(), Rannyu(),i);
        }
        return;
    }
    
    void accept_move (double new_prob){
    // simple example for symmetric transition probability
    // nevertheless, quite general, does not assume in any way prob.
        
        if (new_prob < conf_probability)
            accept_bound = new_prob/conf_probability;
        else
            accept_bound = 1;
        
        attempt_probabilty = Rannyu();
        if (attempt_probabilty < accept_bound){
            assign_prob(new_prob);             // move accepted
            no_acc++;

        }
        else
            assign_position(save_position);
            // otherwise previous position is restored
            // and prob isn't changed
        
        return;
    }
    
    // probability distribution to be followed is assigned from the
    // outside, therefore it can be reused any time
    void assign_prob(double new_prob){
        conf_probability = new_prob;
        return;
    }
    
    void save_current_position(){
        save_position.clear();
        for (int i=0; i<N; i++)
            save_position.push_back(current_position(i));
        return;
    }
    
    int no_acc;                     // number of accepted attempts
};

#endif
