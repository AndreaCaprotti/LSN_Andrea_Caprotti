//
//  random_walk.cpp
//  Class to describe position evolution of a 3D random walk
//  with possible discrete or continuous steps
//
//  Expanded to describe a Metropolis walk (adds the possibility of
//  rejecting the move). Only works in 3D, though
//
//  Created by Andrea Caprotti on 31/03/21.
//

#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"

class Random_Walk : public Random {
// MD random walk
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
            r_2 +=pow(pos[i],2);
        return r_2;
    }
    
    double discrete_step_increment (double sl, int direction){
                                        // receives step length
        double verse;                   // verse of step (+/- 1)
        
        verse = Rannyu();
        if (verse < 0.5)
            pos[direction] -= sl;
        else if (verse > 0.5)
            pos[direction] += sl;
        
        return vect_length();
    }
    
    double continuous_step_increment (double sl){
        double phi;                  // azimuthal angle
        double theta;                // polar angle
        
        phi = Rannyu(0, 2*M_PI);
        theta = Rannyu(0, M_PI);
        
        pos[0] += sl*cos(phi)*sin(theta); // spherical coords to
        pos[1] += sl*sin(phi)*sin(theta); // cartesian (for 3D walk)
        pos[2] += sl*cos(phi);
        
        return vect_length();
    }
    
    double single_dim_step(double fs){
        double step_l;
        step_l = fs*Rannyu(-1,1);    // generates step in either dir
        pos[0]+= step_l;
        return vect_length();
    }
    
    // accesses to position directly
    void print_position (std::ofstream& data_file){
        for (int i = 0; i < dim; i++ )
            data_file << pos[i] << std::setw(12);
        data_file << std::endl;
        return;
    }
    
    double current_position(int k){
        // returns single coordinate in k-th direction
        // since position is obstinately a private member
        return pos[k];
    }
    
    void assign_position (std::vector<double> coords){
        pos = coords;
        return;
    }
};


// Metropolis class
class Metropolis_Walk : public Random_Walk {

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
    ~Metropolis_Walk(){};                   // destructor
    
    void new_move (double fs){          // step lenght given
                                        // (generated outside)
        double step_l;                  // random.h methods only in
        save_current_position();        // class

        if (N==1)
            single_dim_step(fs);
        else {
        step_l = fs*Rannyu();
        continuous_step_increment(step_l);  // evolves position
        }
        return;
    }
    
    void new_gauss_move(double fs){
        save_current_position();
        
        double step_l;                // gaussian step corresponds to
                                      // movement in all three dirs
        for (int i=0; i<N; i++){
            step_l = fs*Gauss(0,1);
            discrete_step_increment(step_l,i);
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

