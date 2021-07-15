//
//  01_3.cpp
//  Simulates Buffon's needle-dropping experiment to estimate \pi
//
//  Created by Andrea Caprotti on 23/03/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "../../libraries/random.h"
#include "../../libraries/block_stat.h"

int main (int argc, char *argv[]){
    
    Random rnd;
    
    int M = 100000;               // Total number of random throws
    int N = 100;                  // Number of blocks
    int N_throws = int(M/N);      // Number of throws per block
    int N_hit;                    // Number of hits per block
    double d = 1;                 // Distance between lines
    double needle_L = 0.5;        // Needle lenght
    double needle_x,needle_y;     // Needle extremes
    
    bool angle_generator = false; // all this to generate second
    double x_angle, y_angle;      // needle head
    double r_angle;
    
    int i,j;                        // cycle index
    
    // estimates cumulative average on increasing number of blocks
    std::vector<double> avg;
    
    std::ofstream avg_data ("pi_extimation.txt");
    
    for (i = 0; i < N; i++) {
        N_hit = 0;
        for (j = 0; j < N_throws; j++){
            needle_x = d * rnd.Rannyu(); // number between 0 and d
            
            angle_generator=false;       // generates second pos
            while (!angle_generator){    // following sin distrib
                                         // (without pi)
                x_angle = rnd.Rannyu() - 0.5;
                y_angle = rnd.Rannyu() - 0.5;
                r_angle = sqrt(pow(x_angle,2)+pow(y_angle,2));
                if (r_angle<0.5)
                    angle_generator = true;
            }
            needle_y = needle_L * y_angle / r_angle;
            
            if (needle_x+needle_y >= d or needle_x+needle_y <= 0)
                N_hit++;
        }
        avg.push_back((2*needle_L*N_throws)/(N_hit*d));
        block_average(avg,avg_data);
    }
    avg_data.close();
    
    return 0;
}
