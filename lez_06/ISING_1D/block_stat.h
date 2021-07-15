//
//  block_stat.h
//  functions to handle vectors of data in order to get a
//  statistics on blocks instead of single values
//
//  Created by Andrea Caprotti on 06/04/21.
//

#ifndef block_stat_h
#define block_stat_h

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <array>
#include <vector>
#include <iomanip>

// Average of elements of std::vector
double average_vec (std::vector<double> &vec){
    double sum = 0;
    int N = vec.size();
    int i;              // cycle index
    for (i=0; i<N; i++)
        sum += vec[i];
    
    if (N == 0)
        return 0;
    else
        return sum/N;
}

// Standard deviation computed from average
double std_dev_vec(std::vector<double> &vec){
    int N = vec.size();
    double avg, avg_2;
    std::vector<double> vec_2;
    int i;
    
    for (i=0; i<N; i++)
        vec_2.push_back(vec[i]*vec[i]);
    
    avg = average_vec(vec);
    avg_2 = average_vec(vec_2);
    return avg_2 - avg*avg;
}

// Error (as std.dev of the average) of elements of std::vector
double error_vec (std::vector<double> &vec){
    int N = vec.size();
    double std_dev;
    
    std_dev = std_dev_vec(vec);
    
    if (N != 1)
        return sqrt(std_dev / (N-1)); // std dev of average
    else
        return 0;
}

// average for vector of averages of single block
void block_average (std::vector<double> &vec, std::ofstream &datafile){
    datafile << average_vec(vec) << std::setw(12) << error_vec(vec) << std::endl;
    return;
}

// block average for vector of raw data
void full_block_average(std::vector<double> &vec, std::ofstream &datafile, int block_no, bool progressive_avg){
    
    int L = int(vec.size()/block_no);
    int k=0;
    int i,j;                              // cycle index
    std::vector<double> partial_avg;      // separates data in blocks
    std::vector<double> block_avg;        // used to evaluate full average
    
    block_avg.clear();
    for (i=0; i<block_no; i++){
        partial_avg.clear();
        for (j=0; j<L;j++){
            partial_avg.push_back(vec[k]);
            k++;
        }
        block_avg.push_back(average_vec(partial_avg));
        // in order to not complicate my life too much,
        // this flag makes activates the progressive computation
        if (progressive_avg)
            block_average(block_avg, datafile);
    }

    // to ensure last value doesn't appear twice
    if (!progressive_avg)
        block_average(block_avg, datafile);

    return;
}

#endif /* block_stat_h */
