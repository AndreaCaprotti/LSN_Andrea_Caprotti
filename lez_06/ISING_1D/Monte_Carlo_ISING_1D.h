/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

//Random numbers
#include <string>
#include "../../libraries/block_stat.h"
#include "../../libraries/random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];

// averages
double blk_norm;
double stima_u,stima_c,stima_m,stima_x,stima_chi;

// vectors that stores a number of vectors for each relevant quantity
std::vector<std::vector<double>> block_avg_vector;
std::vector<std::vector<double>> avg_vector;

//configuration
int arg;
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;
std::string temp_name;

// simulation
int nstep, nblk, metro;
double accepted, total_accepted, attempted;
std::string m_type;
bool restart_bool, measure_mode;
int therm_no;                   // number of therm steps
double bias;                    // bias of spin for low temperatures

//functions
void Input(char**);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Final_Remarks(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
