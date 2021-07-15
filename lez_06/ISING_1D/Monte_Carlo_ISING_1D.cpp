/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"
using namespace std;

//#define PRINT                 // Prints all info to screen
#define RAPID                   // gets indications from command line
                                // instead of input file
#define ACC_RATE

int main(int argc, char *argv[])
{
    arg = argc;
    Input(argv); //Inizialization
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
            Move(metro);
            
#ifndef ACC_RATE
            Measure();
            Accumulate();
#endif
            
        }
        Averages(iblk);   // Print results for current block
    }
    ConfFinal();          // Write final configuration
    
#ifdef ACC_RATE
    if (measure_mode)
        Final_Remarks(); // Prints final characteristic values
                             // depending on temp and h
#endif
    return 0;
}


void Input(char* argv[])
{
    ifstream ReadInput;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    //Read input informations
    ReadInput.open("input.dat");
    
    ReadInput >> temp;
    
    ReadInput >> nspin;

    ReadInput >> J;

    ReadInput >> h;
    
    ReadInput >> metro; // if = 1 Metropolis else Gibbs
    
    ReadInput >> nblk;

    ReadInput >> nstep;
    ReadInput >> restart_bool;   // decides if Metropolis or Gibbs
    ReadInput >> measure_mode;   // decides whether to save final results
    ReadInput >> bias;
    ReadInput >> therm_no;
    ReadInput.close();
    
#ifdef RAPID
    temp = 0.5 + 1.5 * atof(argv[1])/atof(argv[2]);
    
    restart_bool = atoi(argv[3]);
    measure_mode = atoi(argv[4]);
    
    metro = atoi(argv[6]);
    h = atof(argv[5]);
#endif
    
    temp_name = to_string(temp);
    beta = 1.0/temp;
    
    if (metro == 1)
        m_type = "Metropolis";
    else
        m_type = "Gibbs";
    
        cout << "Temperature = " << temp << endl;
    
#ifdef PRINT
        if(metro==1) cout << "The program performs Metropolis moves" << endl;
        else cout << "The program performs Gibbs moves" << endl;
    
        cout << "Classic 1D Ising model             " << endl;
        cout << "Monte Carlo simulation             " << endl << endl;
        cout << "Nearest neighbour interaction      " << endl << endl;
        cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
        cout << "The program uses k_B=1 and mu_B=1 units " << endl;
        
        cout << "Number of spins = " << nspin << endl;
        cout << "Exchange interaction = " << J << endl;
        cout << "External field = " << h << endl << endl;

        cout << "Number of blocks = " << nblk << endl;
        cout << "Number of steps in one block = " << nstep << endl << endl;
    
        if (measure_mode)
            std::cout << "Measure mode activated!" << std::endl;
#endif
    
    // Prepare arrays for measurements
    // Specific variable index
    iu = 0; // Energy
    ic = 1; // Heat capacity
    im = 2; // Magnetization
    ix = 3; // Magnetic susceptibility

    n_props = 4; //Number of observables
    
    std::vector <double> foo_bar;   // initialises the correct
    for (int i=0; i<n_props; ++i){   // no of storage vectors
        block_avg_vector.push_back(foo_bar);
        avg_vector.push_back(foo_bar);
    }

    //initial configuration
    if(restart_bool==1){
        // gets previous configuration (if it exists)
        //std::cout << "Restarting from previous configuration" << std::endl;
        ReadInput.open("config.final");
        for (int i=0; i<nspin; ++i)
            ReadInput >> s[i];
    }
    else{
        // if there isn't a previous configuration a random
        // one is generated
        double foo;
    
        for (int i=0; i<nspin; ++i){
            foo = rnd.Rannyu();
            if( foo >= 0.5) s[i] = 1;
            else s[i] = -1;
        }
    }

    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    //cout << "Initial energy = " << block_avg_vector[iu][0]/(double)nspin << endl;
    accepted=0;
    attempted=0;
    total_accepted = 0;
}


void Move(int metro)    // single Monte Carlo step
{
    int o,sm;
    double p, energy_old, energy_new;
    double energy_up, energy_down;
    double prob_diff, accept_prob;
    for(int i=0; i<nspin; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd.Rannyu()*nspin);

        if(metro==1) //Metropolis
        {
            sm = -1 * s[o];         // flips spin
            energy_old = Boltzmann(s[o],o);
            energy_new = Boltzmann(sm,o);
            prob_diff = exp(-beta * (energy_new - energy_old));
            
            if (prob_diff < 1)
                accept_prob = prob_diff;
            else
                accept_prob = 1;
            
            p = rnd.Rannyu();
            if (p < accept_prob){    // changed only if move is
                s[o] = sm;           // accepted
                accepted++;
                total_accepted++;
            }
        }
        else{ //Gibbs sampling
            energy_up = Boltzmann(1,o);
            energy_down = Boltzmann(-1,o);
            prob_diff = exp(beta * (energy_up - energy_down));
            // probability of passing from down to up
            sm = s[o];
            accept_prob = 1./(1+prob_diff);
            p = rnd.Rannyu();
            
            if (p < accept_prob)    // probability of passing from
                s[o] = 1;           // down to up accepted
            else
                s[o] = -1;          // otherwise it's assigned down
            
            if(sm != s[o]){         // counts the times spin is
                accepted++;         // actually flipped
                total_accepted++;
            }
        }
        attempted++;
    }
}


double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure()
{
    double m = 0.0, u = 0.0;
    
    //cycle over spins
    for (int i=0; i<nspin; ++i)
    {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        m += s[i];
    }
    
    // walker saves value of energy and magnetisation
    walker[iu] = u;
    walker[ic] = u*u;
    walker[im] = m;
    walker[ix] = m*m;
}

void Accumulate(){
    for (int i=0; i < n_props; i++)
        block_avg_vector[i].push_back(walker[i]);
    // average vectors increased with energy and magnetisation

}

void Reset(int iblk) //Reset block averages
{
    for(int i=0; i < n_props; i++){
        walker[i] = 0;              // in order to avoid overflow
        block_avg_vector[i].clear();
        if (iblk==1)
            avg_vector[i].clear();
        // initialises at the beginning also the full average vector (just to be sure)
    }
    
    blk_norm = 0;
    accepted = 0;
    attempted = 0;
}


void Averages(int iblk) //Print results for current block
{
    ofstream Ene, Heat, Mag, Chi;
   
    // total lattice quantities
    stima_u = average_vec(block_avg_vector[iu]);
    stima_x = average_vec(block_avg_vector[ic]);
    stima_c = beta*beta*(stima_x - stima_u*stima_u)/(double)nspin;
    stima_u = stima_u /(double)nspin;   // redefined as en per spin
    stima_m = average_vec(block_avg_vector[im]);
    stima_chi = beta*average_vec(block_avg_vector[ix]);   // <M> = 0 fixed (for h=0)
    
#ifdef PRINT
    if(arg < 5){                    // ignored in rapid simulations
        cout << "Block number " << iblk << endl;
        if (metro==1)
            cout << "Acceptance rate " << accepted/attempted << endl << endl;
        
        std::cout << "Average spin energy: "
                  << stima_u
                  << "; specific heat: "
                  << stima_c
                  << "; magnetisation: "
                  << stima_m /(double)nspin
                  << "; magnetic susceptibility: "
                  << stima_chi/(double)nspin << std::endl;
        cout << "----------------------------" << endl << endl;
    }
#endif
    
    // the following estimates the progressive average, considering
    // the previous blocks, and saves the values in external file
    // for each characteristic quantity
    avg_vector[iu].push_back(stima_u);
    avg_vector[ic].push_back(stima_c);
    avg_vector[ix].push_back(stima_chi/(double)nspin);
    avg_vector[im].push_back(stima_m/(double)nspin);

    if (!measure_mode){   // ignored in measure mode
        if ( h == 0.){
        // Energy
        Ene.open("data_files/output.energy_"+temp_name+"_"+m_type+".txt", ios::app);
        Ene << iblk << std::setw(12) << stima_u << std::setw(12);
        block_average (avg_vector[iu], Ene);
        // estimates progressive average with error and prints
        Ene.close();
            
        // Specific heat
        Heat.open("data_files/output.sp_heat_"+temp_name+"_"+m_type+".txt", ios::app);
        Heat << iblk << std::setw(12) << stima_c << std::setw(12);
            block_average (avg_vector[ic], Heat);
        Heat.close();
            
        // Susceptibility
        Chi.open("data_files/output.chi_"+temp_name+"_"+m_type+".txt", ios::app);
        Heat << iblk << std::setw(12) << stima_chi << std::setw(12);
        block_average (avg_vector[ix], Chi);
        Chi.close();
        }
        else {
        // Magnetisation
        Mag.open("data_files/output.mag_"+temp_name+"_"+m_type+".txt", ios::app);
        Mag << iblk << std::setw(12) << stima_m << std::setw(12);
            block_average (avg_vector[im], Mag);
        Mag.close();
        }
    }
}


void ConfFinal(void)
{
    double total_attempts = nstep*nblk*nspin;
    
    ofstream WriteConf;
    if (arg < 5){
        std::cout << "Acceptance rate: " << total_accepted /total_attempts << std::endl;

        std::cout << "Print final configuration to file config.final " << endl << endl;
    }
    
    WriteConf.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
        WriteConf << s[i] << std::endl;
    }
    WriteConf.close();
    
#ifdef ACC_RATE
    std::ofstream AR;
    AR.open("acceptance_rate."+m_type+".txt", std::ios::app);
    AR << temp << std::setw(12) << total_accepted /total_attempts << std::endl;
    AR.close();
#endif
    
    rnd.SaveSeed();
    
}

void Final_Remarks(void){
    // prints final results of energy, specific heat and magnetic
    // susceptibility (for external field h=0)
    std::ofstream Out;
    if (h==0.){
        Out.open("data_files/temp_compare_"+m_type+".txt",ios::app);
        Out << temp << std::setw(12)
            << average_vec(avg_vector[iu]) << std::setw(12)
            << error_vec(avg_vector[iu]) << std::setw(12)
            << average_vec(avg_vector[ic]) << std::setw(12)
            << error_vec(avg_vector[ic]) << std::setw(12)
            << average_vec(avg_vector[ix]) << std::setw(12)
            << error_vec(avg_vector[ix]) << std::endl;
    }
    else{
        Out.open("data_files/temp_mag_"+m_type+".txt",ios::app);
        Out << temp << std::setw(12)
        << average_vec(avg_vector[im]) << std::setw(12)
        << error_vec(avg_vector[im]) << std::endl;
    }
    
    Out.close();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
