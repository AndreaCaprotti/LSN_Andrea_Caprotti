//
//      Let's start over.
//

#ifndef __NVE__
#define __NVE__


//parameters, observables
const int m_props=1000;
const int wd=12;
int n_props, arg;
int iv,ik,it,ie,ip,igofr;
int ind, jndex;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double walker[m_props];

std::string sys_type;
int sys;

// averages
int blk_len;
double acc,att;
double block_avg[m_props];
double glob_av[m_props];
double glob_av2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double fs;                      // scale factor

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;
double v_tail, p_tail;

// simulation
int nstep, nblocks;
int iprint, seed, istep;
double delta;
bool measure_mode, restart;

//g(r) calculation
int nbins;
double bin_size,norm;

//functions
void Input(char**);
void FirstStart(void);
void NextStart(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void Averages(int);
void Reset(int);
void GofR(void);
double Force(int, int);
double Pbc(double);
double Error (double, double, int);

#endif
/****************************************************************
*****************************************************************/
