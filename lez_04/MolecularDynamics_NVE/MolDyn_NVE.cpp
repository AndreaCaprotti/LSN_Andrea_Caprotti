//
//      Let's start over.
//

#include "../../libraries/random.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>            // rint, pow
#include "MolDyn_NVE.h"

//#define PRINT
#define G_of_R              // g(r) computed only when needed
                            // alternative to other output

int main(int argc, char* argv[]){
    arg = argc;
    Input(argv);             // Inizialization
    int nconf = 1;
    for (int iblk = 1; iblk <= nblocks; ++iblk){
        Reset(iblk);
        for(istep=1; istep <= nstep; ++istep){
            Move();           //Move particles with Verlet algorithm
            if(istep%10 == 0){
                Measure();     //Properties measurement
        //      ConfXYZ(nconf);//Write actual configuration in XYZ format
                nconf += 1;
            }
        }
#ifdef PRINT
        if (iblk*nstep % iprint == 0)
            std::cout << "Number of time-steps: " << iblk * nstep << std::endl;
#endif
        if (measure_mode){
            Averages(iblk);
        }
    }
    
    ConfFinal();         //Write final configuration to restart

    return 0;
}


void Input(char* argv[]){ //Prepare all stuff for the simulation
    std::ifstream ReadInput,ReadConf;
    
#ifdef PRINT
    std::cout << "Classic Lennard-Jones fluid        " << std::endl;
    std::cout << "Molecular dynamics simulation in NVE ensemble  " << std::endl << std::endl;
    std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << std::endl << std::endl;
    std::cout << "The program uses Lennard-Jones units " << std::endl;
#endif
    
    ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;

    ReadInput >> npart;

    ReadInput >> rho;
    vol = (double)npart/rho;
    box = pow(vol,1.0/3.0);
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;         // total block lenght
    ReadInput >> iprint;
    ReadInput >> nblocks;
    ReadInput >> restart;
    ReadInput >> measure_mode;
    ReadInput >> sys;
    
    if (arg > 1){       // passed by command line for long sims
        restart = atoi(argv[1]);
        measure_mode = atoi (argv[2]);
        sys = atoi (argv[3]);
    }
    if (measure_mode)
        std::cout << "Measure mode activated!" << std::endl;
    
    switch (sys) {              // determines file name
        case 0:
            sys_type = "solid";
            break;
        case 1:
            sys_type = "liquid";
            break;
        case 2:
            sys_type = "gas";
            break;
        default:
            sys_type = "_";
            break;
    }
    blk_len = 0;      // sample block lenght

#ifdef PRINT
    std::cout << "Number of particles = " << npart << std::endl;

    std::cout << "Density of particles = " << rho << std::endl;

    std::cout << "Volume of the simulation box = " << vol << std::endl;

    std::cout << "Edge of the simulation box = " << box << std::endl;

    std::cout << "The system is considered to be a " << sys_type << std::endl;
    
    std::cout << "The program integrates Newton equations with the Verlet method " << std::endl;
    std::cout << "Time step = " << delta << std::endl;
    std::cout << "Number of steps = " << nstep << std::endl << std::endl;
#endif
    
    ReadInput.close();

    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    ip = 4; //Pressure
    igofr = 5; //Number of observables

    nbins = 100;
    n_props = igofr + nbins;
    bin_size = (box/2.0)/(double)nbins;
    norm = 4*M_PI*rho*(double)npart/3; // g(r) normalisation
    
    // tail corrections for potential energy and pressure
    v_tail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
    p_tail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
    
    // chooses which starting configuration to follow
    if (restart)
        NextStart();
    else
        FirstStart();
    
    if(!measure_mode){               // if already at temperature
        for (int i=0; i<npart; ++i){ // rscaling in not impellent
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    }
    return;
}

void FirstStart(void){                // First-ever initialisation
    std::ifstream ReadConf;

#ifdef PRINT
    std::cout << "Read initial configuration from file config.0 " << std::endl << std::endl;
    std::cout << "Prepare random velocities with center of mass velocity equal to zero " << std::endl << std::endl;
#endif
    
    //Read initial configuration
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();
        
    //Prepare (random) initial velocities
    double sumv[3] = {0.0, 0.0, 0.0};
    Random rnd; // better random number generator instead of std lib
        
    for (int i=0; i<npart; ++i){
        vx[i] = rnd.Rannyu() - 0.5;
        vy[i] = rnd.Rannyu() - 0.5;
        vz[i] = rnd.Rannyu() - 0.5;
        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
    }
        
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0;
    for (int i=0; i<npart; ++i){
        vx[i] = vx[i] - sumv[0];
        vy[i] = vy[i] - sumv[1];
        vz[i] = vz[i] - sumv[2];
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    return;
}

void NextStart(void){
    // if the bool is true then the old config is loaded
    // in vec_old[] instead of creating a random previous config
    std::ifstream ReadConf;
    
    system ("cp config.final config.0");
    //Read initial configuration
#ifdef PRINT
    std::cout << "Read initial configuration from file config.0 " << std::endl;
    std::cout << "Read previous configuration from file config.old " << std::endl << std::endl;
#endif
    
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i]=x[i]*box;   // coords need to be scaled to box size!
        y[i]=y[i]*box;
        z[i]=z[i]*box;
    }
    ReadConf.close();
    
    // load previous configuration
        ReadConf.open("config.old");
    for (int i=0; i<npart; ++i){
        ReadConf >> xold[i] >> yold[i] >> zold[i];
        
        xold[i]=xold[i]*box;
        yold[i]=yold[i]*box;
        zold[i]=zold[i]*box;
    }
    ReadConf.close();
    
    Move();                      // finds following config
    
    if (measure_mode)            // in measure mode the right
        fs = 1;                  // temperature is already set
    else{                        // rescaling needed only before
        double t=0.0;
        for (int i=0; i<npart; ++i)
            t += 0.5 * ((Pbc(x[i] - xold[i])/delta)*(Pbc(x[i] - xold[i])/delta) + (Pbc(y[i] - yold[i])/delta)*(Pbc(y[i] - yold[i])/delta)  + (Pbc(z[i] - zold[i])/delta)*(Pbc(z[i] - zold[i])/delta));
        // velocity estimated in t+dt/2 instead of using estimation
        // through Move(), as old configuration is a "ghost" conf.
        // needed solely to move first step
        
        t = (2.0 / 3.0) * t/(double)npart; // estimates start temperature
        fs = sqrt(temp/t);               // new velocity scale factor
    }
    return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(ind=0; ind<npart; ++ind){ //Force acting on particle i
    fx[ind] = Force(ind,0);
    fy[ind] = Force(ind,1);
    fz[ind] = Force(ind,2);
  }

  for(ind=0; ind<npart; ++ind){ //Verlet integration scheme

    xnew = Pbc(2.0 * x[ind] - xold[ind] + fx[ind] * pow(delta,2));
    ynew = Pbc(2.0 * y[ind] - yold[ind] + fy[ind] * pow(delta,2));
    znew = Pbc(2.0 * z[ind] - zold[ind] + fz[ind] * pow(delta,2));

    vx[ind] = Pbc(xnew - xold[ind])/(2.0 * delta);
    vy[ind] = Pbc(ynew - yold[ind])/(2.0 * delta);
    vz[ind] = Pbc(znew - zold[ind])/(2.0 * delta);

    xold[ind] = x[ind];
    yold[ind] = y[ind];
    zold[ind] = z[ind];

    x[ind] = xnew;
    y[ind] = ynew;
    z[ind] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }

  return f;
}

void Measure(void){ //Properties measurement
    int bin;
    double v, w, t, vij, wij;
    double dx, dy, dz, dr;
    std::ofstream Epot, Ekin, Etot, Temp, Pres;
  
    v = 0.0; //reset observables
    w = 0.0;
    t = 0.0;
    
    //reset the hystogram of g(r)
    for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

    //cycle over pairs of particles
    for (ind=0; ind < npart-1; ++ind){
        for (jndex = ind+1; jndex<npart; ++jndex){

            dx = Pbc( xold[ind] - xold[jndex] );
            dy = Pbc( yold[ind] - yold[jndex] );
            dz = Pbc( zold[ind] - zold[jndex] );

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);

            // update of the histogram of g(r)
            bin = floor(dr/bin_size) + igofr;
            walker[bin] += 2.;
            if(dr < rcut){
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

                //Potential energy
                v += vij;
                w += wij;
            }
        }
    }

    //Kinetic energy
    for (ind=0; ind < npart; ++ind) t += 0.5 * (vx[ind]*vx[ind] + vy[ind]*vy[ind] + vz[ind]*vz[ind]);
   
    walker[iv] = v/(double)npart + v_tail; //Potential energy per particle
    walker[ik] = t/(double)npart; //Kinetic energy per particle
    walker[it] = (2.0 / 3.0) * t/(double)npart; //Temperature
    walker[ie] = (t+v)/(double)npart; //Total energy per particle
    walker[ip] = (48.0 * w / 3.0 + p_tail * (double)npart)/ vol + rho * walker[it];
    // temperature isn't fixed from the outside but depends on kinetic energy
    
#ifndef G_of_R
    if (!measure_mode){
        Epot.open("output_epot."+sys_type+".txt",std::ios::app);
        Ekin.open("output_ekin."+sys_type+".txt",std::ios::app);
        Temp.open("output_temp."+sys_type+".txt",std::ios::app);
        Etot.open("output_etot."+sys_type+".txt",std::ios::app);
        Pres.open("output_pres."+sys_type+".txt",std::ios::app);
        Epot << walker[iv]  << std::endl;
        Ekin << walker[ik]  << std::endl;
        Temp << walker[it]  << std::endl;
        Etot << walker[ie]  << std::endl;
        Pres << walker[ip]  << std::endl;
        
        Epot.close();
        Ekin.close();
        Temp.close();
        Etot.close();
        Pres.close();
    }
#endif
    
    for (ind = 0; ind < n_props; ++ind){   // accumulates values
        block_avg[ind] += walker[ind];
    }
    blk_len += 1;
    
    return;
}

void Reset(int iblk) {          // Reset block averages
    for(ind=0; ind < n_props; ++ind){
        block_avg[ind] = 0;
        if(iblk == 1){          // initialises global vector
            glob_av[ind] = 0;
            glob_av2[ind] = 0;
        }
    }
    blk_len = 0;
    return;
}

void Averages(int iblk){
#ifndef G_of_R
    std::ofstream Output;
    std::vector<std::string> file_name ={"epot", "ekin","etot","temp","pres"};
    double val, err;
    
    for (ind = 0; ind < igofr; ++ind){
        Output.open("ave_"+file_name[ind]+"."+sys_type+".txt", std::ios::app);
        val = block_avg[ind]/(double)blk_len;
        glob_av[ind] += val;
        glob_av2[ind] += val*val;
        err = Error(glob_av[ind],glob_av2[ind],iblk);
        Output << std::setw(wd) << glob_av[ind]/(double)iblk
               << std::setw(wd) << err << std::endl;
        Output.close();
    }
#endif
#ifdef G_of_R
    double r, gdir, dV, err_gdir;
    std::ofstream Gave;
    Gave.open("g_of_r."+sys_type+".txt",std::ios::app);
    //g(r) progressive average
    for(int i=igofr; i < n_props; ++i){
        r = bin_size/2 + i*bin_size;    // bin representative
        gdir = block_avg[i]/(double)blk_len;
        glob_av[i] += gdir;
        glob_av2[i] += gdir*gdir;
        err_gdir = Error(glob_av[i],glob_av2[i],iblk);
        
        if (iblk == nblocks){   // saves final g(r) w/ errors
            dV = pow(r+bin_size/2,3)-pow(r-bin_size/2,3);
            Gave << std::setw(wd) << r
            << std::setw(wd) << glob_av[i]/(double)nblocks/norm/dV
            << std::setw(wd) << err_gdir/dV/norm << std::endl;
        }
    }
    Gave.close();
#endif
    return;
}

// Write final configuration
// generalised also for second-to-last configuation
void ConfFinal(void){
    std::vector<std::string> extension = {"old","final"};
    for (jndex=0; jndex<=1; ++jndex){
        std::ofstream WriteConf;
        std::cout << sys << " Print " + extension[jndex] + " configuration to file config." + extension[jndex] << std::endl;
        WriteConf.open("config." + extension[jndex]);
        
        for (ind=0; ind < npart; ++ind){
            if (jndex == 0)
                WriteConf << std::setw(wd) << xold[ind]/box << std::setw(wd) <<  yold[ind]/box << std::setw(wd) << zold[ind]/box << std::endl;
            if (jndex==1)
                WriteConf << std::setw(wd) << x[ind]/box << std::setw(wd) <<  y[ind]/box << std::setw(wd) << z[ind]/box << std::endl;
        }
        WriteConf.close();
    }
    
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
    std::ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
    WriteXYZ << npart << std::endl;
    WriteXYZ << "This is only a comment!" << std::endl;
  for (int i=0; i<npart; ++i){
      WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << std::endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************/
