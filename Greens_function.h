#ifndef GREENS_FUNCTION_H
#define GREENS_FUNCTION_H
#include"lanczos_hamiltonian.h"
#include"basis.h"

class Greens_func{
public:
        lhamil gs_config,hole_config,particle_config;
        basis gs_sector,hole_sector,particle_sector;
        vector<double> hole_phi_0,particle_phi_0;
        Greens_func();
        Greens_func(lhamil & _gs_config,basis & _gs_sector);
        ~Greens_func();
        void creation(long r,double coeff,int  alpha);
        void annihilation(long r,double coeff,int alpha);
        void spectral_function_ii_hole(long r_i, double eta,vector<double> &E, vector<double> &A,double &mu, int alpha);
        void spectral_function_ii_particle(long r_i, double eta,vector<double> &E, vector<double> &A,double &mu, int alpha);
        void spectral_function_ij_hole(long r_i,long r_j,double eta,vector<double> &E, vector<double> &A,double &mu,int alpha);
        void spectral_function_ij_particle(long r_i,long r_j,double eta,vector<double> &E, vector<double> &A,double &mu, int alpha);
        void spectral_function_kk_hole(double k, double eta,vector<double> &E,vector<double> &A,double &mu,int alpha);
        void spectral_function_kk_particle(double k, double eta,vector<double> &E,vector<double> &A,double &mu,int alpha);

};
#endif
