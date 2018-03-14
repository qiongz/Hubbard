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
        Greens_func(lhamil & _gs_config);
        ~Greens_func();
        void creation_u(long r,double coeff);
        void annihilation_u(long r,double coeff);
        void spectral_function_ii_uu(int r_i, double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_ij_uu(int r_i,long r_j,double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_kk_uu(double k, double eta,vector<double> &E,vector<double> &A,double &mu);
};
#endif
