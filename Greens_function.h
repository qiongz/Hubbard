#ifndef GREENS_FUNCTION_H
#define GREENS_FUNCTION_H
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include"basis.h"

class Greens_func{
public:
        hamil gs_hconfig,hole_hconfig,particle_hconfig;
        lhamil gs_config,hole_config,particle_config;
        basis gs_sector,hole_sector,particle_sector;
        vector<double> hole_phi_0,particle_phi_0;
        Greens_func();
        Greens_func(lhamil & _gs_config);
        Greens_func(hamil & _gs_config,basis & _gs_sector);
        ~Greens_func();
        void creation_u(long r,double coeff);
        void annihilation_u(long r,double coeff);
        void creation_u_full_hamil(long r,double coeff);
        void annihilation_u_full_hamil(long r,double coeff);
        void spectral_function_ii_uu_hole(long r_i, double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_ii_uu_particle(long r_i, double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_ij_uu_hole(long r_i,long r_j,double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_ij_uu_particle(long r_i,long r_j,double eta,vector<double> &E, vector<double> &A,double &mu);
        void spectral_function_kk_uu_hole(double k, double eta,vector<double> &E,vector<double> &A,double &mu);
        void spectral_function_kk_uu_particle(double k, double eta,vector<double> &E,vector<double> &A,double &mu);

};
#endif
