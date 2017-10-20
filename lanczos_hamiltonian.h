#ifndef LANCZOS_HAMILTONIAN_H
#define LANCZOS_HAMILTONIAN_H
#include"matrix.h"
#include"basis.h"
void swap(Vec *a,Vec *b,Vec *c);

class lhamil {
public:
    unsigned seed;
    long nHilbert,lambda;
    double E0;
    std::vector<double> norm,overlap;
    std::vector<double> psir_0,psi_0,psi_n0,eigenvalues;
    basis sector;
    Mat H,O;

    lhamil();
    lhamil(const Mat &,long,long,unsigned);
    lhamil(basis &,double, double,long,unsigned);
    ~lhamil();

    void coeff_update();
    void coeff_explicit_update();
    // initialize operator O and run the Lanczos update again
    void coeff_update_wopt();
    void set_onsite_optc(int r,int alpha,bool annil);
    void diag();
    void diag(int);
    void eigenstates_reconstruction();

    void print_hamil(int);
    void print_eigen(int);
    double ground_state_energy();
    // spectral moments, diagonal part of Green's function
    double spectral_function(double omega,double eta);
    // real-space Green's function with spin
    // G_{ij}^{\alpha\beta}(\omega) 
    complex<double> Greens_function_ij_ab(int i,int j,int alpha,int beta,double E,double eta);
    // k-space Green's function with spin
    // G_k^{\alpha\beta}(\omega)
    complex<double> Greens_function_k(int k,int alpha, int beta,double E,double eta);
    

    void save_to_file(const char*);
    void read_from_file(const char*);
};
#endif
