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
    // psir_0 size nHilbert, psi_0,psi_n0,eigenvalues size nHilbert_ph
    std::vector<double> psir_0,psi_0,psi_n0,eigenvalues;
    basis sector;
    Mat H,O;

    lhamil();
    lhamil(const Mat &,long,long,unsigned);
    lhamil(basis &,double, double,long,unsigned);
    ~lhamil();

    void set_hamil(basis &,double,double);
    void psir0_creation_el_up(vector<double> &,long);
    void psir0_creation_el_down(vector<double> &,long);
    void psir0_annihilation_el_up(vector<double> &,long);
    void psir0_annihilation_el_down(vector<double> &,long);
    // set up operator O matrix
    void set_onsite_optc(int r,int alpha,int annil);
    
    /* Lanczos update */
    void coeff_update();
    void coeff_explicit_update();
    void coeff_update_wopt(vector<double>);

    // diagonalization the Lanczos hamiltonian
    void diag();
    void diag(int);

    // transform |psi_0> to |psir_0>
    void eigenstates_reconstruction();

    /* calculate physical quantities */
    double ground_state_energy();
    // spectral moments, diagonal part of Green's function
    double spectral_function(double omega,double eta,int annil);
    complex<double> Greens_function(double omega,double eta,int annil);
    // real-space Green's function with spin
    // G_{ij}^{\alpha\beta}(\omega) 
    complex<double> Greens_function_ij_ab(int i,int j,int alpha,int beta,double E,double eta);
    // k-space Green's function with spin
    // G_k^{\alpha\beta}(\omega)
    complex<double> Greens_function_k(int k,int alpha, int beta,double E,double eta);
    
    void print_hamil();
    void print_eigen(int);
    void save_to_file(const char*);
    void read_from_file(const char*);
};
#endif
