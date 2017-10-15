#ifndef LANCZOS_HAMILTONIAN_H
#define LANCZOS_HAMILTONIAN_H
#include"matrix.h"
#include"basis.h"
void swap(Vec *a,Vec *b,Vec *c);

class lhamil {
public:
    unsigned seed;
    long nHilbert,lambda;
    std::vector<double> norm,overlap;
    std::vector<double> psi_0,psi_n0,eigenvalues;
    Vec Psi_0;
    Mat H;
    std::vector< Vec > Vec_list;

    lhamil();
    lhamil(const Mat &,long,long,unsigned);
    lhamil(basis &,double, double,long,unsigned);
    ~lhamil();

    void keep_basis_update();
    void coeff_update();
    void coeff_explicit_update();
    void diag();
    void eigenstates_reconstruction();

    void print_hamil(int);
    double ground_state_energy();
    complex<double> Greens_function(double,double);

    void save_to_file(const char*);
    void read_from_file(const char*);
};
#endif
