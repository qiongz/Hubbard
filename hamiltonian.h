#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"basis.h"
#include"matrix.h"
class hamil {
public:
    long nHilbert;
    unsigned seed;
    Mat H;
    std::vector<double> eigenvalues;
    std::vector<double> psi_0,psi_n0;

    hamil();
    ~hamil();
    hamil(basis &,double, double);
    double ground_state_energy();
    void diag();
    complex<double> Greens_function(double,double);
    void print_hamil();
};
#endif
