#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"basis.h"
#include"matrix.h"
class hamil {
public:
    /** Size of the Hilbert space */
    long nHilbert;
    /** seed for the RNGs */
    unsigned seed;
    /** Hamiltonian matrix in CSR format */
    Mat H;
    /** Eigenvalues of the hamiltonian */
    std::vector<double> eigenvalues;
    /** Ground state wave function */
    std::vector<double> psi_0;
    /** First element of all wave functions */
    std::vector<double> psi_n0;

    hamil();
    ~hamil();
    /** Initialize the hamiltonian matrix elements.
     \param sector basis sector,
     \param t hopping strength,
     \param U onsite replusive interaction strength
     */
    hamil(basis & sector,double t, double U);
    /** Return the ground state energy of the system */
    double ground_state_energy();
    /** Diagonalize the full hamiltonian */
    void diag();

    complex<double> Greens_function(double,double);
    /** Print the hamiltonian matrix */
    void print_hamil();
    /** Print the eigenvalues of the system */
    void print_eigen();
};
#endif
