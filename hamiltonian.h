#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"basis.h"
#include"matrix.h"
class hamil {
public:
    /** Size of the Hilbert space */
    long nHilbert;
    /** seed for the RNGs */
    double t,U,V;
    /** Hamiltonian matrix in CSR format */
    double *hamiltonian;
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
    hamil(basis & sector,double _V,double _t, double _U);
    void init(basis & sector,double _V,double _t, double _U);
    const hamil & operator=(const hamil &);
    /** Return the ground state energy of the system */
    double ground_state_energy();
    /** Diagonalize the full hamiltonian */
    void diag();

    double spectral_function(vector<double> & O_phi_0,double omega,double _E0,double eta, int annil); //!< Spectral moments with spin
    /** Print the hamiltonian matrix */
    void print_hamil(int range);
    /** Print the eigenvalues of the system */
    void print_eigen(int range);
};
#endif
