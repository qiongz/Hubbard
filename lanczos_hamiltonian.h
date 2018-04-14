#ifndef LANCZOS_HAMILTONIAN_H
#define LANCZOS_HAMILTONIAN_H
#include"matrix.h"
#include"basis.h"
void swap(Vec *a,Vec *b,Vec *c);

class lhamil {
public:
    unsigned seed;  //!< Seed for RNGs
    long nHilbert;  //!< Hilbert space size
    long lambda;    //!< Lanczos update steps
    double t,U,V;
    double E0;      //!< Ground state eigen energy
    basis sector;   //!< Basis
    Mat H;  //!< Hamiltonian matrix in CSR format
    //Mat O;  //!< Operator matrix in CSR format
    std::vector<double> norm; //!< Normalization coefficients vector in Lanczos update
    std::vector<double> overlap; //!< Overlap coefficients vector in Lanczos update
    std::vector<double> psir_0; //!< Ground state wave function in real-space
    std::vector<double> psi_0; //!<Ground state eigenvector in Krylov subspace
    std::vector<double> psi_n0; //!<First element of eigenvectors in Krylov subspace
    std::vector<double> eigenvalues; //!< Eigenvalues

    lhamil();  //!< Empty constructor
    /**
     \param _H hamiltonian matrix
     \param _nHilbert Hilbert space size
     \param _lambda Lanczos update steps
     \param _seed Seed for RNGs
    */
    lhamil(const Mat & _H,long _nHilbert,long _lambda,unsigned _seed);  //!< Constructor with hamiltonian matrix as input
    /**
     \param _sector Basis sector
     \param t Hopping strength
     \param U Onsite replusive interaction strength
     \param _lambda Lanczos update steps
     \param _seed Seed for RNGs
    */
    lhamil(basis & _sector,double V,double t, double U,long _lambda,unsigned _seed); //!< Constructor with basis sector as input
    ~lhamil(); //!< Destructor
    void init(basis &_sector,double V,double t, double U, long _lambda,unsigned _seed);
    const lhamil & operator=(const lhamil &);
    /** \param _sector Basis sector
        \param t Hopping strength
        \param U Onsite replusive interaction strength
    */
    void set_hamil(basis & _sector ,double V,double t ,double U);  //!< Initialize hamiltonian matrix
    void coeff_update(); //!< Lanczos update implemenation utilizing the Mat class
    void coeff_explicit_update(); //!< Lanczos update implemenation written in explicit arrays
    void coeff_update_wopt(vector<double> O_phi_0);
    void diag();  //!< Diagonalize the full Lanczos hamiltonian
    void diag(int l); //!< Diagonalize the Lanczos hamiltonain with first lxl elements

    void eigenstates_reconstruction(); //!< Transform |psi_0> to |psir_0>
    double ground_state_energy();    //!< Ground state energy
    double spectral_function(double omega,double eta); //!< Spectral function with spin, continued fraction version
    void print_hamil(int range); //!< print the full hamiltonian matrix
    void print_lhamil(int range);  //!< print the Lanczos hamiltonian matrix with first n x n elements
    void print_eigen(int range);  //!< print the first n eigenvalues
    void save_to_file(const char* filename);  //!< save object to file "filename"
    void read_from_file(const char*);        //!< load object from file "filename"
};
#endif
