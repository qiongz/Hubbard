#include<iostream>
#include<bitset>
#include<iomanip>
#include<string.h>
#include<cstdlib>
using namespace std;

class basis{
private:
    int n_site,n_elu,n_eld;  // No. of sites, spin-up/down electrons
    int nb_up,nb_down;        // No. of basis for spin-up/down electrons
    map<int,int> elu,eld;     // basis set of spin-up/down electrons

public:
    explicit basis();
    basis(int,int,int);
    const basis & operator=(const basis &);
    ~basis();
    basis hopping(int,int);
    basis potential(int);
    void generate_up(int);
    void generate_down(int);
    void factorial(int,int);
    void print();
};
    
