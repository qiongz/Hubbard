/************************************************
Basis generation library for 1D Hubbard model
   @author qiongzhu
   20/1/2017
   Email:qiongzhunow@gmail.com
*************************************************/
#ifndef BASIS_H
#define BASIS_H
#include<iostream>
#include<string.h>
#include<iomanip>
#include<bitset>
#include<cstdlib>
#include<map>
#include<vector>
#include<algorithm>
using namespace std;

class basis {
public:
    long nsite,nel_up,nel_down;  // No. of sites, spin-up/down electrons
    map<long,long> basis_up,basis_down; // basis set of spin-up/down electrons, I-J table

    long nbasis_up,nbasis_down;     // No. of basis for spin-up/down electrons
    vector<long> id_up,id_down;     // reversal table, J->I, Lin's Table is a 2D array
    explicit basis();
    basis(long,long,long);
    const basis & operator=(const basis &);
    ~basis();
    long hopping_up(long,long,long);
    long hopping_down(long,long,long);
    long potential(long,long,long);
    long onsite_up(long,long);
    long onsite_down(long,long);
    long factorial(long,long);
    void init();
    void init(long,long,long);
    void generate_up(long);
    void generate_down(long);
    long creation(long,long);
    long annihilation(long,long);
    void prlong();
    friend ostream & operator<<(ostream & os, const basis &);
};

#endif
