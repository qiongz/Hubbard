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

class basis{
public:
    int nsite,nel_up,nel_down;          // No. of sites, spin-up/down electrons
    map<int,int> basis_up,basis_down;   // basis set of spin-up/down electrons, I-J table

    int nbasis_up,nbasis_down;          // No. of basis for spin-up/down electrons
    vector<int> id_up,id_down;          // reversal table, J->I, Lin's Table is a 2D array 
    explicit basis();
    basis(int,int,int);
    const basis & operator=(const basis &);
    ~basis();
    int hopping_up(int,int);
    int hopping_down(int,int);
    int potential(int,int,int);
    int factorial(int,int);
    void init();
    void generate_up(int);
    void generate_down(int);
    void print();
    friend ostream & operator<<(ostream & os, const basis &); 
};

#endif 
