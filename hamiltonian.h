#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include"lanczos.h"
#include"basis.h"

double normalize(double *,int);
double overlap(double*,double*,int);
void hoperation(double*,double*,std::vector<double>&,std::vector<int> &,std::vector<int> &,int);
void diag_hamil(basis* ,double ,double ,double, double *);

#endif
