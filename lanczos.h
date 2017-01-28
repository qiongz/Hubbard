#ifndef LANCZOS_H
#define LANCZOS_H
#include"basis.h"
class lbasis{
public:
      std::vector<double> coeff;
      int size;
      double norm;
      explicit lbasis();
      lbasis(const int n=0):size(n){};
      lbasis(const lbasis &);
      ~lbasis();
      void init(const int);
      double normalize();
      lbasis & operator=(const lbasis &);
      const lbasis operator+(const lbasis &);
      const lbasis operator*(const double &);
      double operator*(const lbasis &);
      lbasis hoperation(const std::vector<double> &,const std::vector<int> &, const std::vector<int> &);
};

extern "C" int dsyev_(char *, char *, int *, double *, int*, double *, double *, int *, int *);
void diag(double *h, double *e, int l);

#endif
