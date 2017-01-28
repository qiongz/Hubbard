#ifndef LANCZOS_H
#define LANCZOS_H
#include"basis.h"
class lbasis{
public:
      std::vector<double> coeff;
      int size;
      double norm;
      lbasis();
      lbasis(const int);
      lbasis(const lbasis &);
      ~lbasis();
      void init_zeros();
      void init_random();
      double normalize();
      lbasis & operator=(const lbasis &);
      const lbasis operator+(const lbasis &) const;
      const lbasis operator-(const lbasis &) const;
      const lbasis operator*(const double &)const ;
      const double operator*(const lbasis &)const;
      lbasis hoperation(const vector<double> &,const vector<int> &, const vector<int> &);
};

extern "C" int dsyev_(char *, char *, int *, double *, int*, double *, double *, int *, int *);
void diag(double *h, double *e, int l);

#endif
