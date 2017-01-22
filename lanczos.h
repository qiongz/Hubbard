#ifndef LANCZOS_H
#define LANCZOS_H
class lbasis{
public:
      vector<double> coeff;
      int size;
      double norm;
      explicit lbasis();
      lbasis(int n=0):size(n);
      lbasis(const lbasis &);
      ~lbasis();
      void init(int n):size(n);
      double normalize();
      lbasis &operator=(const lbasis &);
      const lbasis operator+(const lbasis &);
      const lbasis operator*(const double &);
      double &operator*(const lbasis &);
      const lbasis hoperation=(const vector<double> &,const vector<int> &, const vector<int> &);
};

extern "C" int dsyev_(char *, char *, int *, double *, int*, double *, double *, int *, int *);
void diag(double *h, double *e, int l);
void hoperation(double*,double*);

#endif
