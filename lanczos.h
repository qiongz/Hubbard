#ifndef LANCZOS_H
#define LANCZOS_H
extern "C" int dsyev_(char *, char *, int *, double *, int*, double *, double *, int *, int *);
void diag(double *h, double *e, int l);

#endif