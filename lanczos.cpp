#include"lanczos.h"
lbasis::lbasis(){}

lbasis::lbasis(const int _size){
    size=_size;
}

lbasis::lbasis(const lbasis &rhs){
    size=rhs.size;
    for(int i=0;i<size;i++)
        coeff.assign((rhs.coeff).begin(),(rhs.coeff).end());
}

lbasis::~lbasis(){}

lbasis & lbasis::operator=(const lbasis & rhs){
   size=rhs.size;
   coeff.assign((rhs.coeff).begin(),(rhs.coeff).end());
   return *this;
}

const double lbasis::operator*(const lbasis &rhs)const{
    if(size!=rhs.size) return 0 ;
    double overlap=0;
    for(int i=0;i<size;i++)
        overlap+=coeff[i]*(rhs.coeff)[i];
    return overlap;
}

const lbasis lbasis::operator+(const lbasis & rhs)const {
    lbasis sum(size);
    sum.init_zeros();
    for(int i=0;i<size;i++)
       (sum.coeff)[i]=coeff[i]+(rhs.coeff)[i];
    return sum;
}

const lbasis lbasis::operator-(const lbasis & rhs)const {
    lbasis minus(size);
    minus.init_zeros();
    for(int i=0;i<size;i++)
       (minus.coeff)[i]=coeff[i]-(rhs.coeff)[i];
    return minus;
}

const lbasis lbasis::operator*(const double &rhs)const {
    lbasis prod(size);
    for(int i=0;i<size;i++)
        (prod.coeff).push_back(coeff[i]*rhs);
    return prod;
}

ostream & operator<<(ostream & os, const lbasis & _lb){
  os<<"[ ";
  for(int i=0;i<_lb.size-1;i++)
    os<<setprecision(6)<<setw(10)<<(_lb.coeff)[i]<<", ";
  os<<setprecision(6)<<setw(10)<<(_lb.coeff)[_lb.size-1]<<" ]";
}

void lbasis::init_zeros(){
    coeff.assign(size,0);
    norm=0; 
}


void lbasis::init_random(){
    coeff.reserve(size);
    norm=0;
    for(int i=0;i<size;i++){
        coeff.push_back(rand()*1.0/RAND_MAX-0.5);
        //coeff.push_back(rand()*1.0/RAND_MAX);
        norm+=coeff.back()*coeff.back();
    }
    norm=sqrt(norm);
    for(int i=0;i<size;i++)
        coeff[i]/=norm;
}

lbasis lbasis::hoperation(const vector<double> &hamil,const vector<int> & row, const vector<int> &col){
    lbasis wf_2(size);
    wf_2.init_zeros();
    for(int i=0;i<size;i++)
        wf_2.coeff[row[i]]+=hamil[i]*coeff[col[i]];
    return wf_2;
}


double lbasis::normalize() {
    int i;
    norm=0;
    for(i=0; i<size; i++)
        norm+=coeff[i]*coeff[i];
    norm=sqrt(norm);
    for(i=0; i<size; i++)
        coeff[i]/=norm;
    return norm;
}

void diag(double *h, double *e, int l){
    char jobz,uplo;
    int info;
    jobz = 'V';
    uplo = 'U';
    int lda=l;
    int lwork = 3*l-1;
    double *work=new double[lwork];
    dsyev_(&jobz, &uplo, &l, h, &lda, e, work, &lwork, &info);
    delete [] work;
}



