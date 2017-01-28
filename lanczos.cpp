#include"lanczos.h"
lbasis::lbasis(){}

lbasis::lbasis(const int _size):size(_size){
    init(size);
}

lbasis::lbasis(const lbasis &rhs){
    size=rhs.size;
    for(int i=0;i<size;i++)
        coeff.assign((rhs.coeff).begin(),(rhs.coeff).end());
}

lbasis::~lbasis(){}

lbasis & lbasis::operator=(const lbasis & rhs){
   size=rhs.size;
   norm=rhs.norm;
   coeff.assign((rhs.coeff).begin(),(rhs.coeff).end());
   return *this;
}

double lbasis::operator*(const lbasis &rhs){
    if(this->size!=rhs.size) return 0 ;
    double overlap=0;
    for(int i=0;i<this->size;i++)
        overlap+=coeff[i]*(rhs.coeff)[i];
    return overlap;
}

const lbasis lbasis::operator+(const lbasis & rhs){
    lbasis sum(*this);
    for(int i=0;i<size;i++)
       (sum.coeff)[i]=(sum.coeff)[i]+(rhs.coeff)[i];
    return sum;
}

const lbasis lbasis::operator*(const double &rhs){
    lbasis sum(*this);
    for(int i=0;i<size;i++)
        (sum.coeff)[i]*=rhs;
    return sum;
}

void lbasis::init(const int _n){
    int i=0;
    size=_n;
    coeff.reserve(size);
    norm=0;
    for(i=0;i<size;i++){
        coeff.push_back(rand()*1.0/RAND_MAX-0.5);
        norm+=coeff.back()*coeff.back();
    }
    norm=sqrt(norm);
    for(i=0;i<size;i++)
        coeff[i]/=norm;
}

lbasis lbasis::hoperation(vector<double> &hamil,vector<int> & row, vector<int> &col){
    lbasis wf_2(this);
    wf_2.coeff.assign(wf_2.size(),0);
    for(int i=0;i<size;i++)
        wf2.coeff[row[i]]+=hamil[i]*(*this).coeff[col[i]];

    return wf_2;
}


void lbasis::normalize() {
    int i;
    for(i=0; i<size; i++)
        norm+=coeff[i]*coeff[i];
    norm=sqrt(norm);
    for(i=0; i<size; i++)
        coeff[i]/=norm;
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



