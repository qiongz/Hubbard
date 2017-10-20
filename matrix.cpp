#include"matrix.h"
Vec::Vec() {}

Vec::Vec(long _size) {
    size = _size;
    value.reserve(size);
}

Vec::Vec(long _size,const double _init) {
    size = _size;
    value.assign(size, _init);
}

Vec::Vec(const Vec &rhs) {
    size = rhs.size;
    for(int i = 0; i < size; i++)
        value.assign((rhs.value).begin(), (rhs.value).end());
}

Vec::~Vec() {
    if(size!=0) {
        value.clear();
        size=0;
    }
}

void Vec::assign(long  _size,const double _init) {
    size=_size;
    value.assign(size,_init);
}

void Vec::init_random(unsigned seed) {
    double norm=0;
    std::mt19937 rng(seed);
    #pragma omp parallel for reduction(+:norm)
    for(int i = 0; i < size; i++) {
        value[i] = rng() * 1.0 / rng.max()-0.5 ;
        norm += value[i] * value[i];
    }
    norm = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= norm;
}

void Vec::init_random(long _size,unsigned seed) {
    size=_size;
    value.resize(size);
    std::mt19937 rng(seed);
    double norm=0;
    #pragma omp parallel for reduction(+:norm)
    for(int i = 0; i < size; i++) {
        value[i] = rng() * 1.0 / rng.max() - 0.5;
        norm += value[i] * value[i];
    }
    norm = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= norm;
}

void Vec::clear() {
    if(size!=0) {
        value.clear();
        size=0;
    }
}

double Vec::normalize() {
    int i;
    double norm = 0;
    #pragma omp parallel for reduction(+:norm)
    for(i = 0; i < size; i++)
        norm += value[i] * value[i];
    double normsq = sqrt(norm);
    #pragma omp parallel for schedule(static)
    for(i = 0; i < size; i++)
        value[i] /= normsq;
    return normsq;
}

Vec & Vec::operator=(const Vec & rhs) {
    if(this==&rhs) return *this;
    size = rhs.size;
    value.assign((rhs.value).begin(), (rhs.value).end());
    return *this;
}

Vec & Vec::operator-=(const Vec & rhs) {
    if(this==&rhs)
        return *this;
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] -= rhs.value[i];
    return *this;
}

Vec & Vec::operator+=(const Vec & rhs) {
    if(this==&rhs)
        return *this;
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] += rhs.value[i];
    return *this;
}

Vec & Vec::operator*=(const double & rhs) {
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] *=rhs;
    return *this;
}

Vec & Vec::operator/=(const double & rhs) {
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        value[i] /= rhs;
    return *this;
}

Vec Vec::operator+(const Vec & rhs) {
    Vec rt(rhs.size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] + (rhs.value)[i];
    return rt;
}

Vec Vec::operator-(const Vec & rhs) {
    Vec rt(rhs.size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] - (rhs.value)[i];
    return rt;
}

Vec Vec::operator/(const double &rhs) {
    Vec rt(size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i] / rhs;
    return rt;
}

Vec Vec::operator*(const double &rhs) {
    Vec rt(size);
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < size; i++)
        rt.value[i] = value[i]* rhs;
    return rt;
}

double Vec::operator*(const Vec &rhs) {
    if(size != rhs.size) return 0 ;
    double overlap = 0;
    #pragma omp parallel for reduction(+:overlap)
    for(int i = 0; i < size; i++)
        overlap += value[i] * (rhs.value)[i];
    return overlap;
}

ostream & operator<<(ostream & os, const Vec & _b) {
    os << "[ ";
    for(int i = 0; i < _b.size - 1; i++)
        os << setprecision(6) << setw(10) << (_b.value)[i] << ", ";
    os << setprecision(6) << setw(10) << (_b.value)[_b.size - 1] << " ]";
    return os;
}

Mat::Mat() {}

Mat::Mat(const Mat &rhs) {
    inner_indices.assign(rhs.inner_indices.begin(),rhs.inner_indices.end());
    value.assign(rhs.value.begin(),rhs.value.end());
    outer_starts.assign(rhs.outer_starts.begin(),rhs.outer_starts.end());
}

Mat::~Mat() {
    if(value.size()!=0) {
        value.clear();
        inner_indices.clear();
        outer_starts.clear();
    }
}

void Mat::clear() {
    if(value.size()!=0) {
        value.clear();
        inner_indices.clear();
        outer_starts.clear();
    }
}

Mat & Mat::operator=(const Mat & rhs) {
    if(this==&rhs) return *this;
    value.assign((rhs.value).begin(), (rhs.value).end());
    inner_indices.assign(rhs.inner_indices.begin(),rhs.inner_indices.end());
    outer_starts.assign(rhs.outer_starts.begin(),rhs.outer_starts.end());
    return *this;
}

Vec Mat::operator*(const Vec &rhs)const {
    Vec phi(rhs.size,0);
    if(rhs.size!=outer_starts.size()-1) return phi;
    #pragma omp parallel for schedule(guided,4)
    for(int i=0; i<outer_starts.size()-1; i++) {
#pragma ivdep
        for(int idx=outer_starts[i]; idx<outer_starts[i+1]; idx++)
            phi.value[i]+=value[idx]*rhs.value[inner_indices[idx]];
    }
    return phi;
}

vector<double> Mat::operator*(const vector<double> &rhs)const {
    vector<double> phi;
    phi.assign(rhs.size(),0);
    if(rhs.size()!=outer_starts.size()-1) return phi;
    #pragma omp parallel for schedule(guided,4)
    for(int i=0; i<outer_starts.size()-1; i++) {
#pragma ivdep
        for(int idx=outer_starts[i]; idx<outer_starts[i+1]; idx++)
            phi[i]+=value[idx]*rhs[inner_indices[idx]];
    }
    return phi;
}

void Mat::init(const vector<long> & _outer,const vector<long> & _inner, const vector<double> &_value) {
    outer_starts.assign(_outer.begin(),_outer.end());
    inner_indices.assign(_inner.begin(),_inner.end());
    value.assign(_value.begin(),_value.end());
}

void Mat::print() {
    std::cout<<"value:=         [";
    for(int i=0; i<value.size(); i++) std::cout<<setw(4)<<setprecision(2)<<value[i]<<" ";
    std::cout<<" ]"<<std::endl;
    std::cout<<"inner_indices:= [";
    for(int i=0; i<inner_indices.size(); i++) std::cout<<setw(4)<<setprecision(2)<<inner_indices[i]<<" ";
    std::cout<<" ]"<<std::endl;
    std::cout<<"outer_starts:=  [";
    for(int i=0; i<outer_starts.size(); i++) std::cout<<setw(4)<<setprecision(2)<<outer_starts[i]<<" ";
    std::cout<<" ]"<<std::endl;
}

void diag_dsyev(double *h, double *e, int l) {
    char jobz, uplo;
    int info;
    jobz = 'V';
    uplo = 'U';
    int lda = l;
    int lwork = 3 * l - 1;
    double *work = new double[lwork];
    dsyev_(&jobz, &uplo, &l, h, &lda, e, work, &lwork, &info);
    delete [] work;
}
