#include"lanczos_hamiltonian.h"
inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,long _nHilbert,long _lambda, unsigned _seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(basis &sector,double t, double U,long _lambda,unsigned _seed) {
    long nsite,nbasis_up,nbasis_down;
    lambda=_lambda;
    seed=_seed; 
    nsite=sector.nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    // CSR format of matrix storage
    std::vector<long> inner_indices, outer_starts;
    std::vector<double> matrix_elements;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    /*
    std::map<int,double> hamil_nonzero;
    std::map<int,double>::iterator it;
    std::vector<int> row_index,col_index;
    int matrix_index;
    */
    long n,i,j,k,l;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            for(n=0; n<nsite-1; n++) {
                k=sector.hopping_up(i,n);
                l=sector.hopping_down(j,n);
                //matrix_index=(i*nbasis_down+j)*nHilbert+k*nbasis_down+l;
                if(k!=i) {
                    //matrix_index=(i*nbasis_down+j)*nHilbert+k*nbasis_down+j;
                    /*
                    it=hamil_nonzero.find(matrix_index);
                    
                    if(it==hamil_nonzero.end())
                        hamil_nonzero.insert(pair<int,double>(matrix_index,-t));
                    else
                        it->second=it->second-t;
                    */
                    row++;
                    inner_indices.push_back(k*nbasis_down+j);
                    matrix_elements.push_back(-t);
                }

                if(l!=j) {
                    /*matrix_index=(i*nbasis_down+j)*nHilbert+i*nbasis_down+l;
                    it=hamil_nonzero.find(matrix_index);
                    if(it==hamil_nonzero.end())
                        hamil_nonzero.insert(pair<int,double>(matrix_index,-t));
                    else
                        it->second=it->second-t;
                    */
                    row++;
                    inner_indices.push_back(i*nbasis_down+l);
                    matrix_elements.push_back(-t);
                }
                if(sector.potential(i,j,n)){
                   row++;
                   inner_indices.push_back(i*nbasis_down+j);
                   matrix_elements.push_back(U);
                }
            }
            if(sector.potential(i,j,nsite-1)){
                   row++;
                   inner_indices.push_back(i*nbasis_down+j);
                   matrix_elements.push_back(U);
            }
            outer_starts.push_back(row);
            /*
            for(n=0; n<nsite; n++)
                if(sector.potential(i,j,n)) {
                    matrix_index=(i*nbasis_down+j)*nHilbert+i*nbasis_down+j;
                    full_hamil[matrix_index]+=U;

                    it=hamil_nonzero.find(matrix_index);
                    if(it==hamil_nonzero.end())
                        hamil_nonzero.insert(pair<int,double>(matrix_index,U));
                    else
                        it->second=it->second+U;
                }
            */
        }
    }

    H.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}

lhamil::~lhamil() {
}

void lhamil::keep_basis_update() {
    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    Vec phi_0,phi_1,phi_2;
    Vec_list.reserve(lambda);

    phi_0.init_random(nHilbert,seed);
    Vec_list.push_back(phi_0);

    phi_1 = H * phi_0;
    overlap[0] = phi_0 * phi_1;
    phi_1 -= phi_0 * overlap[0];

    norm[0]= 1;
    norm[1]= phi_1.normalize();
    /* reorthogonalize for degenerate eigenstates
       prevent loos of orthogonality is essential for Lanczos */
   // double q = phi_0 * phi_1;
   // phi_1 = (phi_1 - phi_0 * q) / (1 - q * q);

    Vec_list.push_back(phi_1);
    for(int i = 1; i < lambda; i++) {
        phi_2 = H * phi_1;
        overlap[i] = phi_1*phi_2;
        phi_2 -= phi_1 * overlap[i] + phi_0 * norm[i];
        norm[i+1] = phi_2.normalize();
        phi_0 = phi_1;
        phi_1 = phi_2;

        Vec_list.push_back(phi_2);
        /*
        for(int j = 0; j <= i; j++) {
            q = Vec_list[j] * Vec_list[i + 1];
            Vec_list[i + 1] = (Vec_list[i + 1] - Vec_list[j] * q) / (1 - q * q);
        }
        */
    }
}

void lhamil::coeff_update() {
    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    Vec phi_0,phi_1,phi_2;
    phi_0.init_random(nHilbert,seed);
    norm[0]=1;
    phi_1 = H*phi_0;
    overlap[0] = phi_0 * phi_1;
    phi_1 -= phi_0 * overlap[0];
    norm[1] =phi_1.normalize();
    for(int i = 1; i < lambda; i++) {
        phi_2= H * phi_1;
        overlap[i] = phi_1 * phi_2;
        #pragma ivdep
        phi_2 -= phi_1 * overlap[i] + phi_0 * norm[i];
        norm[i+1] = phi_2.normalize();
        swap(&phi_0,&phi_1,&phi_2);
    }
}

void lhamil::coeff_explicit_update()
{
    int i,j,idx;
    double norm_factor,overlap_factor;
    double *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new double[nHilbert];
    phi_1=new double[nHilbert];
    phi_2=new double[nHilbert];
    phi_t=new double[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    //phi_0.init_random(nHilbert,seed);
    std::mt19937 rng(seed);
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]=rng()*1.0/rng.max()-0.5;

    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_0[i]*phi_0[i];
    norm_factor=sqrt(norm_factor);

    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_0[i]/=norm_factor;
    norm[0]=1;

    // phi_1=H*phi_0;
    memset(phi_1,0,sizeof(double)*nHilbert);
    #pragma omp parallel for schedule(static)
    for(i=0; i<H.outer_starts.size()-1; i++) {
        #pragma ivdep
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            phi_1[i]+=H.value[idx]*phi_0[H.inner_indices[idx]];
    }

    //overlap[0] = phi_0 * phi_1;
    overlap_factor=0;
    #pragma omp parallel for reduction(+:overlap_factor)
    for(i=0; i<nHilbert; i++)
        overlap_factor+=phi_1[i]*phi_0[i];
    overlap[0]=overlap_factor;

    //phi_1 -= phi_0 * overlap[0];
    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]-=phi_0[i]*overlap_factor;

    //norm[1] = phi_1.normalize();
    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_1[i]*phi_1[i];
    norm_factor=sqrt(norm_factor);

    #pragma omp parallel for schedule(static)
    for(i=0; i<nHilbert; i++)
        phi_1[i]/=norm_factor;
    norm[1]=norm_factor;

    for(j= 1; j < lambda; j++) {
        memset(phi_2,0,sizeof(double)*nHilbert);
        #pragma omp parallel for schedule(static)
        for(i=0; i<H.outer_starts.size()-1; i++) {
            #pragma ivdep
            for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
                phi_2[i]+=H.value[idx]*phi_1[H.inner_indices[idx]];
        }

        //overlap[j] = phi_1 * phi_2;
        overlap_factor=0;
        #pragma omp parallel for reduction(+:overlap_factor)
        for(i=0; i<nHilbert; i++)
            overlap_factor+=phi_1[i]*phi_2[i];
        overlap[j]=overlap_factor;

        //phi_2 -= phi_1 * overlap[j] + phi_0 * norm[j];
        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]=phi_0[i]*norm[j];

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_t[i]+=phi_1[i]*overlap_factor;

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]-=phi_t[i];

        //norm[j+1] = phi_2.normalize();
        norm_factor=0;
        #pragma omp parallel for reduction(+:norm_factor)
        for(i=0; i<nHilbert; i++)
            norm_factor+=phi_2[i]*phi_2[i];
        norm_factor=sqrt(norm_factor);

        #pragma omp parallel for schedule(static)
        for(i=0; i<nHilbert; i++)
            phi_2[i]/=norm_factor;
        norm[j+1]=norm_factor;

        phi_s=phi_0;
        phi_0=phi_1;
        phi_1=phi_2;
        phi_2=phi_s;
    }
    delete phi_0,phi_1,phi_2,phi_t;
}

void lhamil::diag() {
    if(norm.size()==0) coeff_update();
    int l=lambda;
    double *e = new double[l];
    double *h = new double [l * l];
    memset(h, 0, sizeof(double)*l*l);
    for(int i = 0; i < l-1; i++) {
        h[i *l + i + 1] = norm[i+1];
        h[(i + 1)*l + i] = norm[i+1];
        h[i * l + i] = overlap[i];
    }
    h[(l - 1)*l + l - 1] = overlap[l - 1];
    diag_dsyev(h,e,l);

    psi_0.assign(l,0);
    psi_n0.assign(l,0);
    eigenvalues.assign(l,0);

    for(int i=0; i<l; i++) {
        psi_0[i]=h[i];
        psi_n0[i]=h[i*l];
        eigenvalues[i]=e[i];
    }
    delete h,e;
}

// given the seed of phi_0, generate the Lanczos basis again, and transform the eigenstate to the desired eigenstate representation
void lhamil::eigenstates_reconstruction() {
    int l,n;
    double Norm=0;
    
    // if the Vec_list is prestored
    if(Vec_list.size()!=0){
        Psi_0.assign(nHilbert,0);
        for(l=0; l<lambda; l++)
            for(n=0; n<nHilbert; n++)
                Psi_0.value[n]+=psi_0[l]*Vec_list[l].value[n];
        for(n=0; n<nHilbert; n++)
            Norm+=Psi_0.value[n]*Psi_0.value[n];
        for(n=0; n<nHilbert; n++)
            Psi_0.value[n]/=sqrt(Norm);
        return;
    }
   
    Vec phi_0,phi_1,phi_2;
    // repeat the iteration, but with normalization and overlap values 
    // at hand
    phi_0.init_random(nHilbert,seed);
    phi_1 = H * phi_0;
    phi_1 -= phi_0*overlap[0];
    phi_1 /= norm[1];
    Psi_0.assign(nHilbert,0);
    for(n=0; n<nHilbert; n++)
        Psi_0.value[n]+=psi_0[0]*phi_0.value[n]+psi_0[1]*phi_1.value[n];

    for(l=2; l<lambda; l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++)
            Psi_0.value[n]+=psi_0[l]*phi_2.value[n];
        swap(&phi_0,&phi_1,&phi_2);
    }
    for(n=0; n<nHilbert; n++)
        Norm+=Psi_0.value[n]*Psi_0.value[n];
    for(n=0; n<nHilbert; n++)
        Psi_0.value[n]/=sqrt(Norm);
}

void lhamil::print_hamil( int range) {
    H.print();
}

double lhamil::ground_state_energy() {
    if(Psi_0.value.size()!=0) return Psi_0*(H*Psi_0);
       else   
         return 0; 
}


/*
complex<double> lhamil::Greens_function(double E_real,double epsilon) {
    complex<double> E=complex<double>(E_real,epsilon);
    complex<double> t=0.5*(E-overlap.back()-sqrt(pow(E-overlap.back(),2)-4*norm[lambda-1]*norm[lambda-1]));
    for(int n=lambda-2; n>0; n--) {
        t=E-overlap[n]-norm[lambda-1]*norm[n+1]/t;
    }
    return 1.0/(E-overlap[0]-norm[1]*norm[1]/t);
}
*/

complex<double> lhamil::Greens_function(double E_real,double epsilon) {
    complex<double> E=complex<double>(E_real,epsilon);
    complex<double> G=0;
    for(int i=0; i<lambda; i++)
        G+=psi_n0[i]*psi_n0[i]/(E-eigenvalues[i]);
    return G;
}


void lhamil::save_to_file(const char* filename){
    if(eigenvalues.size()==0) return;
    ofstream odf;
    odf.open(filename,ofstream::out);
    odf<<"nHilbert:="<<setw(10)<<" "<<nHilbert<<endl;
    odf<<"lambda:="<<setw(10)<<" "<<lambda<<endl;
    odf<<"seed:="<<setw(10)<<" "<<seed<<endl;
    odf<<"eigenvalues:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0;i<lambda;i++)
       odf<<eigenvalues[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"norm:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0;i<=lambda;i++)
       odf<<norm[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"overlap:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0;i<lambda;i++)
       odf<<overlap[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0;i<lambda;i++)
       odf<<psi_0[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_n0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0;i<lambda;i++)
       odf<<psi_n0[i]<<" , ";
    odf<<" ] "<<endl;
    odf.close();
}

void lhamil::read_from_file(const char* filename){
    ifstream idf;
    idf.open(filename,ifstream::in);
    string buffer;
    idf>>buffer;
    idf>>nHilbert;
    idf>>buffer;
    idf>>lambda;
    idf>>buffer;
    idf>>seed;
    double value;
    idf>>buffer;
    idf>>buffer;
    for(int i=0;i<lambda;i++){
        idf>>value;
        idf>>buffer;
        eigenvalues.push_back(value);
        }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0;i<=lambda;i++){
        idf>>value;
        idf>>buffer;
        norm.push_back(value);
        }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0;i<lambda;i++){
        idf>>value;
        idf>>buffer;
        overlap.push_back(value);
        }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0;i<lambda;i++){
        idf>>value;
        idf>>buffer;
        psi_0.push_back(value);
        }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0;i<lambda;i++){
        idf>>value;
        idf>>buffer;
        psi_n0.push_back(value);
       }   
    idf.close();
}
