#include"lanczos_hamiltonian.h"
inline void swap(Vec *a,Vec *b,Vec *c) {
    *a=*b;
    *b=*c;
}

lhamil::lhamil() {}

lhamil::lhamil(const Mat &_H,long _nHilbert,long _lambda, unsigned _seed):H(_H),nHilbert(_nHilbert),lambda(_lambda),seed(_seed) {}

lhamil::lhamil(basis &_sector,double t, double U,long _lambda,unsigned _seed) {
    sector=_sector;
    lambda=_lambda;
    seed=_seed;
    set_hamil(sector,t,U);
}

lhamil::~lhamil() {
}

void lhamil::set_hamil(basis & _sector,double t,double U)
{
    long nsite,nbasis_up,nbasis_down,signu,signd;
    sector=_sector;
    nsite=sector.nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    std::vector<long> inner_indices, outer_starts;
    std::vector<double> matrix_elements;
    std::map<long,double>::iterator it;
    inner_indices.reserve(nHilbert*nsite);
    matrix_elements.reserve(nHilbert*nsite);
    outer_starts.reserve(nHilbert+1);
    long n,m,i,j,k,l,s,nsignu,nsignd;
    long row=0;
    outer_starts.push_back(0);
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            std::map<long,double> col_indices;
            for(n=0; n<nsite; n++) {
                if(n+1>=nsite) {
                    signu=pow(-1,sector.nel_up-1);
                    signd=pow(-1,sector.nel_down-1);
                }
                else {
                    signu=1;
                    signd=1;
                }
                m=n+1;
                k=sector.hopping_up(i,n,m);
                if(k!=i) {
                    it=col_indices.find(k*nbasis_down+j);
                    if(it==col_indices.end())
                        col_indices.insert(std::pair<long,double>(k*nbasis_down+j,-t*signu));
                    else
                        it->second+=-t*signu;

                /*
                cout<<"----------------------"<<endl;
                cout<<"(n,m)=("<<n<<","<<m<<"): "<<endl; 
                cout<<"i:="<<bitset<4>(sector.id_up[i])<<endl;
                cout<<"k:="<<bitset<4>(sector.id_up[k])<<endl;
                cout<<"matrix element:="<<it->second<<endl;
                 */

                }
                l=sector.hopping_down(j,n,m);
                if(l!=j) {
                    it=col_indices.find(i*nbasis_down+l);
                    if(it==col_indices.end())
                        col_indices.insert(std::pair<long,double>(i*nbasis_down+l,-t*signd));
                    else
                        it->second+=-t*signd;
                }
                if(sector.potential(i,j,n)) {
                    it=col_indices.find(i*nbasis_down+j);
                    if(it==col_indices.end())
                        col_indices.insert(std::pair<long,double>(i*nbasis_down+j,U));
                    else
                        it->second+=U;
                }
            }
            for(it=col_indices.begin(); it!=col_indices.end(); it++) {
                inner_indices.push_back(it->first);
                matrix_elements.push_back(it->second);
                
                /*
                cout<<"----------------------"<<endl;
                cout<<"(i,j)=("<<i<<","<<j<<"): "<<endl; 
                cout<<bitset<5>(sector.id_up[i])<<" "<<bitset<5>(sector.id_down[j])<<endl;
                cout<<bitset<5>(sector.id_up[it->first/nbasis_down])<<" "<<bitset<5>(sector.id_down[it->first%nbasis_down])<<endl;
                cout<<"matrix element:="<<it->second<<endl;
                */
                
            }

            row+=col_indices.size();
            outer_starts.push_back(row);
            col_indices.clear();
        }
    }
    H.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}

void lhamil::set_onsite_optc(int r,int alpha,int annil)
{
    vector<long> inner_indices, outer_starts;
    vector<double> matrix_elements;
    inner_indices.reserve(nHilbert*sector.nsite);
    matrix_elements.reserve(nHilbert*sector.nsite);
    outer_starts.reserve(nHilbert+1);
    long n,i,j,k,l;
    long row=0;
    double potential_spin_up,potential_spin_down;
    outer_starts.push_back(0);
    for(i=0; i<sector.nbasis_up; i++)
        for(j=0; j<sector.nbasis_down; j++) {
            potential_spin_up=0;
            potential_spin_down=0;
             
            if(annil==1) {
                if(alpha==0 && sector.onsite_up(i,r)==1)
                    potential_spin_up++;

                if(alpha==1 && sector.onsite_down(j,r)==1)
                    potential_spin_down++;
            }
            else {
                if(alpha==0 && sector.onsite_up(i,r)==0)
                    potential_spin_up++;

                if(alpha==1 && sector.onsite_down(j,r)==0)
                    potential_spin_down++;
            }
            row++;
            inner_indices.push_back(i*sector.nbasis_down+j);
            if(alpha==0)
                matrix_elements.push_back(potential_spin_up);
            else
                matrix_elements.push_back(potential_spin_down);
            outer_starts.push_back(row);
            /*
            cout<<"----------------------"<<endl;
            cout<<bitset<4>(sector.id_up[i])<<" "<<bitset<4>(sector.id_down[j])<<endl;
            cout<<"matrix element:="<<potential_spin_up<<endl;
            cout<<"matrix element:="<<potential_spin_down<<endl;
            cout<<"alpha="<<alpha<<endl;
            cout<<"annil="<<annil<<endl;
            */
        }
    O.init(outer_starts,inner_indices,matrix_elements);
    outer_starts.clear();
    inner_indices.clear();
    matrix_elements.clear();
}


void lhamil::coeff_update() {
    double eigenvalues_0=1;
    double epsilon=1e-8;

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

        // checking if the iteration is converged
        if(i>10 and i%5==0) {
            diag(i);
            if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon) {
                lambda=i+1;
                break;
            }
            else
                eigenvalues_0=eigenvalues[0];
        }
    }
}

void lhamil::coeff_explicit_update()
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-8;

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

        // checking if the iteration is converged
        // the overhead of diagonalization is small
        /* 
        if(j>10 and j%5==0){
          diag(j);
          if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon){
              lambda=j+1;
              break;
          }
          else
             eigenvalues_0=eigenvalues[0];
        }
        */
        
    }
    delete phi_0,phi_1,phi_2,phi_t;
}

//Lanczos update version for spectral function calculation
void lhamil::coeff_update_wopt(vector<double> O_psir_0)
{
    int i,j,idx;
    double eigenvalues_0=1;
    double epsilon=1e-8;

    double norm_factor,overlap_factor;
    double *phi_0,*phi_1,*phi_2,*phi_t,*phi_s;
    phi_0=new double[nHilbert];
    phi_1=new double[nHilbert];
    phi_2=new double[nHilbert];
    phi_t=new double[nHilbert];

    norm.assign(lambda+1,0);
    overlap.assign(lambda,0);

    for(i=0;i<nHilbert;i++)
        phi_0[i]=O_psir_0[i];

    norm_factor=0;
    #pragma omp parallel for reduction(+:norm_factor)
    for(i=0; i<nHilbert; i++)
        norm_factor+=phi_0[i]*phi_0[i];
    norm_factor=sqrt(norm_factor)+1e-24;
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

        // checking if the iteration is converged
        // the overhead of diagonalization is small
        /*
        if(j>10 and j%5==0){
          diag(j);
          if(abs((eigenvalues[0]-eigenvalues_0)/(abs(eigenvalues_0)+1e-8))<epsilon){
              lambda=j+1;
              break;
          }
          else
             eigenvalues_0=eigenvalues[0];
        }
        */
    }

    delete phi_0,phi_1,phi_2,phi_t;
}


void lhamil::diag()
{
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

void lhamil::diag(int l)
{
    if(norm.size()==0) coeff_update();
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
    eigenvalues.assign(l,0);
    psi_0.assign(l,0);
    psi_n0.assign(l,0);
    for(int i=0; i<l; i++) {
        eigenvalues[i]=e[i];
        psi_0[i]=h[i];
        psi_n0[i]=h[i*l];
    }
    delete h,e;
}

// given the seed of phi_0, generate the Lanczos basis again, and transform the eigenstate to the desired eigenstate representation
void lhamil::eigenstates_reconstruction() {
    int l,n;
    double Norm=0;
    Vec phi_0,phi_1,phi_2;
    E0=eigenvalues[0];
    // repeat the iteration, but with normalization and overlap values at hand
    phi_0.init_random(nHilbert,seed);
    phi_1 = H * phi_0;
    phi_1 -= phi_0*overlap[0];
    phi_1 /= norm[1];
    psir_0.assign(nHilbert,0);
    for(n=0; n<nHilbert; n++)
        psir_0[n]+=psi_0[0]*phi_0.value[n]+psi_0[1]*phi_1.value[n];

    for(l=2; l<overlap.size(); l++) {
        phi_2 = H * phi_1;
        phi_2 -= phi_1 * overlap[l-1] + phi_0*norm[l-1];
        phi_2/= norm[l];
        for(n=0; n<nHilbert; n++)
            psir_0[n]+=psi_0[l]*phi_2.value[n];
        swap(&phi_0,&phi_1,&phi_2);
    }
    for(n=0; n<nHilbert; n++)
        Norm+=psir_0[n]*psir_0[n];
    for(n=0; n<nHilbert; n++)
        psir_0[n]/=sqrt(Norm);
}

void lhamil::psir0_creation_el_up(basis & sector_i,basis &sector_O,vector<double> & O_psir_0, long n)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(psir_0.size()==0) return;
    for(i=0;i<sector_i.nbasis_up;i++){
       //apply operator O on basis id_up[i]
       Ob=sector_i.creation(sector_i.id_up[i],n);
       // if the operator could be applied
       if(Ob!=sector_i.id_up[i]){
        // find the corresponding index in the _sector basis set
        it=sector_O.basis_up.find(Ob);
        if(it!=sector_O.basis_up.end()){
          k=it->second;   
          //cout<<"i,id_up[i],k,new_id_up[k]:=";
          //cout<<i<<","<<sector_i.id_up[i]<<","<<k<<","<<sector_O.id_up[k]<<endl; 
       for(j=0;j<sector_i.nbasis_down;j++) 
          O_psir_0[k*sector_i.nbasis_down+j]=psir_0[i*sector_i.nbasis_down+j];
        }
      } 
    }
}

void lhamil::psir0_creation_el_down(basis & sector_i,basis &sector_O,vector<double> &O_psir_0,long n)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(psir_0.size()==0) return;
    for(j=0;j<sector_i.nbasis_down;j++) {
       //apply operator O on basis id_down[j]
       Ob=sector_i.creation(sector_i.id_down[j],n);
       // if the operator could be applied
       if(Ob!=sector_i.id_down[j]){
           // find the corresponding index in the _sector basis set
           it=sector_O.basis_down.find(Ob);
           if(it!=sector_O.basis_down.end()){
              k=it->second;   
          for(i=0;i<sector_i.nbasis_up;i++)
           O_psir_0[i*sector_O.nbasis_down+k]=psir_0[i*sector_i.nbasis_down+j];
           }
       }
    }
}

void lhamil::psir0_annihilation_el_up(basis &sector_i,basis &sector_O,vector<double> &O_psir_0,long n)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(psir_0.size()==0) return;
    for(i=0;i<sector_i.nbasis_up;i++){
       //apply operator O on basis id_up[i]
       Ob=sector_i.annihilation(sector_i.id_up[i],n);
       // if the operator could be applied
       if(Ob!=sector_i.id_up[i]){
        // find the corresponding index in the _sector basis set
        it=sector_O.basis_up.find(Ob);
        if(it!=sector_O.basis_up.end())
          k=it->second;   
        for(j=0;j<sector_i.nbasis_down;j++) 
          O_psir_0[k*sector_i.nbasis_down+j]=psir_0[i*sector_i.nbasis_down+j];
       }
    }
}

void lhamil::psir0_annihilation_el_down(basis & sector_i,basis &sector_O,vector<double> &O_psir_0,long n)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(psir_0.size()==0) return;
    for(j=0;j<sector_i.nbasis_down;j++) {
       //apply operator O on basis id_down[j]
       Ob=sector_i.annihilation(sector_i.id_down[j],n);
       // if the operator could be applied
       if(Ob!=sector_i.id_down[j]){
         // find the corresponding index in the _sector basis set
         it=sector_O.basis_down.find(Ob);
         if(it!=sector_O.basis_down.end())
           k=it->second;   
         for(i=0;i<sector_i.nbasis_up;i++)
           O_psir_0[i*sector_O.nbasis_down+k]=psir_0[i*sector_i.nbasis_down+j];
       }
    }
}


void lhamil::print_lhamil(int range){
    for(int i=0;i<range;i++){
       if(i==0)
         cout<<"[[";
       else cout<<" [";
       for(int j=0;j<nHilbert;j++){
          if(j==i+1)
            cout<<norm[j]<<",";
          else if(j==i-1)
            cout<<norm[i]<<",";
          else if(j==i)
            cout<<overlap[i]<<",";
          else
            cout<<0<<",";
         }
       if(i==range-1)
         cout<<"]]"<<endl; 
       else cout<<"]"<<endl;
    }
}


void lhamil::print_hamil() {
    int i,j,count;
    for(i=0;i<nHilbert;i++){
       if(i==0)
         cout<<"[[";
       else cout<<" [";
       count=0;
       for(j=0;j<nHilbert;j++){
          if(j==H.inner_indices[H.outer_starts[i]+count])
            cout<<H.value[H.outer_starts[i]+count++]<<",";
          else
            cout<<0<<",";
         }
       if(i==nHilbert-1)
         cout<<"]]"<<endl; 
       else cout<<"]"<<endl;
    } 
}

void lhamil::print_eigen( int range) {
    if(range>=lambda) range=lambda;
    cout<<"Eigenvalues:= [";
    for(int i=0; i<range; i++)
        cout<<eigenvalues[i]<<", ";
    cout<<",...]"<<endl;
}

double lhamil::ground_state_energy() {
    vector<double> H_psir0;
    double overlap=0;
    if(psir_0.size()!=0) {
        H_psir0=H*psir_0;
        for(int i=0; i<nHilbert; i++)
            overlap+=psir_0[i]*H_psir0[i];
    }
    return overlap;
}


// continued fraction version
/*
double lhamil::spectral_function(double omega, double eta, int annihil) {
    // calculation continued fraction using modified Lentz method
    complex<double> E;
    if(annihil==1)
       E=complex<double>(omega+E0,eta);
    else
       E=complex<double>(omega-E0,eta);
    vector<complex<double>>  f,c,d,delta;
    complex<double> a,b,G;
    double I;
    b=E-overlap[0];
    if(abs(b)<1e-30)
        f.push_back(1e-30);
    else
        f.push_back(b);
    c.push_back(f[0]);
    d.push_back(0);
    delta.push_back(0);
    for(int n=1; n<lambda; n++) {
        b=E-overlap[n];
        a=-norm[n]*norm[n];
        d.push_back(b+a*d[n-1]);
        if(abs(d[n])<1e-30)
            d[n]=1e-30;
        c.push_back(b+a/c[n-1]);
        if(abs(c[n])<1e-30)
            c[n]=1e-30;
        d[n]=1.0/d[n];
        delta.push_back(c[n]*d[n]);
        f.push_back(f[n-1]*delta[n]);
        if(abs(delta.back()-1)<1e-15)
            break;
    }
    G=1.0/f.back();

    I=-G.imag()/M_PI;

    c.clear();
    d.clear();
    f.clear();
    delta.clear();
    return I;
}
*/

// spectrum decomposition version
double lhamil::spectral_function(double omega, double eta, int annil) {
    complex<double> E(omega,eta);
    complex<double> G=0;
    for(int i=0; i<lambda; i++)
        if(annil==1)
            G+=psi_n0[i]*psi_n0[i]/(E-(E0-eigenvalues[i]));
        else
            G+=psi_n0[i]*psi_n0[i]/(E-(-E0+eigenvalues[i]));

    return -G.imag()/M_PI;
}


complex<double> lhamil::Greens_function(double omega, double eta,int annil) {
    complex<double> E(omega,eta);
    complex<double> G=0;
    for(int i=0; i<lambda; i++)
        if(annil==1)
            G+=psi_n0[i]*psi_n0[i]/(E-E0+eigenvalues[i]);
        else
            G+=psi_n0[i]*psi_n0[i]/(E+E0-eigenvalues[i]);
    return G;
}

void lhamil::save_to_file(const char* filename) {
    if(eigenvalues.size()==0) return;
    ofstream odf;
    odf.open(filename,ofstream::out);
    odf<<"nHilbert:="<<setw(10)<<" "<<nHilbert<<endl;
    odf<<"lambda:="<<setw(10)<<" "<<lambda<<endl;
    odf<<"seed:="<<setw(10)<<" "<<seed<<endl;
    odf<<"eigenvalues:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<eigenvalues[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"norm:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<=lambda; i++)
        odf<<norm[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"overlap:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<overlap[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<psi_0[i]<<" , ";
    odf<<" ] "<<endl;
    odf<<"psi_n0:="<<setw(10)<<setprecision(8)<<" [ ";
    for(int i=0; i<lambda; i++)
        odf<<psi_n0[i]<<" , ";
    odf<<" ] "<<endl;
    odf.close();
}

void lhamil::read_from_file(const char* filename) {
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
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        eigenvalues.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<=lambda; i++) {
        idf>>value;
        idf>>buffer;
        norm.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        overlap.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        psi_0.push_back(value);
    }
    idf>>buffer;
    idf>>buffer;
    idf>>buffer;
    for(int i=0; i<lambda; i++) {
        idf>>value;
        idf>>buffer;
        psi_n0.push_back(value);
    }
    idf.close();
}
