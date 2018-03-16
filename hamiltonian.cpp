#include"hamiltonian.h"
hamil::hamil() {}

hamil::hamil(basis &sector,double _t, double _U) {
    long nsite,nbasis_up,nbasis_down,signu,signd;
    t=_t;
    U=_U;
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

hamil::~hamil() {}

void hamil::init(basis &sector,double _t,double _U){
    long nsite,nbasis_up,nbasis_down,signu,signd;
    t=_t;
    U=_U;
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

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this !=&_gs_hconfig) {
        seed=_gs_hconfig.seed;
        nHilbert=_gs_hconfig.nHilbert;
        H=_gs_hconfig.H;
        t=_gs_hconfig.t;
        U=_gs_hconfig.U;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(),_gs_hconfig.eigenvalues.end());
        psi_0.assign(_gs_hconfig.psi_0.begin(),_gs_hconfig.psi_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(),_gs_hconfig.psi_n0.end());
    }
    return *this;
}


double hamil::spectral_function(vector<double> &O_phi_0,double omega,double _E0, double eta, int annil) {
    complex<double> E;
    complex<double> G=0;
    for(int i=0; i<nHilbert; i++)
        // set annil==1, which gives hole-sector
        if(annil==1){
            E=complex<double>(omega,eta);
            G+=pow(psi_n0[i]*O_phi_0[i],2)/(E+eigenvalues[i]-_E0);
          }
        // else particle-sector
        else{
            E=complex<double>(omega,eta);
            G+=pow(psi_n0[i]*O_phi_0[i],2)/(E+_E0-eigenvalues[i]);
          }

    return -G.imag()/M_PI;
}

double hamil::ground_state_energy() {
    if(psi_0.size()==0) return 0;
    double E_gs=0;
    vector<double> psi_t;
    psi_t=H*psi_0;
    for(int i=0; i<nHilbert; i++)
        E_gs+=psi_t[i]*psi_0[i];
    return E_gs;
}

void hamil::diag() {
    int i,idx;
    double *hamiltonian=new double[nHilbert*nHilbert];
    double *en=new double[nHilbert];
    memset(hamiltonian,0,sizeof(double)*nHilbert*nHilbert);
    for(i=0; i<H.outer_starts.size()-1; i++)
        for(idx=H.outer_starts[i]; idx<H.outer_starts[i+1]; idx++)
            hamiltonian[i*nHilbert+H.inner_indices[idx]]=H.value[idx];
    diag_dsyev(hamiltonian,en,nHilbert);
    psi_0.assign(nHilbert,0);
    psi_n0.assign(nHilbert,0);
    eigenvalues.assign(nHilbert,0);
    for(i=0; i<nHilbert; i++) {
        eigenvalues[i]=en[i];
        psi_0[i]=hamiltonian[i];
        psi_n0[i]=hamiltonian[i*nHilbert];
    }
    delete hamiltonian,en;
}


void hamil::print_hamil() {
    std::cout<<"hamiltonian in CSR format: "<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    H.print();
}

void hamil::print_eigen(){
  std::cout<<"Eigenvalues:=[ ";
   for(int i=0;i<nHilbert;i++)
      if(i!=nHilbert-1)
        std::cout<<eigenvalues[i]<<", ";
      else
         std::cout<<eigenvalues[i]<<" ]"<<std::endl;
}
