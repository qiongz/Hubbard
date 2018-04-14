#include"hamiltonian.h"
hamil::hamil() {}

hamil::hamil(basis &sector,double _V,double _t, double _U){
    V=_V;
    t=_t;
    U=_U;
    init(sector,V,t,U);
}

hamil::~hamil() {}

void hamil::init(basis &sector,double _V,double _t,double _U){
    long nsite,nbasis_up,nbasis_down,signu,signd;
    V=_V;
    t=_t;
    U=_U;
    nsite=sector.nsite;
    nbasis_up=sector.nbasis_up;
    nbasis_down=sector.nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    hamiltonian=new double[nHilbert*nHilbert];
    memset(hamiltonian,0,nHilbert*nHilbert*sizeof(double));
    long n,m,i,j,k,l,s;
    for(i=0; i<nbasis_up; i++)
        for(j=0; j<nbasis_down; j++) {
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
                if(k!=i)
                    hamiltonian[(i*nbasis_down+j)*nHilbert+k*nbasis_down+j]+=-t*signu;

                l=sector.hopping_down(j,n,m);
                if(l!=j)
                    hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+l]+=-t*signd;

                if(sector.potential(i,j,n))
                    hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=U;
                if(sector.onsite_up(i,n))
                    hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=V*pow(-1,n);
                if(sector.onsite_down(j,n))
                    hamiltonian[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=V*pow(-1,n);
            }
    }
}

const hamil & hamil::operator =(const hamil & _gs_hconfig) {
    if(this !=&_gs_hconfig) {
        nHilbert=_gs_hconfig.nHilbert;
        t=_gs_hconfig.t;
        V=_gs_hconfig.V;
        U=_gs_hconfig.U;
        eigenvalues.assign(_gs_hconfig.eigenvalues.begin(),_gs_hconfig.eigenvalues.end());
        psir_0.assign(_gs_hconfig.psir_0.begin(),_gs_hconfig.psir_0.end());
        psi_n0.assign(_gs_hconfig.psi_n0.begin(),_gs_hconfig.psi_n0.end());
    }
    return *this;
}

double hamil::ground_state_energy() {
    if(psir_0.size()==0) return 0;
    double E_gs=0;
    vector<double> psi_t;
    psi_t.assign(nHilbert,0);
    for(int i=0;i<nHilbert;i++)
      for(int j=0;j<nHilbert;j++)
         psi_t[i]+=hamiltonian[i*nHilbert+j]*psir_0[j];
    for(int i=0; i<nHilbert; i++)
        E_gs+=psi_t[i]*psir_0[i];
    psi_t.clear();
    return E_gs;
}

void hamil::diag() {
    int i;
    double *h=new double[nHilbert*nHilbert];
    double *en=new double[nHilbert];
    memset(h,0,sizeof(double)*nHilbert*nHilbert);
    for(i=0;i<nHilbert*nHilbert;i++)
        h[i]=hamiltonian[i];
    diag_dsyev(h,en,nHilbert);
    psir_0.assign(nHilbert,0);
    psi_n0.assign(nHilbert,0);
    eigenvalues.assign(nHilbert,0);
    for(i=0; i<nHilbert; i++) {
        eigenvalues[i]=en[i];
        psir_0[i]=h[i];
        psi_n0[i]=h[i*nHilbert];
    }
    delete h,en;
}

void hamil::print_hamil(int range){
    int i, j, count;
    if(range>nHilbert)
        range=nHilbert;
    for(i = 0; i < range; i++) {
        if(i == 0)
            cout <<setw(2)<< "[[";
        else cout <<setw(2)<< " [";
        // count is the No. of nonzero elements in the row
        for(j=0;j<range;j++)
            cout<<setw(5)<<setprecision(2)<<hamiltonian[i*nHilbert+j]<<", ";
        if(i == range - 1)
            cout << ",...]]" << endl;
        else cout << ",...]" << endl;
    }
}

void hamil::print_eigen(int range){
    if(range>=nHilbert)
       range=nHilbert;
    std::cout << "Eigenvalues:=[ ";
    for(int i = 0; i < range; i++)
        if(i != range - 1)
            std::cout << eigenvalues[i] << ", ";
        else
            std::cout << eigenvalues[i] << " , ...]" << std::endl;
}
