#include"hamiltonian.h"
hamil::hamil() {}

hamil::hamil(basis &sector,double t, double U) {
    long nsite,nbasis_up,nbasis_down;
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
                k=sector.hopping_up(i,n,n+1);
                l=sector.hopping_down(j,n,n+1);
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
                if(sector.potential(i,j,n)) {
                    row++;
                    inner_indices.push_back(i*nbasis_down+j);
                    matrix_elements.push_back(U);
                }
            }
            if(sector.potential(i,j,nsite-1)) {
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

hamil::~hamil() {}

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

complex<double> hamil::Greens_function(double E_real,double epsilon) {
    complex<double> E=complex<double>(E_real,epsilon);
    complex<double> G=0;
    for(int i=0; i<nHilbert; i++)
        G+=psi_n0[i]*psi_n0[i]/(E-eigenvalues[i]);
    return G;
}


void hamil::print_hamil() {
    std::cout<<"hamiltonian in CSR format: "<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    H.print();
}

