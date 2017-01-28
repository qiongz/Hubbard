#include"hamiltonian.h"
void diag_hamil(basis *sector,double t, double U, double energy,double *wf) {
    int nsite,nbasis_up,nbasis_down,nHilbert;
    nsite=sector->nsite;
    nbasis_up=sector->nbasis_up;
    nbasis_down=sector->nbasis_down;
    nHilbert=nbasis_up*nbasis_down;
    srand(1);

    // only the nonzero hamiltonian elements with their indices are stored
    double *full_hamil,*full_eigenvalues;
    full_hamil=new double[nHilbert*nHilbert];
    full_eigenvalues=new double[nHilbert];
    memset(full_hamil,0,sizeof(double)*nHilbert*nHilbert);

    std::map<int,double> hamil_nonzero;
    std::vector<int> row_index,col_index;
    std::vector<double> hamil;
    int n,i,j,k,l;
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            for(n=0; n<nsite-1; n++) {
                k=sector->hopping_up(i,n);
                l=sector->hopping_down(j,n);
                if(k!=i) {
                    hamil_nonzero[(i*nbasis_down+j)*nHilbert+k*nbasis_down+j]+=-t;
                    full_hamil[(i*nbasis_down+j)*nHilbert+k*nbasis_down+j]+=-t;
                }
                if(l!=j) {
                    hamil_nonzero[(i*nbasis_down+j)*nHilbert+i*nbasis_down+l]+=-t;
                    full_hamil[(i*nbasis_down+j)*nHilbert+i*nbasis_down+l]+=-t;
                }
            }
            for(n=0; n<nsite; n++)
                if(sector->potential(i,j,n)) {
                    hamil_nonzero[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=U;
                    full_hamil[(i*nbasis_down+j)*nHilbert+i*nbasis_down+j]+=U;
                }
        }
    }
    //std::cout<<"hamiltonain nonzero elements:"<<std::endl;
    for(auto &x:hamil_nonzero) {
        // std::cout<<x.first/nHilbert<<" "<<x.first%nHilbert<<" "<<x.second<<std::endl;
        row_index.push_back(x.first/nHilbert);
        col_index.push_back(x.first%nHilbert);
        hamil.push_back(x.second);
    }
    hamil_nonzero.clear();

    /* Lanczos basis generation */
    int lambda,lambda_max=100;
    vector<double> norm_factor,overlap_factor; 
    vector<lbasis> basis_list;
    
    // initialize phi_0
    lbasis phi_0(nHilbert);
    phi_0.init_random();
    basis_list.push_back(phi_0);
    norm_factor.push_back(1);
    
    // initialize H|phi_0> 
    lbasis phi_1=phi_0.hoperation(hamil,row_index,col_index);
    // a_0=<phi_0|H|phi_0>
    overlap_factor.push_back(phi_0*phi_1);
    // |phi_1>=H|phi_0>-a_0|phi_0>
    phi_1=phi_1-phi_0*overlap_factor[0];
    // |phi_1>=|phi_1>/N_1
    norm_factor.push_back(phi_1.normalize());
    basis_list.push_back(phi_1);

    // iterative generation of basis of phi_1,phi_2,...,phi_lambda
    for(i=1; i<lambda_max; i++) {
        // initialize H|phi_i>
        phi_1=basis_list[i].hoperation(hamil,row_index,col_index);
        // a_i=<phi_i|H|phi_i>
        overlap_factor.push_back(basis_list[i]*phi_1); 
        // |phi_i+1>=H|phi_i>-a_i|phi_i>-N_i|phi_i-1>
        phi_1=phi_1-basis_list[i]*overlap_factor[i]-basis_list[i-1]*norm_factor[i];
        // |phi_i+1>=|phi_i+1>/N_i+1
        norm_factor.push_back(phi_1.normalize());
        basis_list.push_back(phi_1);
        if(fabs(overlap_factor[i]/norm_factor[i])<1e-8){
           lambda=i-1;
           break;
        }
    }
    lambda=4; 
    cout<<"# Lanczos basis, norm_factor,  overlap_factor"<<endl;
    for(i=0;i<lambda;i++)
        cout<<basis_list[i]<<" "<<setprecision(6)<<setw(12)<<norm_factor[i]<<" "<<overlap_factor[i]<<endl;
    
  
    // hamiltonian in the Lanczos basis
    double *lanczos_eigenvalues=new double[lambda];
    double *lanczos_hamil=new double [lambda*lambda];
    memset(lanczos_hamil,0,sizeof(double)*lambda*lambda);

    for(i=0;i<lambda-1;i++){
       lanczos_hamil[i*lambda+i+1]=sqrt(norm_factor[i+1]);
       lanczos_hamil[(i+1)*lambda+i]=sqrt(norm_factor[i+1]);
       lanczos_hamil[i*lambda+i]=overlap_factor[i];
    }
    lanczos_hamil[lambda*(lambda-1)+lambda-1]=overlap_factor[lambda-1];
    
    // print the lanczos hamiltonian
    std::cout<<"#lanczos hamiltonian"<<std::endl;
    for(i=0; i<lambda; i++) {
        if(i==0)
            std::cout<<"[[ ";
        else
            std::cout<<" [ ";
        for(j=0; j<lambda-1; j++)
            std::cout<<setprecision(3)<<setw(8)<<lanczos_hamil[i*lambda+j]<<", ";
        std::cout<<setprecision(3)<<setw(8)<<lanczos_hamil[i*lambda+lambda-1]<<" ";
        if(i==lambda-1)
            std::cout<<"]]"<<std::endl;
        else
            std::cout<<"] "<<std::endl;
    }
    // print the full hamiltonian
    std::cout<<"#full hamiltonian"<<std::endl;
    for(i=0;i<nHilbert;i++){
        if(i==0)
            std::cout<<"[[ ";
        else
            std::cout<<" [ ";
        for(j=0; j<nHilbert-1; j++)
            std::cout<<setprecision(3)<<setw(8)<<full_hamil[i*nHilbert+j]<<", ";
        std::cout<<setprecision(3)<<setw(8)<<full_hamil[i*nHilbert+nHilbert-1]<<" ";
        if(i==nHilbert-1)
            std::cout<<"]]"<<std::endl;
        else
            std::cout<<"] "<<std::endl;
    }

    // diagonalization and compare
    diag(full_hamil,full_eigenvalues,nHilbert);
    diag(lanczos_hamil,lanczos_eigenvalues,lambda);
    // print and compare the eigenvalues
    std::cout<<"# Full-Eigenvalues:= [";
    for(i=0; i<nHilbert-1; i++)
        std::cout<<full_eigenvalues[i]<<", ";
    std::cout<<full_eigenvalues[i]<<" ]"<<std::endl;

    std::cout<<"# Lanczos-Eigenvalues:= [";
    for(i=0; i<lambda-1; i++)
        std::cout<<lanczos_eigenvalues[i]<<", ";
    std::cout<<lanczos_eigenvalues[i]<<" ]"<<std::endl;

    delete [] full_hamil,full_eigenvalues,norm_factor,overlap_factor;
    delete [] lanczos_eigenvalues,lanczos_hamil;
}


