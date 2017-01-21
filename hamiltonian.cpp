#include"hamiltonian.h"
double normalize(double *wf_1,int n) {
    double sum=0;
    for(int i=0; i<n; i++)
        sum+=wf_1[i]*wf_1[i];
    for(int i=0; i<n; i++)
        wf_1[i]/=sum;
    return sum;
}

double overlap(double *wf_1,double *wf_2,int n) {
    double sum=0;
    for(int i=0; i<n; i++)
        sum+=wf_1[i]*wf_2[i];
    return sum;
}

void hoperation(double *wf_1, double * wf_2, std::vector<double> &hamil,std::vector<int> & row_index, std::vector<int> & col_index,int n) {
    memset(wf_2,0,sizeof(double)*n);
    for(int i=0; i<hamil.size(); i++)
        wf_2[col_index[i]]+=hamil[i]*wf_1[row_index[i]];
}

void diag_hamil(basis *sector,double t, double U, double energy,double *wf) {
    int nsite,nbasis_up,nbasis_down,nbasis;
    nsite=(*sector).nsite;
    nbasis_up=(*sector).nbasis_up;
    nbasis_down=(*sector).nbasis_down;
    nbasis=nbasis_up*nbasis_down;
    srand(1);

    // only the nonzero hamiltonian elements with their indices are stored
    double *full_hamil,*full_eigenvalues;
    full_hamil=new double[nbasis*nbasis];
    full_eigenvalues=new double[nbasis];
    memset(full_hamil,0,sizeof(double)*nbasis*nbasis);

    std::map<int,double> hamil_nonzero;
    std::vector<int> row_index,col_index;
    std::vector<double> hamil;
    int n,i,j,k,l;
    for(i=0; i<nbasis_up; i++) {
        for(j=0; j<nbasis_down; j++) {
            for(n=0; n<nsite-1; n++) {
                k=(*sector).hopping_up(i,n);
                l=(*sector).hopping_down(j,n);
                if(k!=i) {
                    hamil_nonzero[(i*nbasis_down+j)*nbasis+k*nbasis_down+j]+=-t;
                    full_hamil[(i*nbasis_down+j)*nbasis+k*nbasis_down+j]+=-t;
                }
                if(l!=j) {
                    hamil_nonzero[(i*nbasis_down+j)*nbasis+i*nbasis_down+l]+=-t;
                    full_hamil[(i*nbasis_down+j)*nbasis+i*nbasis_down+l]+=-t;
                }
            }
            for(n=0; n<nsite; n++)
                if((*sector).potential(i,j,n)) {
                    hamil_nonzero[(i*nbasis_down+j)*nbasis+i*nbasis_down+j]+=U;
                    full_hamil[(i*nbasis_down+j)*nbasis+i*nbasis_down+j]+=U;
                }
        }
    }
    //std::cout<<"hamiltonain nonzero elements:"<<std::endl;
    for(auto &x:hamil_nonzero) {
        //    std::cout<<x.first/nbasis<<" "<<x.first%nbasis<<" "<<x.second<<std::endl;
        row_index.push_back(x.first/nbasis);
        col_index.push_back(x.first%nbasis);
        hamil.push_back(x.second);
    }
    hamil_nonzero.clear();

    /* Lanczos basis generation */
    int lambda=20;
    double *norm_factor,*overlap_factor;
    double **lanczos_basis;
    norm_factor=new double[lambda+1];
    overlap_factor=new double[lambda+1];
    lanczos_basis=new double *[lambda+1];
    for(i=0; i<lambda+1; i++)
        lanczos_basis[i]=new double[nbasis];

    // initialize phi_0
    for(j=0; j<nbasis; j++)
        lanczos_basis[0][j]=(rand()*1.0/RAND_MAX-0.5);
    norm_factor[0]=normalize(lanczos_basis[0],nbasis);

    // initialize phi_1
    hoperation(lanczos_basis[0],lanczos_basis[1],hamil,row_index,col_index,nbasis);
    overlap_factor[0]=overlap(lanczos_basis[0],lanczos_basis[1],nbasis);
    for(j=0; j<nbasis; j++)
        lanczos_basis[1][j]=lanczos_basis[1][j]-overlap_factor[0]*lanczos_basis[0][j];
    norm_factor[1]=normalize(lanczos_basis[1],nbasis);

    // iterative generation of basis of phi_1,phi_2,...,phi_lambda
    int lambda_cut=1;
    for(i=1; i<lambda; i++) {
        if(lambda_cut>1)
           break;
        hoperation(lanczos_basis[i],lanczos_basis[i+1],hamil,row_index,col_index,nbasis);
        overlap_factor[i]=overlap(lanczos_basis[i],lanczos_basis[i+1],nbasis);
        for(j=0; j<nbasis; j++)
            lanczos_basis[i+1][j]=lanczos_basis[i+1][j]-overlap_factor[i]*lanczos_basis[i][j]-norm_factor[i]*lanczos_basis[i-1][j];
        norm_factor[i+1]=normalize(lanczos_basis[i+1],nbasis);
        if(fabs(overlap_factor[i]/norm_factor[i])<1e-6)
           lambda_cut=i;
    }
 //   lambda=lambda_cut;
    // hamiltonian in the Lanczos basis
    double *lanczos_eigenvalues=new double[lambda];
    double *lanczos_hamil=new double [lambda*lambda];
    memset(lanczos_hamil,0,sizeof(double)*lambda*lambda);

    for(i=0;i<lambda;i++){
       lanczos_hamil[i*lambda+i+1]=sqrt(norm_factor[i+1]/norm_factor[i]);
       lanczos_hamil[(i+1)*lambda+i]=lanczos_hamil[i*lambda+i+1];
       lanczos_hamil[i*lambda+i]=overlap_factor[i]/norm_factor[i];
    }

    
    // print the lanczos hamiltonian
    std::cout<<"#lanczos hamiltonian"<<std::endl;
    std::cout<<setprecision(2)<<setw(6);
    for(i=0; i<lambda; i++) {
        if(i==0)
            std::cout<<"[[ ";
        else
            std::cout<<" [ ";
        for(j=0; j<lambda; j++)
            std::cout<<lanczos_hamil[i*lambda+j]<<", ";
        if(i==lambda-1)
            std::cout<<"]]"<<std::endl;
        else
            std::cout<<"] "<<std::endl;
    }

    // print the full hamiltonian
    /*
    std::cout<<"#full hamiltonian"<<std::endl;
    for(i=0;i<nbasis;i++){
        if(i==0)
            std::cout<<"[[ ";
        else
            std::cout<<" [ ";
        for(j=0; j<nbasis; j++)
            std::cout<<full_hamil[i*nbasis+j]<<", ";
        if(i==nbasis-1)
            std::cout<<"]]"<<std::endl;
        else
            std::cout<<"] "<<std::endl;
    }
    */

    // diagonalization and compare
    diag(full_hamil,full_eigenvalues,nbasis);
    diag(lanczos_hamil,lanczos_eigenvalues,lambda);
    // print and compare the eigenvalues

    for(i=0; i<lambda; i++)
        std::cout<<full_eigenvalues[i]<<" "<<lanczos_eigenvalues[i]<<std::endl;

    delete [] full_hamil,full_eigenvalues,norm_factor,overlap_factor;
    for(i=0; i<lambda+1; i++)
        delete [] lanczos_basis[i];
    delete [] lanczos_basis;
    delete [] lanczos_eigenvalues,lanczos_hamil;

}


