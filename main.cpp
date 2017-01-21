#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include<ctime>
#include<cstdlib>

using namespace std;

int main(int argc,char *argv[]) {
    int nsite,nel;
    int nel_up,nel_down;
    int nbasis,nbasis_up,nbasis_down;
    int Sz;
    double t,U;
    double *wf,energy;

    /* initialize parameters */
    nsite=2;
    nel=2;
    t=1.0;
    U=0.5;
    init_argv(&nsite,&nel,&t,&U,argc,argv);

    /* generating basis */
    int n,i,j,k,l;
    for(nel_up=1; nel_up<=nel/2; nel_up++) {
        nel_down=nel-nel_up;
        Sz=(nel_up-nel_down)/2;
        cout<<"Sz=: "<<Sz<<endl;
        /* generating basis */
        basis sector(nsite,nel_up,nel_down);
        sector.init();
        /* print basis set */
        sector.print();

        diag_hamil(&sector,t,U,energy,wf);
        /* calculating hamiltonian matrix elements */
        /* 
        nb_up=sector.nb_up;
        nb_down=sector.nb_up;
        n_basis=nb_up*nb_down;
        double *hamiltonian=new double[n_basis*n_basis];
       // energy=new double[n_basis];
        memset(hamiltonian,0,sizeof(double)*n_basis*n_basis);
        for(i=0;i<nb_up;i++){
            for(j=0;j<nb_down;j++){
                for(n=0;n<n_site-1;n++){
                    k=sector.hopping_up(i,n);
                    l=sector.hopping_down(j,n);
                    if(k!=i)
                        hamiltonian[(i*nb_down+j)*n_basis+k*nb_down+j]+=-t;
                    if(l!=j)
                        hamiltonian[(i*nb_down+j)*n_basis+i*nb_down+l]+=-t;
                }
                for(n=0;n<n_site;n++)
                    if(sector.potential(i,j,n))
                        hamiltonian[(i*nb_down+j)*n_basis+i*nb_down+j]+=U;
            }
        }
        
        */
        /* print hamiltonian matrix */
        /* 
        for(i=0;i<n_basis;i++){
            cout<<setw(3);
            if(i==0)
               cout<<"[[ ";
            else
               cout<<" [ ";
            for(j=0;j<n_basis;j++)
               cout<<hamiltonian[i*n_basis+j]<<", ";
            if(i==n_basis-1)
                cout<<"]]"<<endl;
            else
                cout<<"] "<<endl;
        }
        */ 

        /* print eigenvalues */
        //diag(hamiltonian,energy,n_basis);
        /*
        cout<<"Eigenenergies:"<<endl;
        for(n=0;n<n_basis;n++)
            cout<<setprecision(3)<<energy[n]<<endl;
        */
        //cout<<"E0=:"<<energy[0]<<endl;
        // print the wavefunction corresponding to the lowest eigenenergy in the sector
        /*
        cout<<"# basis "<<"spin-up"<<" "<<"spin-down"<<" "<<"coefficients"<<endl;
        for(i=0;i<nb_up;i++)
            for(j=0;j<nb_down;j++)
                cout<<i*nb_down+j<<" "<<bitset<8>(sector.relu[i]).to_string()<<" "<<bitset<8>(sector.reld[j]).to_string()<<" "<<hamiltonian[i*nb_down+j]<<endl;

        delete [] hamiltonian;
        delete [] energy;
        */
    }
    return 0;
}
