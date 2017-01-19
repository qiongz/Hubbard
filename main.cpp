#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<unistd.h>
#include<vector>
#include<list>
#include<map>
#include<algorithm>
#include"basis.h"

using namespace std;

int main(int argc,char *argv[]){
    int n_site,n_el;
    int n_elu,n_eld;
    int n_basis,nb_up,nb_down;
    int Sz_tot;
    double t,U;
    double *hamiltonian;
    double *energy;

    /* initialize parameters */
    n_site=2; 
    n_el=2;
    t=1.0;
    U=0.5;
    void usage(char*);
    int getopt(int argc,char * const argv[],const char *optstring);
    extern char *optarg;
    extern int optind, opterr, optopt;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:n:t:U:h:"))!=-1){
    switch(ch){
        case 'l':
            n_site=atoi(optarg);
            break;
        case 'n':
            n_el=atoi(optarg);
            break;
        case 't':
            t=atof(optarg);
            break;
        case 'U':
            U=atof(optarg);
            break;
        case 'h':
            errFlag++;
            break;
        default:
            errFlag++;
            break;
       }
    }
    if(errFlag) {
        usage(argv[0]);
        exit(2);
    }

    /* generating basis */
    int n,i,j,k,l;
    // Sz=0 sector 
    n_elu=n_el/2;
    n_eld=n_el/2;
    basis sector(n_site,n_elu,n_eld);
    /* print basis set */
    sector.print();


    /* calculating hamiltonian matrix elements */
    n_basis=sector.n_basis;
    nb_up=sector.nb_up;
    nb_down=sector.nb_up;
    hamiltonian=new double[n_basis*n_basis];
    energy=new double[n_basis];
    memset(hamiltonian,0,16*n_basis*n_basis);
    for(i=0;i<nb_up;i++){
        for(j=0;j<nb_down;j++){
            for(n=0;n<n_site-1;n++){
                k=sector.hopping_up(i,n);
                l=sector.hopping_down(j,n);
                if(k!=i)
                    hamiltonian[(i*nb_up+j)*n_basis+k*nb_up+j]+=-t;
                if(l!=j)
                    hamiltonian[(i*nb_up+j)*n_basis+i*nb_up+l]+=-t;
            }
            for(n=0;n<n_site;n++)
                if(sector.potential(i,j,n))
                    hamiltonian[(i*nb_up+j)*n_basis+i*nb_up+j]+=U;
        }
    }

    /* print hamiltonian matrix */
    for(i=0;i<n_basis;i++){
        if(i==0)
           cout<<"[[ ";
        else
           cout<<"[ ";
        for(j=0;j<n_basis;j++)
           cout<<hamiltonian[i*n_basis+j]<<", ";
        if(i==n_basis-1)
            cout<<"]]"<<endl;
        else
            cout<<"]"<<endl;
    }
         
    return 0;
}
void usage(char *target){
    cout<<"Usage: "<<target<<" [Options]\n";
    cout<<"Options:\n";
    cout<<"  -l                       Number of sites\n";
    cout<<"  -n                       Number of electrons\n";
    cout<<"  -t                       Hopping strength\n";
    cout<<"  -U                       Onsite replusion energy\n";
    cout<<"Default: (l,n,t,U) = (2,2,1.0,0.5)"<<endl;
}
