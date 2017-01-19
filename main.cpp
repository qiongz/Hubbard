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
    int Sz_tot;
    double t,U;

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
    n_elu=n_el/2;
    n_eld=n_el/2;
    basis sector(n_site,n_elu,n_eld);
    sector.print();
    
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
