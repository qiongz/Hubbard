#include<iostream>
#include<sstream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<unistd.h>
#include<vector>
#include<algorithm>
//#include"mkl.h"

using namespace std;

int main(int argc,char *argv[]){
    int n_site,n_spin;
    n_site=2; 
    n_spin=2;
    
    void usage(char*);
    int getopt(int argc,char * const argv[],const char *optstring);
    extern char *optarg;
    extern int optind, opterr, optopt;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:n:h:"))!=-1){
    switch(ch){
        case 'l':
            n_site=atoi(optarg);
            break;
        case 'n':
            n_spin=atoi(optarg);
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
    

    return 0;
}
void usage(char *target){
    cout<<"Usage: "<<target<<" -[Options]\n";
    cout<<"Options:\n";
    cout<<"   -l    Number of sites\n";
    cout<<"   -n    Number of spins\n";
    cout<<"Default: l=2  n=2"<<endl;
}
