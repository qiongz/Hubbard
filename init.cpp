#include<unistd.h>
#include"init.h"
void usage(char *target){
    std::cout<<"Usage: "<<target<<" [Options]\n";
    std::cout<<"Options:\n";
    std::cout<<"  -l                       Number of sites\n";
    std::cout<<"  -n                       Number of electrons\n";
    std::cout<<"  -t                       Hopping strength\n";
    std::cout<<"  -U                       Onsite replusion energy\n";
    std::cout<<"Default: (l,n,t,U) = (2,2,1.0,0.5)\n";
}
void init(int* n_site, int* n_el, double *t, double *U, int argc,char *argv[]){
    int getopt(int argc,char * const argv[],const char *optstring);
    extern char *optarg;
    extern int optind, opterr, optopt;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:n:t:U:h:"))!=-1){
    switch(ch){
        case 'l':
            *n_site=atoi(optarg);
            break;
        case 'n':
            *n_el=atoi(optarg);
            break;
        case 't':
            *t=atof(optarg);
            break;
        case 'U':
            *U=atof(optarg);
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
}
