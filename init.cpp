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

void init_argv(int& nsite, int& nel, double &t, double &U, int argc,char *argv[])
{
    extern char *optarg;
    int ch,errFlag;
    errFlag=0;
    while((ch=getopt(argc,argv,"l:n:t:U:h:"))!=-1){
        switch(ch){
        case 'l':
            nsite=atoi(optarg);
            break;
        case 'n':
            nel=atoi(optarg);
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
     try{  
     if(nsite<=0)
        throw std::logic_error("-l: positive value required !");
     if(nel<=0)
        throw std::logic_error("-n: positive value required !");
    if(nsite*2<nel)
        throw std::logic_error("-l -n: electron filling nu=(n/2l)>1 !");
     if(fabs(t)<1e-8 && fabs(U)<1e-8)
        throw std::logic_error("-t -U: at least one finite coupling constant required !");
     }catch(std::logic_error &e){
        std::cout<<e.what()<<std::endl;
        usage(argv[0]);
        exit(2);
     }
    
     if(errFlag){
        usage(argv[0]);
        exit(0);
     }
}

