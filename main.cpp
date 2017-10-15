#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<ctime>
#include<chrono>
#include<sstream>
#include<cstdlib>

using namespace std;

int main(int argc,char *argv[]) {
    int nsite,nel,lambda;
    int nel_up,nel_down;
    int Sz;
    unsigned seed;
    double t,U;

    /* initialize parameters */
    nsite=2;
    nel=2;
    t=1.0;
    U=5;
    lambda=100;
    
    init_argv(nsite,nel,t,U,argc,argv);

    /* generating basis */
   // for(nel_up=1; nel_up<=nel/2; nel_up++) {
    nel_up=nel/2;
    
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        nel_down=nel-nel_up;
        Sz=(nel_up-nel_down)/2;
        //cout<<"Sz=: "<<Sz<<endl;
        /* generating basis */
        basis sector(nsite,nel_up,nel_down);
        sector.init();
        /* print basis set */
        //sector.print();
        //cout<<"sector size=: "<<sector.nbasis_up*sector.nbasis_down<<endl;
        lambda=sqrt(sector.nbasis_up*sector.nbasis_down)*2;


        stringstream sf;
        string filename;
        sf<<"sector_"<<Sz;
        sf>>filename;

        lhamil config(sector,t,U,lambda,seed);
        config.coeff_explicit_update();
        config.diag();
        //config.save_to_file(filename.c_str());
        double E;
        complex<double> G;
        for(int i=0;i<200;i++){
           E=-10+i/10.0;
           G=config.Greens_function(E,0.02); 
           cout<<E<<" "<<G.real()<<" "<<G.imag()<<endl;
        }
    //}
    return 0;
}
