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

    nsite=2;
    nel=2;
    t=1.0;
    U=1;
    lambda=100;
    init_argv(nsite,nel,t,U,lambda,argc,argv);

    nel_up=nel/2;
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    nel_down=nel-nel_up;
    basis sector(nsite,nel_up,nel_down);
    sector.init();

    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    lhamil config(sector,t,U,lambda,seed);
    config.coeff_explicit_update();
    config.diag();
    config.eigenstates_reconstruction();

    //config.save_to_file(filename.c_str());
    //config.eigenstates_reconstruction();

    // for spin-up
    config.set_onsite_optc(0,0,1);
    config.coeff_update_wopt();
    config.diag();
    vector<double> E,I;
    for(int i=0; i<200; i++) {
        E.push_back(i/10.0-10);
        I.push_back(config.spectral_function(E[i],0.05));
    }
    config.set_onsite_optc(0,0,0);
    config.coeff_update_wopt();
    config.diag();
    for(int i=0; i<200; i++) {
        I[i]+=config.spectral_function(E[i],0.05);
    }

    // for spin-down
    config.set_onsite_optc(0,1,1);
    config.coeff_update_wopt();
    config.diag();
    for(int i=0; i<200; i++) {
        I[i]+=config.spectral_function(E[i],0.05);
    }
    config.set_onsite_optc(0,1,0);
    config.coeff_update_wopt();
    config.diag();
    for(int i=0; i<200; i++) {
        I[i]+=config.spectral_function(E[i],0.05);
    }

    for(int i=0; i<200; i++)
        cout<<E[i]<<" "<<I[i]<<endl;


    /*
    //test for eigenvalues convergence
    for(int l=5;l<lambda;l++){
       config.diag(l);
       //config.eigenstates_reconstruction();
       cout<<l<<" ";
       for(int i=0;i<5;i++)
         cout<<config.eigenvalues[i]<<" ";
       cout<<endl;
    }
    */

    return 0;
}
