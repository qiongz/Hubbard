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
    U=5;
    lambda=100;
    init_argv(nsite,nel,t,U,lambda,argc,argv);

    nel_up=nel/2;
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    nel_down=nel-nel_up;
    basis sector(nsite,nel_up,nel_down);
    sector.init();
    //sector.prlong();


    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    lhamil config(sector,t,U,lambda,seed);
    config.coeff_explicit_update();
    config.diag();
    config.eigenstates_reconstruction();

    //cout<<"the hamiltonian matrix:="<<endl;  
    //config.print_hamil();
    
    //config.save_to_file(filename.c_str());
    //config.eigenstates_reconstruction();
    /* 
    cout<<"gs energy per site:="<<config.ground_state_energy()/nsite<<endl;
    cout<<"gs energy tot:="<<config.ground_state_energy()<<endl;
    double analytic_gs=-4*t/M_PI+U/4.0-0.017*U*U/t;
    cout<<"analytical gs energy(band limit):="<<analytic_gs<<endl;

    analytic_gs=-4*t*t*log(2)/U;
    cout<<"analytical gs energy(atomic limit):="<<analytic_gs<<endl;
    */

    // for spin-up
    vector<double> I;
    I.assign(200,0);
       
    int r=0;
    basis sector_annihilation=sector;
    sector_annihilation.annihilation_el_up(r);
    vector<double> O_psir_0;
    config.psir0_annihilation_el_up(O_psir_0,r);
    config.set_hamil(sector_annihilation,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    cout<<"mu_hole:="<<config.eigenvalues[0]-config.E0<<endl;
    for(int i=0; i<200; i++) 
        I[i]+=config.spectral_function(i/4.0-25,0.1,1);
    O_psir_0.clear(); 
    
    basis sector_creation=sector;
    sector_creation.creation_el_up(r);
    config.psir0_creation_el_up(O_psir_0,r);
    config.set_hamil(sector_creation,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    cout<<"mu_particle:="<<config.eigenvalues[0]-config.E0<<endl;
    for(int i=0; i<200; i++) 
        I[i]+=config.spectral_function(i/4.0-25,0.1,0);
    O_psir_0.clear(); 
    
    // for spin-down
    /*
    basis sector_annihilation_down=sector;
    sector_annihilation_down.annihilation_el_down(r);
    config.psir0_annihilation_el_down(O_psir_0,r);
    config.set_hamil(sector_annihilation_down,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    for(int i=0; i<200; i++) 
        I[i]+=config.spectral_function(i/4.0-25,0.1,1);
    O_psir_0.clear(); 
    
    basis sector_creation_down=sector;
    sector_creation_down.creation_el_down(r);
    config.psir0_creation_el_down(O_psir_0,r);
    config.set_hamil(sector_creation_down,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    for(int i=0; i<200; i++) 
        I[i]+=config.spectral_function(i/4.0-25,0.1,0);
    O_psir_0.clear(); 
    */ 
    // for(int i=0; i<200; i++)
    //     cout<<i/4.0-12.5+config.E0<<" "<<I[i]<<endl;
    

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
