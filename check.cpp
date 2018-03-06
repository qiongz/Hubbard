#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include<gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include<ctime>
#include<sstream>
#include<cstdlib>
#if __cplusplus > 199711L
#include<chrono>
#endif

using namespace std;

int main(int argc,char *argv[]) {
    int nsite,nel,lambda;
    int nel_up,nel_down;
    int Sz;
    unsigned seed;
    double t,U;
    double func_Lieb_Wu(double x, void *params);
    double Integrate_Lieb_Wu(double U);

    nsite=2;
    nel=2;
    t=1.0;
    U=5;
    lambda=200;
    init_argv(nsite,nel,t,U,lambda,argc,argv);

    nel_up=(nel+1)/2;
    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    nel_down=nel-nel_up;
    Sz=nel_up-nel_down;
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

    
    /* basis check */
    /*
    config.print_hamil();
    config.print_lhamil(6);
    config.print_eigen(6);
    
    
    onfig.save_to_file(filename.c_str());
    */
     
    
    /* ground-state energy check */
    //cout<<"#U  Lanczos Lieb-Wu "<<endl;
    //cout<<U<<" "<<config.ground_state_energy()/nsite<<" "<<-Integrate_Lieb_Wu(U)<<endl;  
    
    
    /*
    cout<<"O_psir_0:=[";
    for(int i=0;i<(config.psir_0).size();i++)
       cout<<config.psir_0[i]<<",";
    cout<<"]"<<endl;
    */ 
   

    // for spin-up
    vector<double> IP,IH,EH,EP;
    vector<double> O_psir_0;
    IP.assign(1000,0);
    IH.assign(1000,0);
    EP.assign(1000,0);
    EH.assign(1000,0);
       
       
    int r=nsite/2;
    basis sector_annihilation(nsite,nel_up-1,nel_down);
    sector_annihilation.init();
    //sector_annihilation.prlong();

    O_psir_0.assign(sector_annihilation.nbasis_up*sector_annihilation.nbasis_down,0);
    config.psir0_annihilation_el_up(sector,sector_annihilation,O_psir_0,r);
    config.set_hamil(sector_annihilation,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    
    //config.print_hamil();
    //config.print_lhamil(6);
    //config.print_eigen(6);
    //cout<<"E0_N:="<<config.eigenvalues[0]<<endl;
    //cout<<"mu_N:="<<config.E0-config.eigenvalues[0]<<endl;
    double mu_N=config.E0-config.eigenvalues[0]; 
    for(int i=0; i<1000; i++) {
        IH[i]+=config.spectral_function(i/50.0-5,0.2,1);
        EH[i]=i/50.0-5;
    }
    O_psir_0.clear(); 
    
    basis sector_creation(nsite,nel_up+1,nel_down);
    sector_creation.init();
    //sector_creation.prlong();
    O_psir_0.assign(sector_creation.nbasis_up*sector_creation.nbasis_down,0);
    config.psir0_creation_el_up(sector,sector_creation,O_psir_0,r);
    config.set_hamil(sector_creation,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();

    
    //config.print_hamil();
    //config.print_lhamil(6);
    //config.print_eigen(6);
    //cout<<"E0_N+1:="<<config.eigenvalues[0]<<endl;
   //cout<<"mu_N+1:="<<config.eigenvalues[0]-config.E0<<endl;
    double mu_Np=config.eigenvalues[0]-config.E0; 
    for(int i=0; i<1000; i++) {
       EP[i]=i/50.0-5;
       IH[i]+=config.spectral_function(i/50.0-5,0.2,0);
    }

    O_psir_0.clear(); 
    
    // for spin-down
    /* 
    basis sector_annihilation_down(nsite,nel_up,nel_down-1);
    sector_annihilation_down.init();

    O_psir_0.assign(sector_annihilation_down.nbasis_up*sector_annihilation_down.nbasis_down,0);
    config.psir0_annihilation_el_down(sector,sector_annihilation_down,O_psir_0,r);
    config.set_hamil(sector_annihilation_down,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    for(int i=0; i<200; i++) 
        I[i]+=config.Greens_function(i/4.0-25,0.1,1);
    O_psir_0.clear(); 
    
    basis sector_creation_down(nsite,nel_up,nel_down+1);
    sector_creation_down.init();
    O_psir_0.assign(sector_creation_down.nbasis_up*sector_creation_down.nbasis_down,0);
    config.psir0_creation_el_down(sector,sector_creation_down,O_psir_0,r);
    config.set_hamil(sector_creation_down,t,U);
    config.coeff_update_wopt(O_psir_0);
    config.diag();
    for(int i=0; i<200; i++) 
        I[i]+=config.Greens_function(i/4.0-25,0.1,0);
    O_psir_0.clear(); 
    */
     
    for(int i=0; i<1000; i++)
     //cout<<EH[i]<<" "<<IH[i]<<" "<<EP[i]<<" "<<IP[i]<<endl;
        cout<<EH[i]-(mu_Np+mu_N)/2.0<<" "<<IH[i]<<endl;

    

    return 0;
}
double func_Lieb_Wu(double x, void *params) {
    double U = *(double *)params;
    double ret_value = 4*gsl_sf_bessel_J0(x)*gsl_sf_bessel_J1(x)/(1+exp(0.5*x*U))/(x+1e-8);
    return ret_value;
}
double Integrate_Lieb_Wu(double U) {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(20000);
    double result,error;

    gsl_function F;
    F.function = &func_Lieb_Wu;
    F.params= & U;

    gsl_integration_qagiu(& F,0,0, 1e-3, 100,w, &result, & error);
    return result;
}
