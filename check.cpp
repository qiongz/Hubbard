#include"init.h"
#include"basis.h"
#include"hamiltonian.h"
#include"lanczos_hamiltonian.h"
#include"Greens_function.h"
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
    double t,U,V,K,mu;
    double func_Lieb_Wu(double x, void *params);
    double Integrate_Lieb_Wu(double U);

    nsite=10;
    nel=10;
    t=1;
    V=0;
    U=5;
    K=0;
    lambda=200;

    init_argv(nsite,nel,V,t,U,lambda,K,argc,argv);

    nel_up=(nel+1)/2;
    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif
    K=K*2.0*M_PI/nsite;
    nel_down=nel-nel_up;
    Sz=nel_up-nel_down;
    basis sector(nsite,nel_up,nel_down);
    sector.init();
    //sector.prlong();
    stringstream sf;
    string filename;
    sf<<"sector_"<<Sz;
    sf>>filename;

    lhamil config(sector,V,t,U,lambda,seed);
    config.coeff_explicit_update();
    config.diag();
    config.eigenstates_reconstruction();


    vector<double> E,A_hu,A_pu,A_hd,A_pd;
    double mu_hu,mu_hd,mu_pu,mu_pd;
    for(int i=0;i<3000;i++)
       E.push_back(i/200.0-5);

    Greens_func lspectral(config,sector);
    //lspectral.spectral_function_ii_uu_hole(0,0.2,E,Ah,mu_h);
    //lspectral.spectral_function_ii_uu_particle(0,0.2,E,Ap,mu_p);


    //spectral.spectral_function_ii_uu_hole_full_hamil(1,0.05,E,A_h,mu_h);
    //spectral.spectral_function_ii_uu_particle_full_hamil(1,0.05,E,A_p,mu_p);

    //lspectral.spectral_function_ij_uu_hole(3,5,0.05,E,Ah,mu_h);
    //lspectral.spectral_function_ij_uu_particle(3,5,0.05,E,Ap,mu_p);
    lspectral.spectral_function_kk_hole(K,0.05,E,A_hu,mu_hu,0);
    lspectral.spectral_function_kk_hole(K,0.05,E,A_hd,mu_hd,1);
    lspectral.spectral_function_kk_particle(K,0.05,E,A_pu,mu_pu,0);
    lspectral.spectral_function_kk_particle(K,0.05,E,A_pd,mu_pd,1);

    mu=(mu_hu+mu_pd+mu_hd+mu_pd)/4.0;
    for(int i=0;i<3000;i++)
       cout<<E[i]-mu<<" "<<-A_hu[i]-A_hd[i]+A_pu[i]+A_pd[i]<<endl;


      // cout<<E[i]-mu_h/2.0<<" "<<Ah[i]<<endl;
     // print dos

    /* basis check */

    //config.print_hamil();
    //config.print_lhamil(6);
    //config.print_eigen(6);
    // config.save_to_file(filename.c_str());

    /* ground-state energy check */
    //cout<<"#U  Lanczos Lieb-Wu "<<endl;
    //cout<<U<<" "<<config.ground_state_energy()/nsite<<" "<<-Integrate_Lieb_Wu(U)<<endl;

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
