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
    double t,U,K,mu;
    double func_Lieb_Wu(double x, void *params);
    double Integrate_Lieb_Wu(double U);

    nsite=10;
    nel=10;
    t=1;
    U=5;
    K=0.5;
    lambda=300;

    init_argv(nsite,nel,t,U,lambda,K,argc,argv);

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

    vector<double> E,A_h,A_p,Ah,Ap;
    double mu_h,mu_p;
    for(int i=0;i<4000;i++)
       E.push_back(i/200.0-10);

    Greens_func lspectral(config);
    lspectral.spectral_function_ii_uu_hole(0,0.1,E,Ah,mu_h);
    lspectral.spectral_function_ii_uu_particle(0,0.1,E,Ap,mu_p);

    //spectral.spectral_function_ii_uu_hole_full_hamil(1,0.05,E,A_h,mu_h);
    //spectral.spectral_function_ii_uu_particle_full_hamil(1,0.05,E,A_p,mu_p);

    //lspectral.spectral_function_ij_uu_hole(3,5,0.05,E,Ah,mu_h);
    //lspectral.spectral_function_ij_uu_particle(3,5,0.05,E,Ap,mu_p);
    //lspectral.spectral_function_kk_uu_hole(K,0.01,E,Ah,mu_h);
    //lspectral.spectral_function_kk_uu_particle(K,0.01,E,Ap,mu_p);

    mu=(mu_h+mu_p)/2.0;
    for(int i=0;i<4000;i++)
       cout<<E[i]-mu<<" "<<-Ah[i]+Ap[i]<<endl;

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
