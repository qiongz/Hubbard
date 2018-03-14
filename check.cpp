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

    nsite=8;
    nel=8;
    t=1;
    U=5;
    K=1.0;
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

    Greens_func spectral(config);
    vector<double> E,A;
    for(int i=0;i<2000;i++)
       E.push_back(i/100.0-10);
    //spectral.spectral_function_ij_uu(1,1,0.1,E,A,mu);
    spectral.spectral_function_kk_uu(K,0.2,E,A,mu);
    for(int i=0;i<2000;i++)
       cout<<E[i]-mu<<" "<<A[i]<<endl;

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
