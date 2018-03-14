#include"Greens_function.h"

Greens_func::Greens_func(){}

Greens_func::Greens_func(lhamil & _gs_config){
    gs_config=_gs_config;
    gs_sector=gs_config.sector;
    hole_sector.init(gs_sector.nsite,gs_sector.nel_up-1,gs_sector.nel_down);
    particle_sector.init(gs_sector.nsite,gs_sector.nel_up+1,gs_sector.nel_down);
    hole_config.init(hole_sector,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
    particle_config.init(particle_sector,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
}

Greens_func::~Greens_func(){}

void Greens_func::creation_u(long r,double coeff)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    for(i=0; i<gs_sector.nbasis_up; i++){
        // gs_sector has N electrons, particle_sector has N+1 electrons
        //apply operator c_r on basis id_up[i]
        Ob=gs_sector.creation(gs_sector.id_up[i],r);
        // if site r has no electron, then the returned value is different
        if(Ob!=gs_sector.id_up[i]){
           // find the corresponding index in the N+1 electron sector
           it=particle_sector.basis_up.find(Ob);
           if(it!=particle_sector.basis_up.end()){
             k=it->second;
             for(j=0;j<gs_sector.nbasis_down;j++)
                particle_phi_0[k*gs_sector.nbasis_down+j]+=coeff*gs_config.psir_0[i*gs_sector.nbasis_down+j];
           }
        }
    }
}

void Greens_func::annihilation_u(long r,double coeff)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    for(i=0; i<gs_sector.nbasis_up; i++){
          // gs_sector has N electrons, hole_sector has N-1 electrons
          Ob=gs_sector.annihilation(gs_sector.id_up[i],r);
          if(Ob!=gs_sector.id_up[i]){
             it=hole_sector.basis_up.find(Ob);
             if(it!=hole_sector.basis_up.end()){
               k=it->second;
               for(j=0;j<gs_sector.nbasis_down;j++)
                  hole_phi_0[k*gs_sector.nbasis_down+j]+=coeff*gs_config.psir_0[i*gs_sector.nbasis_down+j];
             }
          }
    }
}

void Greens_func::spectral_function_ii_uu(int r_i,double eta,vector<double> &E, vector<double> &A,double &mu){
    A.assign(E.size(),0);
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation_u(r_i,1);
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    double mu_h=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=-hole_config.Greens_function(E[i],gs_config.E0,eta,1).imag()/M_PI;
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation_u(r_i,1);
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    double mu_p=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=-particle_config.Greens_function(E[i],gs_config.E0,eta,0).imag()/M_PI;
    mu=(mu_h+mu_p)/2.0;
}

void Greens_func::spectral_function_ij_uu(int r_i,long r_j, double eta,vector<double> &E,vector<double> &A,double &mu)
{
    A.assign(E.size(),0);
    // O_ij=O_{i+j}-1/2(O_i+O_j)
    // operator O_{i+j}=1/sqrt{2}(c_i+c_j)
    hole_phi_0.assign(hole_config.nHilbert,0);
    particle_phi_0.assign(particle_config.nHilbert,0);
    annihilation_u(r_i,1.0/sqrt(2));
    annihilation_u(r_j,1.0/sqrt(2));
    creation_u(r_i,1.0/sqrt(2));
    creation_u(r_j,1.0/sqrt(2));
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    double mu_h=gs_config.E0-hole_config.eigenvalues[0];
    double mu_p=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++){
        A[i]+=-hole_config.Greens_function(E[i],gs_config.E0,eta,1).imag()/M_PI;
        A[i]+=-particle_config.Greens_function(E[i],gs_config.E0,eta,0).imag()/M_PI;
    }
    // calculate O_i=c_i
    vector<double> A_ii,A_jj;
    spectral_function_ii_uu(r_i,eta,E,A_ii,mu);
    spectral_function_ii_uu(r_j,eta,E,A_jj,mu);
    for(int i=0;i<E.size();i++)
        A[i]-=0.5*(A_ii[i]+A_jj[i]);

    mu=(mu_h+mu_p)/2.0;
}


void Greens_func::spectral_function_kk_uu(double k, double eta,vector<double> & E, vector<double> &A,double &mu){
    A.assign(E.size(),0);

    // initialize the operator applied wave function |O_phi_0>
    // C_k perform Fourier-transform and sum over r-space operators
    particle_phi_0.assign(particle_config.nHilbert,0);
    hole_phi_0.assign(hole_config.nHilbert,0);
    for(int r=0;r<gs_sector.nsite;r++){
      creation_u(r,cos(r*k)/sqrt(gs_sector.nsite));
      annihilation_u(r,cos(r*k)/sqrt(gs_sector.nsite));
    }
    particle_config.coeff_update_wopt(particle_phi_0);
    hole_config.coeff_update_wopt(hole_phi_0);
    particle_config.diag();
    hole_config.diag();
    double mu_p=particle_config.eigenvalues[0]-gs_config.E0;
    double mu_h=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++){
        A[i]+=-particle_config.Greens_function(E[i],gs_config.E0,eta,0).imag()/M_PI;
        A[i]+=-hole_config.Greens_function(E[i],gs_config.E0,eta,1).imag()/M_PI;
    }
    mu=(mu_h+mu_p)/2.0;
}

/*
void Greens_func::spectral_function_kk_uu(double k, double eta,vector<double> & E, vector<double> &A, double &mu){
    vector<double> Ar;
    A.assign(E.size(),0);
    for(int r=0;r<gs_sector.nsite;r++){
        spectral_function_ij_uu(0,r,eta,E,Ar,mu);
        for(int i=0;i<E.size();i++)
            A[i]+=cos(k*r)/sqrt(gs_sector.nsite)*Ar[i];
    }
}
*/
