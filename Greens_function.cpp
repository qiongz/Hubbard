#include"Greens_function.h"

Greens_func::Greens_func(){}

Greens_func::Greens_func(lhamil & _gs_config){
    gs_config=_gs_config;
    gs_sector=gs_config.sector;
    hole_sector.init(gs_sector.nsite,gs_sector.nel_up-1,gs_sector.nel_down);
    particle_sector.init(gs_sector.nsite,gs_sector.nel_up+1,gs_sector.nel_down);
    hole_config.init(hole_sector,gs_config.V,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
    particle_config.init(particle_sector,gs_config.V,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
}

Greens_func::Greens_func(hamil & _gs_hconfig,basis &_gs_sector){
    gs_hconfig=_gs_hconfig;
    gs_sector=_gs_sector;
    hole_sector.init(gs_sector.nsite,gs_sector.nel_up-1,gs_sector.nel_down);
    particle_sector.init(gs_sector.nsite,gs_sector.nel_up+1,gs_sector.nel_down);
    hole_hconfig.init(hole_sector,gs_hconfig.V,gs_hconfig.t,gs_hconfig.U);
    particle_hconfig.init(particle_sector,gs_hconfig.V,gs_hconfig.t,gs_hconfig.U);
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

void Greens_func::creation_u_full_hamil(long r,double coeff)
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
                particle_phi_0[k*gs_sector.nbasis_down+j]+=coeff*gs_hconfig.psi_0[i*gs_sector.nbasis_down+j];
           }
        }
    }
}

void Greens_func::annihilation_u_full_hamil(long r,double coeff)
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
                  hole_phi_0[k*gs_sector.nbasis_down+j]+=coeff*gs_hconfig.psi_0[i*gs_sector.nbasis_down+j];
             }
          }
    }
}

void Greens_func::spectral_function_ii_uu_hole(long r_i, double eta, vector<double> &E, vector<double> &A,double &mu){
    A.assign(E.size(),0);
    double norm=0;
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation_u(r_i,1);
    hole_config.coeff_update_wopt(hole_phi_0);
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.diag();
    mu=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
}

void Greens_func::spectral_function_ii_uu_particle(long r_i, double eta, vector<double> &E, vector<double> &A,double &mu){
    A.assign(E.size(),0);
    double norm=0;
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation_u(r_i,1);
    for(int i=0;i<particle_config.nHilbert;i++)
        norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    mu=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
}

void Greens_func::spectral_function_ij_uu_hole(long r_i,long r_j, double eta,vector<double> &E,vector<double> &A,double &mu)
{
    A.assign(E.size(),0);
    // O_ij=O_{i+j}-1/2(O_i+O_j)
    // operator O_{i+j}=1/sqrt{2}(c_i+c_j)
    double norm=0;
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation_u(r_i,1.0/sqrt(2));
    annihilation_u(r_j,1.0/sqrt(2));
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    mu=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=0.5*norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
    // calculate O_i=c_i
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation_u(r_i,1.0/sqrt(2));
    annihilation_u(r_j,-1.0/sqrt(2));
    norm=0;
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    for(int i=0;i<E.size();i++)
        A[i]-=0.5*norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
}

void Greens_func::spectral_function_ij_uu_particle(long r_i,long r_j, double eta,vector<double> &E,vector<double> &A,double &mu)
{
    A.assign(E.size(),0);
    double norm=0;
    // O_ij=O_{i+j}-1/2(O_i+O_j)
    // operator O_{i+j}=1/sqrt{2}(c_i+c_j)
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation_u(r_i,1.0/sqrt(2));
    creation_u(r_j,1.0/sqrt(2));
    for(int i=0;i<particle_config.nHilbert;i++)
       norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    mu=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=0.5*norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
    // operator O_{i-j}=1/sqrt{2}(c_i-c_j)
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation_u(r_i,1.0/sqrt(2));
    creation_u(r_j,-1.0/sqrt(2));
    norm=0;
    for(int i=0;i<particle_config.nHilbert;i++)
       norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    for(int i=0;i<E.size();i++)
        A[i]-=0.5*norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
}


/*
void Greens_func::spectral_function_kk_uu_hole(double k, double eta,vector<double> & E, vector<double> &A,double &mu){
    A.assign(E.size(),0);
    // initialize the operator applied wave function |O_phi_0>
    // C_k perform Fourier-transform and sum over r-space operators
    double norm=0;
    hole_phi_0.assign(hole_config.nHilbert,0);
    for(int r=0;r<gs_sector.nsite;r++)
      annihilation_u(r,cos(r*k)/sqrt(gs_sector.nsite));
    for(int i=0;i<particle_config.nHilbert;i++)
        norm+=particle_phi_0[i]*particle_phi_0[i];
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    mu=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=norm*hole_config.spectral_function(E[i],eta);
}

void Greens_func::spectral_function_kk_uu_particle(double k, double eta,vector<double> & E, vector<double> &A,double &mu){
    A.assign(E.size(),0);
    // initialize the operator applied wave function |O_phi_0>
    // C_k perform Fourier-transform and sum over r-space operators
    particle_phi_0.assign(particle_config.nHilbert,0);
    for(int r=0;r<gs_sector.nsite;r++)
      creation_u(r,cos(r*k)/sqrt(gs_sector.nsite));
    double norm=0;
    for(int i=0;i<particle_config.nHilbert;i++)
        norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    mu=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=norm*particle_config.spectral_function(E[i],eta);
}
*/

void Greens_func::spectral_function_kk_uu_hole(double k, double eta,vector<double> & E, vector<double> &A, double &mu){
    vector<double> Ar;
    A.assign(E.size(),0);
    for(int r=0;r<gs_sector.nsite;r++){
        if(r==0)
           spectral_function_ii_uu_hole(0,eta,E,Ar,mu);
        else
           spectral_function_ij_uu_hole(0,r,eta,E,Ar,mu);
        for(int i=0;i<E.size();i++)
            A[i]+=cos(k*r)*Ar[i];
    }
}

void Greens_func::spectral_function_kk_uu_particle(double k, double eta,vector<double> & E, vector<double> &A, double &mu){
    vector<double> Ar;
    A.assign(E.size(),0);
    for(int r=0;r<gs_sector.nsite;r++){
        if(r==0)
           spectral_function_ii_uu_particle(0,eta,E,Ar,mu);
        else
           spectral_function_ij_uu_particle(0,r,eta,E,Ar,mu);
        for(int i=0;i<E.size();i++)
            A[i]+=cos(k*r)*Ar[i];
    }
}
