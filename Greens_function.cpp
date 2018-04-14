#include"Greens_function.h"

Greens_func::Greens_func(){}

Greens_func::Greens_func(lhamil & _gs_config,basis& _gs_sector){
    gs_config=_gs_config;
    gs_sector=_gs_sector;
    hole_sector.init(gs_sector.nsite,gs_sector.nel_up-1,gs_sector.nel_down);
    particle_sector.init(gs_sector.nsite,gs_sector.nel_up+1,gs_sector.nel_down);
    hole_config.init(hole_sector,gs_config.V,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
    particle_config.init(particle_sector,gs_config.V,gs_config.t,gs_config.U,gs_config.lambda,gs_config.seed);
}

Greens_func::~Greens_func(){}

void Greens_func::creation(long r,double coeff, int alpha)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(alpha==0)
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
    else
     for(j=0; j<gs_sector.nbasis_down; j++){
        // gs_sector has N electrons, particle_sector has N+1 electrons
        //apply operator c_r on basis id_up[i]
        Ob=gs_sector.creation(gs_sector.id_down[j],r);
        // if site r has no electron, then the returned value is different
        if(Ob!=gs_sector.id_down[j]){
           // find the corresponding index in the N+1 electron sector
           it=particle_sector.basis_down.find(Ob);
           if(it!=particle_sector.basis_down.end()){
             k=it->second;
             for(i=0;i<gs_sector.nbasis_up;i++)
                particle_phi_0[i*gs_sector.nbasis_down+k]+=coeff*gs_config.psir_0[i*gs_sector.nbasis_down+j];
           }
        }
    }
}

void Greens_func::annihilation(long r,double coeff, int alpha)
{
    long i,j,k,Ob;
    map<long,long>::iterator it;
    if(alpha==0)
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
    else
    for(j=0; j<gs_sector.nbasis_down; j++){
          // gs_sector has N electrons, hole_sector has N-1 electrons
          Ob=gs_sector.annihilation(gs_sector.id_down[j],r);
          if(Ob!=gs_sector.id_down[j]){
             it=hole_sector.basis_down.find(Ob);
             if(it!=hole_sector.basis_down.end()){
               k=it->second;
               for(i=0;i<gs_sector.nbasis_down;i++)
                  hole_phi_0[i*gs_sector.nbasis_down+k]+=coeff*gs_config.psir_0[i*gs_sector.nbasis_down+j];
             }
          }
    }
}

void Greens_func::spectral_function_ii_hole(long r_i, double eta, vector<double> &E, vector<double> &A,double &mu, int alpha){
    A.assign(E.size(),0);
    double norm=0;
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation(r_i,1,alpha);
    hole_config.coeff_update_wopt(hole_phi_0);
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.diag();
    mu=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
}

void Greens_func::spectral_function_ii_particle(long r_i, double eta, vector<double> &E, vector<double> &A,double &mu,int alpha){
    A.assign(E.size(),0);
    double norm=0;
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation(r_i,1,alpha);
    for(int i=0;i<particle_config.nHilbert;i++)
        norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    mu=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
}

void Greens_func::spectral_function_ij_hole(long r_i,long r_j, double eta,vector<double> &E,vector<double> &A,double &mu,int alpha)
{
    A.assign(E.size(),0);
    // O_ij=O_{i+j}-1/2(O_i+O_j)
    // operator O_{i+j}=1/sqrt{2}(c_i+c_j)
    double norm=0;
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation(r_i,1.0/sqrt(2),alpha);
    annihilation(r_j,1.0/sqrt(2),alpha);
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    mu=gs_config.E0-hole_config.eigenvalues[0];
    for(int i=0;i<E.size();i++)
        A[i]+=0.5*norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
    // calculate O_i=c_i
    hole_phi_0.assign(hole_config.nHilbert,0);
    annihilation(r_i,1.0/sqrt(2),alpha);
    annihilation(r_j,-1.0/sqrt(2),alpha);
    norm=0;
    for(int i=0;i<hole_config.nHilbert;i++)
        norm+=hole_phi_0[i]*hole_phi_0[i];
    hole_config.coeff_update_wopt(hole_phi_0);
    hole_config.diag();
    for(int i=0;i<E.size();i++)
        A[i]-=0.5*norm*hole_config.spectral_function(-E[i]+gs_config.E0,-eta);
}

void Greens_func::spectral_function_ij_particle(long r_i,long r_j, double eta,vector<double> &E,vector<double> &A,double &mu,int alpha)
{
    A.assign(E.size(),0);
    double norm=0;
    // O_ij=O_{i+j}-1/2(O_i+O_j)
    // operator O_{i+j}=1/sqrt{2}(c_i+c_j)
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation(r_i,1.0/sqrt(2),alpha);
    creation(r_j,1.0/sqrt(2),alpha);
    for(int i=0;i<particle_config.nHilbert;i++)
       norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    mu=particle_config.eigenvalues[0]-gs_config.E0;
    for(int i=0;i<E.size();i++)
        A[i]+=0.5*norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
    // operator O_{i-j}=1/sqrt{2}(c_i-c_j)
    particle_phi_0.assign(particle_config.nHilbert,0);
    creation(r_i,1.0/sqrt(2),alpha);
    creation(r_j,-1.0/sqrt(2),alpha);
    norm=0;
    for(int i=0;i<particle_config.nHilbert;i++)
       norm+=particle_phi_0[i]*particle_phi_0[i];
    particle_config.coeff_update_wopt(particle_phi_0);
    particle_config.diag();
    for(int i=0;i<E.size();i++)
        A[i]-=0.5*norm*particle_config.spectral_function(E[i]+gs_config.E0,eta);
}

void Greens_func::spectral_function_kk_hole(double k, double eta,vector<double> & E, vector<double> &A, double &mu, int alpha){
    vector<double> Ar;
    A.assign(E.size(),0);
    for(int r=0;r<gs_sector.nsite;r++){
        if(r==0)
           spectral_function_ii_hole(0,eta,E,Ar,mu,alpha);
        else
           spectral_function_ij_hole(0,r,eta,E,Ar,mu,alpha);
        for(int i=0;i<E.size();i++)
            A[i]+=cos(k*r)*Ar[i];
    }
}

void Greens_func::spectral_function_kk_particle(double k, double eta,vector<double> & E, vector<double> &A, double &mu,int alpha){
    vector<double> Ar;
    A.assign(E.size(),0);
    for(int r=0;r<gs_sector.nsite;r++){
        if(r==0)
           spectral_function_ii_particle(0,eta,E,Ar,mu,alpha);
        else
           spectral_function_ij_particle(0,r,eta,E,Ar,mu,alpha);
        for(int i=0;i<E.size();i++)
            A[i]+=cos(k*r)*Ar[i];
    }
}
