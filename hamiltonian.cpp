#include"hamiltonian.h"

void diag_hamil(basis sector,double t, double U, double energy, double *wf){
   int n_site,nb_up,nb_down,n_basis;
   n_site=sector.n_site;
   nb_up=sector.nb_up; 
   nb_down=sector.nb_down; 
   n_basis=nb_up*nb_down;

   // only the nonzero hamiltonian elements with their indices are stored
   std::map<int,double> hamil_nonzero;
   int n,i,j,k,l;
   for(i=0;i<nb_up;i++){
       for(j=0;j<nb_down;j++){
           for(n=0;n<n_site-1;n++){
               k=sector.hopping_up(i,n);
               l=sector.hopping_down(j,n);
               if(k!=i)
                  hamil_nonzero[(i*nb_down+j)*n_basis+k*nb_down+j]+=-t;
               if(l!=j)
                  hamil_nonzero[(i*nb_down+j)*n_basis+i*nb_down+l]+=-t;
           }
           for(n=0;n<n_site;n++)
               if(sector.potential(i,j,n))
                  hamil_nonzero[(i*nb_down+j)*n_basis+i*nb_down+j]+=U;
       }
   }
   /* printing the hamiltonian matrix */
   /*
   std::cout<<"hamiltonain nonzero elements:"<<std::endl;
   for(auto &x:hamil_nonzero)
       std::cout<<x.first/n_basis<<" "<<x.first%n_basis<<" "<<x.second<<std::endl; 
   */

}

void hoperation(){



}
