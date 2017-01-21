#include"hamiltonian.h"

void diag_hamil(basis *sector,double t, double U, double energy, double *wf){
   int nsite,nbasis_up,nbasis_down,nbasis;
   nsite=(*sector).nsite;
   nbasis_up=(*sector).nbasis_up; 
   nbasis_down=(*sector).nbasis_down; 
   nbasis=nbasis_up*nbasis_down;

   // only the nonzero hamiltonian elements with their indices are stored
   std::map<int,double> hamil_nonzero;
   int n,i,j,k,l;
   for(i=0;i<nbasis_up;i++){
       for(j=0;j<nbasis_down;j++){
           for(n=0;n<nsite-1;n++){
               k=(*sector).hopping_up(i,n);
               l=(*sector).hopping_down(j,n);
               if(k!=i)
                  hamil_nonzero[(i*nbasis_down+j)*nbasis+k*nbasis_down+j]+=-t;
               if(l!=j)
                  hamil_nonzero[(i*nbasis_down+j)*nbasis+i*nbasis_down+l]+=-t;
           }
           for(n=0;n<nsite;n++)
               if((*sector).potential(i,j,n))
                  hamil_nonzero[(i*nbasis_down+j)*nbasis+i*nbasis_down+j]+=U;
       }
   }
   /* printing the hamiltonian matrix */
   /*
   std::cout<<"hamiltonain nonzero elements:"<<std::endl;
   for(auto &x:hamil_nonzero)
       std::cout<<x.first/n_basis<<" "<<x.first%n_basis<<" "<<x.second<<std::endl; 
   */

}

//void hoperation(){
//}
