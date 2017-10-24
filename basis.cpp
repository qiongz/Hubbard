#include"basis.h"
using namespace std;

basis::basis(){
}

basis::basis(long _nsite,long _nel_up, long _nel_down):nsite(_nsite),nel_up(_nel_up),nel_down(_nel_down){
}

const basis & basis::operator =(const basis & _basis){
    if(this !=&_basis){
        nsite=_basis.nsite;
        nel_up=_basis.nel_up;
        nel_down=_basis.nel_down;
        nbasis_up=_basis.nbasis_up;
        nbasis_down=_basis.nbasis_down;
        basis_up=_basis.basis_up;
        basis_down=_basis.basis_down;
        id_up=_basis.id_up;
        id_down=_basis.id_down;
    }
    return *this;
}

basis::~basis(){}

long basis::factorial(long N, long m){
    unsigned long num,denum;
    long i;
    num=1;
    for(i=N-m+1;i<=N;i++)
        num*=i;
    denum=1;
    for(i=1;i<=m;i++)
        denum*=i;
    return num/denum;
}

void basis::init(){
    long i,config_init;
    nbasis_up=factorial(nsite,nel_up);
    nbasis_down=factorial(nsite,nel_down);
    config_init=0;
    for(i=0;i<nel_up;i++)
       config_init+=(1<<i);
    generate_up(config_init);

    config_init=0;
    for(i=0;i<nel_down;i++)
       config_init+=(1<<i);
    generate_down(config_init);
    
    for(auto &x:basis_up)
      id_up.push_back(x.first);
    for(auto &x:basis_down)
      id_down.push_back(x.first);

    sort(id_up.begin(),id_up.end());     
    sort(id_down.begin(),id_down.end());     
    basis_up.clear();
    basis_down.clear();
    for(i=0;i<nbasis_up;i++)
        basis_up[id_up[i]]=i; 
    for(i=0;i<nbasis_down;i++)
        basis_down[id_down[i]]=i; 
}

long basis::hopping_up(long i,long n,long m){
   long mask,K,L,b;
   if(m<0) m+=nsite;
    else if (m>=nsite) m-=nsite;

   mask=(1<<n)+(1<<m);
   K=mask&id_up[i];
   L=K^mask;
   if(L!=0 && L!=mask)
       b=id_up[i]-K+L;
   else
       b=id_up[i];
   return basis_up[b];
}

long basis::hopping_down(long i,long n,long m){
   long mask,K,L,b;
   if(m<0) m+=nsite;
    else if (m>=nsite) m-=nsite;
   mask=(1<<n)+(1<<m);
   K=mask&id_down[i];
   L=K^mask;
   if(L!=0 && L!=mask)
       b=id_down[i]-K+L;
   else
       b=id_down[i];
   return basis_down[b];
}

long basis::potential(long i,long j,long n){
   long mask,bu,bd;
   mask=(1<<n);
   bu=id_up[i]&mask;
   bd=id_down[j]&mask;
   if(bu==mask && bd==mask)
       return 1;
   else
       return 0;
}

// c_i^\dagger c_i terms for spin-up
long basis::onsite_up(long i,long n){
  long mask,bu;
  mask=(1<<n);
  bu=id_up[i]&mask;
  if(bu==mask)
    return 1;
  else
    return 0; 
}

// c_i^\dagger c_i terms for spin-down
long basis::onsite_down(long i,long n){
  long mask,bd;
  mask=(1<<n);
  bd=id_down[i]&mask;
  if(bd==mask)
    return 1;
  else
    return 0; 
}


void basis::creation_el_up(long n)
{
    long mask,bu;
    mask=(1<<n);
    std::map<long,long> basis_up_particle;
    int count=0;
    for( auto &x:basis_up){
      bu=x.first&mask;
      if(bu!=mask) basis_up_particle.insert(pair<long,long>(x.first+mask,count++));
    } 

    id_up.clear();
    for(auto &x:basis_up_particle)
      id_up.push_back(x.first);

    sort(id_up.begin(),id_up.end());     

    basis_up=basis_up_particle;
    basis_up_particle.clear(); 
    nbasis_up=basis_up.size();
}

void basis::creation_el_down(long n){
    long mask,bd;
    mask=(1<<n);
    std::map<long,long> basis_down_particle;
    int count=0;
    for( auto &x:basis_down){
      bd=x.first&mask;
      if(bd!=mask) basis_down_particle.insert(pair<long,long>(x.first+mask,count++));
    } 
    id_down.clear();
    for(auto &x:basis_down_particle)
      id_down.push_back(x.first);
    sort(id_down.begin(),id_down.end());     
    basis_down=basis_down_particle;
    basis_down_particle.clear(); 
    nbasis_down=basis_down.size();
}

void basis::annihilation_el_up(long n)
{
    long mask,bu;
    mask=(1<<n);
    std::map<long,long> basis_up_hole;
    int count=0;
    for( auto &x:basis_up){
      bu=x.first&mask;
      if(bu==mask) basis_up_hole.insert(pair<long,long>(x.first-mask,count++));
    } 
    id_up.clear();
    for(auto &x:basis_up_hole)
      id_up.push_back(x.first);
    sort(id_up.begin(),id_up.end());     
    basis_up=basis_up_hole;
    basis_up_hole.clear(); 
    nbasis_up=basis_up.size();
}

void basis::annihilation_el_down(long n)
{
    long mask,bd;
    mask=(1<<n);
    std::map<long,long> basis_down_hole;
    int count=0;
    for( auto &x:basis_down){
      bd=x.first&mask;
      if(bd==mask) basis_down_hole.insert(pair<long,long>(x.first-mask,count++));
    } 
    id_down.clear();
    for(auto &x:basis_down_hole)
      id_down.push_back(x.first);
    sort(id_down.begin(),id_down.end());     
    basis_down=basis_down_hole;
    basis_down_hole.clear(); 
    nbasis_down=basis_down.size();
}

void basis::generate_up(long a){
   long mask,K,L,b,j;
   basis_up.emplace(a,a);
   for(long i=0;i<nsite;i++){
       j=(i+1>=nsite)?i+1-nsite:i+1; 
       mask=(1<<i)+(1<<j);
       K=mask&a;
       L=K^mask;
       if(L!=0 && L!=mask){
            b=a-K+L;
            if(basis_up.find(b)==basis_up.end())
                   generate_up(b);
            if(basis_up.size()==nbasis_up)
                 return;
       }
   }
   return;
}

void basis::generate_down(long a){
   long mask,K,L,b,j;
   basis_down.emplace(a,a);
   for(long i=0;i<nsite-1;i++){
       j=(i+1>=nsite)?i+1-nsite:i+1; 
       mask=(1<<i)+(1<<j);
       K=mask&a;
       L=K^mask;
       if(L!=0 && L!=mask){
            b=a-K+L;
            if(basis_down.find(b)==basis_down.end())
                   generate_down(b);
            if(basis_down.size()==nbasis_down)
                 return;
       }
   }
}

void basis::prlong(){
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: basis_up)
         cout<<bitset<20>(x.first).to_string()<<" "<<setw(6)<<x.first<<" "<<x.second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: basis_down)
         cout<<bitset<20>(x.first).to_string()<<" "<<setw(6)<<x.first<<" "<<x.second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis for spin-up electrons: "<<setw(6)<<nbasis_up<<endl;
    cout<<"No. basis for spin-down electrons: "<<setw(6)<<nbasis_down<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"Lin's Table:"<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: id_up)
       cout<<x<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: id_down)
       cout<<x<<endl;
    cout<<"---------------------------------------"<<endl;
}


