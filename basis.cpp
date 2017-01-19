#include"basis.h"

using namespace std;

basis::basis(){
}

basis::basis(int _n_site,int _n_elu, int _n_eld):n_site(_n_site),n_elu(_n_elu),n_eld(_n_eld){
    init();
}

const basis & basis::operator =(const basis & _basis){
    if(this !=&_basis){
        n_site=_basis.n_site;
        n_elu=_basis.n_elu;
        n_eld=_basis.n_eld;
        nb_up=_basis.nb_up;
        nb_down=_basis.nb_down;
    }
    return *this;
}

basis::~basis(){}

int basis::factorial(int N, int m){
    unsigned long num,denum;
    int i;
    num=1;
    for(i=N-m+1;i<=N;i++)
        num*=i;
    denum=1;
    for(i=1;i<=m;i++)
        denum*=i;
    return num/denum;
}

void basis::init(){
    int i,config_init;
    nb_up=factorial(n_site,n_elu);
    nb_down=factorial(n_site,n_eld);
    n_basis=nb_up*nb_down;
    config_init=0;
    for(i=0;i<n_elu;i++)
       config_init+=(1<<i);
    generate_up(config_init);

    config_init=0;
    for(i=0;i<n_eld;i++)
       config_init+=(1<<i);
    generate_down(config_init);
    
    for(auto &x:elu)
      relu.push_back(x.first);
    for(auto &x:eld)
      reld.push_back(x.first);

    sort(relu.begin(),relu.end());     
    sort(reld.begin(),reld.end());     
    elu.clear();
    eld.clear();
    for(i=0;i<nb_up;i++)
        elu[relu[i]]=i; 
    for(i=0;i<nb_down;i++)
        eld[reld[i]]=i; 
    
}


int basis::hopping_up(int i,int n){
   int mask,K,L,b;
   mask=(1<<n)+(1<<(n+1));
   K=mask&relu[i];
   L=K^mask;
   if(L!=0 && L!=mask)
       b=relu[i]-K+L;
   else
       b=relu[i];
   return elu[b];
}

int basis::hopping_down(int i,int n){
   int mask,K,L,b;
   mask=(1<<n)+(1<<(n+1));
   K=mask&reld[i];
   L=K^mask;
   if(L!=0 && L!=mask)
       b=reld[i]-K+L;
   else
       b=reld[i];
   return eld[b];
}

int basis::potential(int i,int j,int n){
   int mask,bu,bd;
   mask=(1<<n);
   bu=relu[i]&mask;
   bd=reld[i]&mask;
   if(bu==mask && bd==mask)
       return 1;
   else
       return 0;
}

void basis::generate_up(int a){
   elu.emplace(a,a);
   int mask,K,L,b;
   for(int i=0;i<n_site-1;i++){
       mask=(1<<i)+(1<<(i+1));
       K=mask&a;
       L=K^mask;
       if(L!=0 && L!=mask){
            b=a-K+L;
            if(elu.find(b)==elu.end())
                   generate_up(b);
            if(elu.size()==nb_up)
                 return;
       }
   }
   return;
}

void basis::generate_down(int a){
   eld.emplace(a,a);
   int mask,K,L,b;
   for(int i=0;i<n_site-1;i++){
       mask=(1<<i)+(1<<(i+1));
       K=mask&a;
       L=K^mask;
       if(L!=0 && L!=mask){
            b=a-K+L;
            if(eld.find(b)==eld.end())
                   generate_down(b);
            if(eld.size()==nb_down)
                 return;
       }
   }
}

void basis::print(){
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: elu)
         cout<<bitset<16>(x.first).to_string()<<" "<<setw(6)<<x.first<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(auto &x: eld)
         cout<<bitset<16>(x.first).to_string()<<" "<<setw(6)<<x.first<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis for spin-up electrons: "<<setw(6)<<nb_up<<endl;
    cout<<"No. basis for spin-down electrons: "<<setw(6)<<nb_down<<endl;
    cout<<"---------------------------------------"<<endl;
}


