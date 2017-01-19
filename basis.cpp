#include"basis.h"

using namespace std;

basis::basis(){
}

basis::basis(int _n_site,int _n_elu, int _n_eld):n_site(_n_site),n_elu(_n_elu),n_eld(_n_eld){
    int i,config_init;
    nb_up=factorial(n_site,n_elu);
    nb_down=factorial(n_site,n_eld);
    config_init=0;
    for(i=0;i<n_elu;i++)
       config_init+=(1<<i);
    generate_up(config_init);

    config_init=0;
    for(i=0;i<n_eld;i++)
       config_init+=(1<<i);
    generate_down(config_init);
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

//basis basis::hopping(int i){

//}

//basis basis::potential(int _i){

//}

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
    cout<<"--------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"--------------------------------"<<endl;
    for(auto &x: elu)
         cout<<bitset<16>(x.first).to_string()<<" "<<setw(6)<<x.first<<endl;
    cout<<"--------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"--------------------------------"<<endl;
    for(auto &x: eld)
         cout<<bitset<16>(x.first).to_string()<<" "<<setw(6)<<x.first<<endl;
    cout<<"--------------------------------"<<endl;
}


