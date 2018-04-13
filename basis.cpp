#include"basis.h"
using namespace std;

basis::basis() {
}

basis::basis(long _nsite,long _nel_up, long _nel_down):nsite(_nsite),nel_up(_nel_up),nel_down(_nel_down) {
}

const basis & basis::operator =(const basis & _basis) {
    if(this !=&_basis) {
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

basis::~basis() {}

long basis::factorial(long N, long m) {
    unsigned long num,denum;
    long i;
    num=1;
    for(i=N-m+1; i<=N; i++)
        num*=i;
    denum=1;
    for(i=1; i<=m; i++)
        denum*=i;
    return num/denum;
}

void basis::init() {
    long i,j,count,config;
    std::map<long,long>::iterator it;
    count=j=config=0;
    generate_up(count,j,config);
    count=j=config=0;
    generate_down(count,j,config);

    for(it=basis_up.begin(); it!=basis_up.end(); it++)
        id_up.push_back(it->first);
    for(it=basis_down.begin(); it!=basis_down.end(); it++)
        id_down.push_back(it->first);

    sort(id_up.begin(),id_up.end());
    sort(id_down.begin(),id_down.end());
    basis_up.clear();
    basis_down.clear();
    for(i=0; i<id_up.size(); i++)
        basis_up[id_up[i]]=i;
    for(i=0; i<id_down.size(); i++)
        basis_down[id_down[i]]=i;
}

void basis::init(long _nsite, long _nel_up, long _nel_down){
   nsite=_nsite;
   nel_up=_nel_up;
   nel_down=_nel_down;
   init();
}

long basis::hopping_up(long i,long n,long m) {
    long mask,K,L,b;
    if(m<0) m+=nsite;
    else if (m>=nsite) m-=nsite;

    mask=(1<<n)+(1<<m);
    K=mask&id_up[i];
    L=K^mask;

    if(L!=0 && L!=mask) {
        b=id_up[i]-K+L;
        if(basis_up.find(b)!=basis_up.end())
            return basis_up[b];
        else
            return basis_up[id_up[i]];
    }
    else
        return basis_up[id_up[i]];
}

long basis::hopping_down(long i,long n,long m) {
    long mask,K,L,b;
    if(m<0) m+=nsite;
    else if (m>=nsite) m-=nsite;

    mask=(1<<n)+(1<<m);
    K=mask&id_down[i];
    L=K^mask;
    if(L!=0 && L!=mask) {
        b=id_down[i]-K+L;
        if(basis_down.find(b)!=basis_down.end())
            return basis_down[b];
        else
            return basis_down[id_down[i]];
    }
    else
        return basis_down[id_down[i]];
}

long basis::potential(long i,long j,long n) {
    long mask,bu,bd;
    mask=(1<<n);
    bu=id_up[i]&mask;
    bd=id_down[j]&mask;
    if(bu==mask && bd==mask)
        return 1;
    else
        return 0;
}

long basis::onsite_up(long i,long n) {
    long mask,bu;
    mask=(1<<n);
    bu=id_up[i]&mask;
    if(bu==mask)
        return 1;
    else
        return 0;
}

long basis::onsite_down(long i,long n) {
    long mask,bd;
    mask=(1<<n);
    bd=id_down[i]&mask;
    if(bd==mask)
        return 1;
    else
        return 0;
}


long basis::creation(long s,long n)
{
    long mask,bu;
    mask=(1<<n);
    bu=s&mask;
    // there's no electron on site n
    // which means s_n=0
    if(bu!=mask)
        return s+mask;
    // there's already electron on site n
    else
        return s;
}

long basis::annihilation(long s,long n)
{
    long mask,bu;
    mask=(1<<n);
    bu=s&mask;
    if(bu==mask)
        return s-mask;
    else
        return s;
}


void basis::generate_up(long count,long j, long config){
  long i,id,k,c;
  count++;
  config=(count==1?0:config);
  j=(count==1?0:j+1);
  for(i=j;i<nsite-(nel_up-count);i++){
    c=config+(1<<i);
    if(count<nel_up)
       generate_up(count,i,c);
    else
        basis_up[c]=c;
  }
}

void basis::generate_down(long count,long j, long config) {
  long i,id,k,c;
  count++;
  config=(count==1?0:config);
  j=(count==1?0:j+1);
  for(i=j;i<nsite-(nel_down-count);i++){
    c=config+(1<<i);
    if(count<nel_down)
       generate_down(count,i,c);
    else
        basis_down[c]=c;
  }
}

void basis::prlong() {
    std::map<long,long>::iterator it;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-up electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_up.begin(); it!=basis_up.end(); it++)
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"spin-down electrons:"<<endl;
    cout<<"---------------------------------------"<<endl;
    for(it=basis_down.begin(); it!=basis_down.end(); it++)
        cout<<bitset<20>(it->first).to_string()<<" "<<setw(6)<<it->first<<" "<<it->second<<endl;
    cout<<"---------------------------------------"<<endl;
    cout<<"No. basis for spin-up electrons: "<<setw(6)<<nbasis_up<<endl;
    cout<<"No. basis for spin-down electrons: "<<setw(6)<<nbasis_down<<endl;
    /*
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
    */
}
