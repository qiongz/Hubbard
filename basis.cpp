#include"basis.h"

using namespace std;

basis::basis() {}

basis::basis(int _uc,int _dc){

}

const basis & basis::operator =(const basis & _basis){
    if(this !=_basis){
        u_config=_basis.u_config;
        d_config=_basis.d_config;
        u_id=_basis.u_id;
        d_id=_basis.d_id;
        u_J=_basis.u_J;
        d_J=_basis.d_J;
        J=_basis.J;
    }
    return *this;
}

basis::~basis(){}

basis basis::hopping(int _i,int _j){
    
    
}

basis basis::potential(int _i){

}
