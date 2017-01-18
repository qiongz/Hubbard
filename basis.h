#include<iostream>
using namespace std;

class basis{
private:
    int u_config,d_config;   // up/down-spin configurations
    int u_id,d_id;           // up/down-spin identification numbers
    int u_J,d_J,J;             // up/down-spin J values, total J value

public:
    explicit basis();
    basis(int,int);
    const basis & operator=(const basis &);
    ~basis();
    basis hopping(int,int);
    basis potential(int);
};
    
