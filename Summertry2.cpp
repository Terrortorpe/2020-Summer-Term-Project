#include <vector>
#include "Eigen/Dense"
#include <complex>
using namespace std;

typedef Eigen::MatrixXcd Matrix;
typedef Eigen::VectorXcd Vector;
complex<double> I = {0,1};



struct edge{
componentnode* A;
tyingnode* B;
double Current; //positive if going INTO the tyingnode
int CurrentID;
};

class node{

public:
   vector<edge*> edges;
    virtual void constraint(Matrix& mat){
        throw 1;
    }
};

class componentnode :public node{
    virtual void constraint(Matrix& mat){
        throw 1;
    }
};

class tyingnode :public node{

public:
    double Voltage;
    int VoltageID;
    

    // Currents have to sum to 0
    void constraint(Matrix& mat, Vector& vec){
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;
    }
};

class Capacitornode :public componentnode{
    double Capacitance;
    complex<double> Impedance;

    Capacitornode(double Capa,double omega){
        Capacitance = Capa;
        Impedance = 1.0/(I*Capacitance*omega);
    }

    void constraint(Matrix& mat, Vector& vec){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;


        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = 1;      // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistancce
        mat(mat.rows()-1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[0]->CurrentID) = Impedance;
        vec(vec.size()-1) = 0;


    }
};

class Inductornode :public componentnode{
    double Inductance;
    complex<double> Impedance;

    Inductornode(double Indu,double omega){
        Inductance = Indu;
        Impedance = I*Inductance*omega;
    }

    void constraint(Matrix& mat, Vector& vec){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;


        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[0]->CurrentID) = Impedance;
        vec(vec.size()-1) = 0;


    }
};

class resistornode :public componentnode{
    double Resistance;
    void constraint(Matrix& mat, Vector& vec){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;


        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[0]->CurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
        vec(vec.size()-1) = 0;


    }
};


class voltagesourcenode :public componentnode{
    double Voltage;
    void constraint(Matrix& mat, Vector& vec){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //Voltage source has larger voltage on the edge[1] side
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = 1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
        vec(vec.size()-1) = Voltage;                            //V(1)-V(0)=U


    }
};

class currentsourcenode :public componentnode{
    double Current;
    void constraint(Matrix& mat, Vector& vec){

        //Current source makes current flow from the 0 to the 1
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->CurrentID) = 1;
        vec(vec.size()-1) = -Current;

        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[1]->CurrentID) = 1;
        vec(vec.size()-1) = Current;
    }
};


// edge[0] is the Collector 
// edge[1] is the Base 
// edge[2] is the Emitter
/*class transistornode :public componentnode{
    float Beta;
    transistornode(float bet){
        Beta = bet;
    }
    void constraint(Matrix& mat, Vector& vec){

       //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[1]->CurrentID) = 1;
        vec(vec.size()-1) = Current;
    }
};*/

int main()
{
    return 0;
}