﻿#include <vector>
#include "Eigen/Dense"
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

typedef Eigen::MatrixXcd Matrix;
typedef Eigen::VectorXcd Vector;
complex<double> I = { 0,1 };
map<string, tyingnode*> PairT;
map<string, componentnode*> PairC;
int IDcount;
double omega;



struct edge {
    componentnode* A;
    tyingnode* B;
    double Current; //positive if going INTO the tyingnode
    int CurrentID;

    edge(componentnode* a_, tyingnode* b_) {
        A = a_;
        B = b_;
        CurrentID = IDcount;
        IDcount++;
        A->edges.push_back(this);
        B->edges.push_back(this);
    }
};

class node {

public:
    vector<edge*> edges;
    virtual void constraint(Matrix& mat) {
        throw 1;
    }
};

class componentnode :public node {
    virtual void constraint(Matrix& mat) {
        throw 1;
    }
};

class tyingnode :public node {

public:
    double Voltage;
    int VoltageID;


    // Currents have to sum to 0
    void constraint(Matrix& mat, Vector& vec) {
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        for (auto&& edge : edges) {
            mat(mat.rows() - 1, edge->CurrentID) = 1;
        }
        vec(vec.size() - 1) = 0;
    }
};

class Capacitornode :public componentnode {
public:
    double Capacitance;
    complex<double> Impedance;

    Capacitornode(stringstream& ss) {

    }

    void constraint(Matrix& mat, Vector& vec) {

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        for (auto&& edge : edges) {
            mat(mat.rows() - 1, edge->CurrentID) = 1;
        }
        vec(vec.size() - 1) = 0;


        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[0]->B->VoltageID) = 1;      // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistancce
        mat(mat.rows() - 1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows() - 1, edges[0]->CurrentID) = Impedance;
        vec(vec.size() - 1) = 0;


    }
};

class Inductornode :public componentnode {
public:
    double Inductance;
    complex<double> Impedance;

    Inductornode(stringstream& ss) {

    }

    void constraint(Matrix& mat, Vector& vec) {

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        for (auto&& edge : edges) {
            mat(mat.rows() - 1, edge->CurrentID) = 1;
        }
        vec(vec.size() - 1) = 0;


        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows() - 1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows() - 1, edges[0]->CurrentID) = Impedance;
        vec(vec.size() - 1) = 0;


    }
};

class resistornode :public componentnode {
public:
    double Resistance;
    resistornode(stringstream& ss) {

    };

    void setResistance(double sRes)
    {
        Resistance = sRes;
    };

    void constraint(Matrix& mat, Vector& vec) {

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        for (auto&& edge : edges) {
            mat(mat.rows() - 1, edge->CurrentID) = 1;
        }
        vec(vec.size() - 1) = 0;


        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows() - 1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows() - 1, edges[0]->CurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
        vec(vec.size() - 1) = 0;


    }
};


class voltagesourcenode :public componentnode {
public:
    double Voltage;

    voltagesourcenode(stringstream& ss) {

    }

    void setVoltage(double pVoltage)
    {
        Voltage = pVoltage;
    }

    void constraint(Matrix& mat, Vector& vec) {

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        for (auto&& edge : edges) {
            mat(mat.rows() - 1, edge->CurrentID) = 1;
        }
        vec(vec.size() - 1) = 0;

        //Voltage source has larger voltage on the edge[1] side
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[0]->B->VoltageID) = -1;
        mat(mat.rows() - 1, edges[1]->B->VoltageID) = 1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
        vec(vec.size() - 1) = Voltage;                            //V(1)-V(0)=U


    }
};

class currentsourcenode :public componentnode {
public:
    double Current;

    currentsourcenode(stringstream& ss) {

    }
    void constraint(Matrix& mat, Vector& vec) {

        //Current source makes current flow from the 0 to the 1
        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[0]->CurrentID) = 1;
        vec(vec.size() - 1) = -Current;

        mat.resize(mat.rows() + 1, mat.cols());
        vec.resize(vec.size() + 1);
        mat(mat.rows() - 1, edges[1]->CurrentID) = 1;
        vec(vec.size() - 1) = Current;
    }
};


// edge[0] is the Collector 
// edge[1] is the Base 
// edge[2] is the Emitter
/*class transistornode :public componentnode{
    public:
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

// Function for searching a tyingnode by name or creating it if it doenst exist yet

tyingnode* Gettyingnode(string name) {
    tyingnode* res = PairT[name];
    if (res == nullptr) {
        PairT[name] = res = new tyingnode(); //fancy syntax for a new tyingnode named pair which is equal to the name of Pair
    }
    return res;
}

double ReadValue(string s) {
    double x,y;
    x = stod(s);
    if (s.find("k"))
    {
        y = x * 1000;
        return y;
    }
    else if (s.find("Meg"))
    {
        y = x * 1000000;
        return y;
    }
    else if (s.find("μ"))
    {
        y = x / 1000000;
        return y;
    }
    else if (s.find("m"))
    {
        y = x / 1000;
        return y;
    }
    else if (s.find("p"))
    {
        y = x / 1000000000000;
        return y;
    }
    else if(s.find("n"))
    {
        y = x / 1000000000;
        return y;
    }
    else
    {
        return x;
    }
}

int main()
{

    //Reading in 

    ifstream something("netlist.txt");
    string netlist_line;
    while (getline(something, netlist_line))
    {
        string component, nodefrom, nodeto, value;
        stringstream iss(netlist_line);
        iss >> component >> nodefrom >> nodeto >> value;
        tyingnode* tyingfrom = Gettyingnode(nodefrom);
        tyingnode* tyingto = Gettyingnode(nodeto);
        componentnode* componentnodepointer;
        switch (component[0]) {
        case 'V':
            componentnodepointer = new voltagesourcenode(iss);
            break;
        case 'I':
            componentnodepointer = new currentsourcenode(iss);
            break;
        case 'R':
            componentnodepointer = new resistornode(iss);
            break;
        case 'L':
            componentnodepointer = new Inductornode(iss); // need constructor
            break;
        case 'C':
            componentnodepointer = new Capacitornode(iss); // need constructor
            break;
        }
        new edge(componentnodepointer, tyingfrom);
        new edge(componentnodepointer, tyingto);

        PairC[component] = componentnodepointer;

    }

    //Solving the matrix of doom

    Matrix mat;
    Vector vec;
    Vector res;
    res = mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vec);
    cout << res << endl;
    return 0;
}


/*
A test circuit to demonstrate SPICE syntax
V1 N003 0 SINE(2 1 1000)
R1 N001 N003 1k
C1 N001 0 1μ
I1 0 N004 0.1
D1 N004 N002 D
L1 N002 N001 1m
R2 N002 N001 1Meg
Q1 N003 N001 0 NPN
.tran 0 10ms 0 1us
.end
*/
