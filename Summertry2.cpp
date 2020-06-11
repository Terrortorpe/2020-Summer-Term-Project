#include <vector>
#include "Eigen/Dense"
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <regex>

using namespace std;


class componentnode;
class tyingnode;

typedef Eigen::MatrixXcd Matrix;
typedef Eigen::VectorXcd Vector;

complex<double> I = {0,1};

map<string,tyingnode*> PairT;
map<string,componentnode*> PairC;

vector<edge*> Edges;
//Indices for voltages and currents
int IDcount;
map<int, string> id_name;
//Indices for derivatives
int dIDcount;
//double ac_omega;

struct edge {
    componentnode* A;
    tyingnode* B;
    //Cache from previous result
    double Current = 0; //positive if going INTO the tyingnode
    int CurrentID;
    int dCurrentID;

    edge(componentnode* a_, tyingnode* b_){
        A=a_;
        B=b_;
        CurrentID =  IDcount;
        id_name[CurrentID] = A->name + "->" + B->name;
        IDcount++;
        dCurrentID =  dIDcount;
        dIDcount++;
        A->edges.push_back(this);
        B->edges.push_back(this);

        Edges.push_back(this);
    }

    void cache_prev_result(double dt, Vector& res) {
        Current = res(CurrentID);
    }
};

class node {
public:
    string name;
    vector<edge*> edges;
    node(string _name) : name(_name) {

    }
    virtual void constraint_cons(Matrix& mat, Vector& vec, double t){
        throw 1;
    }

    virtual void constraint_td(Matrix& mat, Vector& vec, double t){
        throw 1;
    }

    virtual void cache_prev_result(double dt, Vector& res) {
        throw 1;
    }
};

class componentnode : public node {
public:
    componentnode(string _name) : node(_name) {

    }

    virtual void constraint_cons(Matrix& mat, Vector& vec, double t){
        throw 1;
    }

    virtual void constraint_td(Matrix& mat, Vector& vec, double t){
        throw 1;
    }

    /*virtual void prepare_ac_analysis(double omega) {
        throw 1;
    }

    virtual void prepare_dc_analysis() {
        throw 1;
    }*/

    virtual void cache_prev_result(double dt, Vector& res) {
        throw 1;
    }
};

class tyingnode : public node {
public:
    double Voltage = 0;
    int VoltageID;
    int dVoltageID;

    tyingnode(string _name) : node(_name) {
        VoltageID = IDcount;
        id_name[VoltageID] = name;
        IDcount++;
        dVoltageID = IDcount;
        IDcount++;
    }
    
    void constraint_cons(Matrix& mat, Vector& vec, double t){
        //Sum current
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){
        //Sum current derivative
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->dCurrentID) = 1;
        }
        vec(vec.size()-1) = 0;
    }

    void cache_prev_result(double dt, Vector& res) {
        Voltage = res(VoltageID);
    }
};

class Capacitornode :public componentnode{
    public:
    double Capacitance;
    //complex<double> Impedance;

    Capacitornode(string _name, stringstream &ss) : componentnode(_name) {
        string s;
        s = ss.str();
        Capacitance = ReadValue(s);
    }

    void constraint_cons(Matrix& mat, Vector& vec, double t){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //Act as voltage source
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = 1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
        vec(vec.size()-1) = edges[1]->B->Voltage - edges[0]->B->VoltageID;
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){
        //Current derivative sums to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //V' = CI
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->dVoltageID) = 1.0;      // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistancce
        mat(mat.rows()-1, edges[1]->B->dVoltageID) = -1.0;
        mat(mat.rows()-1, edges[0]->CurrentID) = Capacitance;
        vec(vec.size()-1) = 0;
    }

    /*void prepare_ac_analysis(double omega) {
        Impedance = 1.0/(I * omega * Capacitance);
    }*/

    void cache_prev_result(double dt, Vector& res) {

    }
};

class Inductornode :public componentnode{
    public:
    double Inductance;
    //complex<double> Impedance;

    Inductornode(string _name, stringstream &ss) : componentnode(_name) {
        string s;
        s = ss.str();
        Inductance = ReadValue(s);
    }

    void constraint_cons(Matrix& mat, Vector& vec, double t){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //Current stays constant
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->CurrentID) = 1;
        vec(vec.size()-1) = edges[0]->Current;
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){
        //Current derivative flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->dCurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //V = L*I'
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[0]->dCurrentID) = Inductance;
        vec(vec.size()-1) = 0;


    }

    /*void prepare_ac_analysis(double omega) {
        Impedance = I * omega * Inductance;
    }*/

    void cache_prev_result(double dt, Vector& res) {

    }
};

class resistornode :public componentnode{
    public:
    double Resistance;
    resistornode(string _name, stringstream &ss) : componentnode(_name) {
        string s;
        s = ss.str();
        Resistance = ReadValue(s);
    };

    void constraint_cons(Matrix& mat, Vector& vec, double t){

        //Current flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->CurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //V = I R
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->VoltageID) = 1;
        mat(mat.rows()-1, edges[1]->B->VoltageID) = -1;
        mat(mat.rows()-1, edges[0]->CurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
        vec(vec.size()-1) = 0;
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){
        //Current derivative flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->dCurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //V' = I' R
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->dVoltageID) = 1;
        mat(mat.rows()-1, edges[1]->B->dVoltageID) = -1;
        mat(mat.rows()-1, edges[0]->dCurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
        vec(vec.size()-1) = 0;
    }

    /*void prepare_ac_analysis(double omega) {
        
    }*/

    void cache_prev_result(double dt, Vector& res) {

    }
};


class voltagesourcenode :public componentnode{
    public:
    bool time_dep = false;
    double Offset;
    double Amplitude;
    double Frequency;

    voltagesourcenode(string _name, stringstream &ss) : componentnode(_name) {
        string s;
        s = ss.str();
        if(s.rfind("SINE") == 0){
            regex rgx("[(](.*)[)]");
            smatch match;
            if (regex_search(s.begin(), s.end(), match, rgx)) {
                ss.clear();
                ss << match[1];
                ss >> Offset >> Amplitude >> Frequency;
                time_dep = true;
            } else {
                throw 1;
            }
        } else {
            Offset = 0;
            Amplitude = ReadValue(s);
            Frequency = 0;
            time_dep = false;
        }
    }


    void constraint_cons(Matrix& mat, Vector& vec, double t){
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
        vec(vec.size()-1) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude; //V(1)-V(0)=U
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){
        //Current derivative flowing in and out has to sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        for(auto && edge:edges){
            mat(mat.rows()-1, edge->dCurrentID) = 1;
        }
        vec(vec.size()-1) = 0;

        //Voltage derivative must match
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->B->dVoltageID) = -1;
        mat(mat.rows()-1, edges[1]->B->dVoltageID) = 1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
        vec(vec.size()-1) = 0; //V'(1)-V'(0)=0
    }

    void cache_prev_result(double dt, Vector& res) {

    }
};

class currentsourcenode :public componentnode{
public:
    bool time_dep = false;
    double Offset;
    double Amplitude;
    double Frequency;

    currentsourcenode(string _name, stringstream &ss) : componentnode(_name) {
        string s;
        s = ss.str();
        if(s.rfind("SINE") == 0){
            regex rgx("[(](.*)[)]");
            smatch match;
            if (regex_search(s.begin(), s.end(), match, rgx)) {
                ss.clear();
                ss << match[1];
                ss >> Offset >> Amplitude >> Frequency;
                time_dep = true;
            } else {
                throw 1;
            }
        } else {
            Offset = 0;
            Amplitude = ReadValue(s);
            Frequency = 0;
            time_dep = false;
        }
    }
    void constraint_cons(Matrix& mat, Vector& vec, double t){

        //Current source makes current flow from the 0 to the 1
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->CurrentID) = 1;
        vec(vec.size()-1) = -time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude;

        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[1]->CurrentID) = 1;
        vec(vec.size()-1) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude;
    }

    void constraint_td(Matrix& mat, Vector& vec, double t){

        //Current changes sum to 0
        mat.resize(mat.rows()+1,mat.cols());
        vec.resize(vec.size()+1);
        mat(mat.rows()-1, edges[0]->dCurrentID) = 1;
        mat(mat.rows()-1, edges[1]->dCurrentID) = -1;
        vec(vec.size()-1) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude;
    }

    void cache_prev_result(double dt, Vector& res) {

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

tyingnode* Gettyingnode(string name){
    tyingnode* res = PairT[name];
    if (res == nullptr){
        PairT[name] = res = new tyingnode(name); //fancy syntax for a new tyingnode named pair which is equal to the name of Pair
    }
    return res;
}

double ReadValue(string s){
    double x,y;
    x = stod(s);
    if (s.find("k"))
    {
        y = x * 1e3;
        return y;
    }
    else if (s.find("Meg"))
    {
        y = x * 1e6;
        return y;
    }
    else if (s.find("μ"))
    {
        y = x * 1e-6;
        return y;
    }
    else if (s.find("m"))
    {
        y = x * 1e-3;
        return y;
    }
    else if (s.find("p"))
    {
        y = x * 1e-12;
        return y;
    }
    else if(s.find("n"))
    {
        y = x * 1e-9;
        return y;
    }
    else
    {
        return x;
    }
}

int main()
{
    bool transient;
    
    //Reading in 

    ifstream something("netlist.txt");
    string netlist_line;
    while (getline(something,netlist_line))
    {
        string component, nodefrom, nodeto;
        stringstream iss(netlist_line);
        iss >> component >> nodefrom >> nodeto;
        tyingnode* tyingfrom = Gettyingnode(nodefrom);
        tyingnode* tyingto = Gettyingnode(nodeto);
        componentnode* componentnodepointer;
        switch(component[0]){
            case 'V':
            componentnodepointer = new voltagesourcenode(component, iss);
            break;
            case 'I':
            componentnodepointer = new currentsourcenode(component, iss);
            break;
            case 'R':
            componentnodepointer = new resistornode(component, iss);
            break;
            case 'L':
            componentnodepointer = new Inductornode(component, iss);
            break;
            case 'C':
            componentnodepointer = new Capacitornode(component, iss);
            break;
        }
        new edge(componentnodepointer, tyingfrom);
        new edge(componentnodepointer, tyingto);

        PairC[component] = componentnodepointer;

    }

    fstream textfile("out.txt");
 
    double stoptime, timestep;
    textfile.open("simulation_results.txt");
 
    for(double time = 0; time = stoptime; time = time+timestep)
    {
        //Step 1: Time-Indep solution
        {
            Matrix mat = Matrix(0, IDcount);
            Vector vec = Vector(0);
            Vector res = Vector(IDcount);

            for(auto&& tnode : PairT) {
                tnode.second->constraint_cons(mat, vec, time);
            }

            for(auto&& cnode : PairC) {
                cnode.second->constraint_cons(mat, vec, time);
            }

            //Solving the matrix of doom
            res = mat.bdcSvd(Eigen::ComputeThinU|Eigen::ComputeThinV).solve(vec);

            //Cache results
            for(auto&& tnode : PairT) {
                tnode.second->cache_prev_result(timestep, res);
            }

            for(auto&& cnode : PairC) {
                cnode.second->cache_prev_result(timestep, res);
            }

            for(auto&& edge : Edges) {
                edge->cache_prev_result(timestep, res);
            }

            textfile << time << ":" << res <<endl;
        }
        //Time se
        {

        }

    }

    return 0;

}


/* 
A test circuit to demonstrate SPICE syntax 
V1 N003 0 SINE(2 1 1000)
V2 N005 0 5
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




