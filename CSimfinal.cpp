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
struct edge;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

//complex<double> I = {0,1};

map<string, tyingnode*> PairT;
map<string, componentnode*> PairC;

vector<edge*> Edges;
//Indices for voltages and currents
int IDcount;
map<int, string> id_name;
//Indices for derivatives
int dIDcount;
//double ac_omega;

double ReadValue(string s) {
	double x, y;
	x = stod(s);
	if (s.find("k") != string::npos) {
		y = x * 1e3;
		return y;
	} else if (s.find("Meg") != string::npos) {
		y = x * 1e6;
		return y;
	} else if (s.find("u") != string::npos) {
		y = x * 1e-6;
		return y;
	} else if (s.find("m") != string::npos) {
		y = x * 1e-3;
		return y;
	} else if (s.find("p") != string::npos) {
		y = x * 1e-12;
		return y;
	} else if (s.find("n") != string::npos) {
		y = x * 1e-9;
		return y;
	} else {
		return x;
	}
}

struct edge {
	componentnode* A;
	tyingnode* B;
	//Cache from previous result
	double Current = 0; //positive if going INTO the tyingnode
	int CurrentID;
	int dCurrentID;

	edge(componentnode* a_, tyingnode* b_);

	void cache_prev_result(double dt, Vector& res, Vector& rest);
};

class node {
public:
	string name;
	vector<edge*> edges;
	node(string _name);
	virtual int constraint_cons_N();
	virtual void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	virtual void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual void cache_prev_result(double dt, Vector& res, Vector& rest);
};

class componentnode : public node {
public:
	componentnode(string _name);

	virtual int constraint_cons_N();
	virtual void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	virtual void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	/*virtual void prepare_ac_analysis(double omega) {
			throw 1;
	}

	virtual void prepare_dc_analysis() {
			throw 1;
	}*/

	virtual void cache_prev_result(double dt, Vector& res, Vector& rest);
};

class tyingnode : public node {
public:
	double Voltage = 0;
	int VoltageID;
	int dVoltageID;

	tyingnode(string _name);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	void cache_prev_result(double dt, Vector& res, Vector& rest);
};

class Capacitornode :public componentnode {
public:
	double Capacitance;
	//complex<double> Impedance;

	Capacitornode(string _name, stringstream& ss);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	/*void prepare_ac_analysis(double omega) {
			Impedance = 1.0/(I * omega * Capacitance);
	}*/
};

class Inductornode :public componentnode {
public:
	double Inductance;
	//complex<double> Impedance;

	Inductornode(string _name, stringstream& ss);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	/*void prepare_ac_analysis(double omega) {
			Impedance = I * omega * Inductance;
	}*/
};

class resistornode :public componentnode {
public:
	double Resistance;
	resistornode(string _name, stringstream& ss);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);

	/*void prepare_ac_analysis(double omega) {

	}*/
};

class voltagesourcenode :public componentnode {
public:
	bool time_dep = false;
	double Offset;
	double Amplitude;
	double Frequency;

	voltagesourcenode(string _name, stringstream& ss);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);
};

class currentsourcenode :public componentnode {
public:
	bool time_dep = false;
	double Offset;
	double Amplitude;
	double Frequency;

	currentsourcenode(string _name, stringstream& ss);

	virtual int constraint_cons_N();
	void constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt);

	virtual int constraint_td_N();
	void constraint_td(Matrix& mat, Vector& vec, double t, int& cnt);
};

//////////////////////////////////////////////////////////////////////////////////////////

edge::edge(componentnode* a_, tyingnode* b_) {
	A = a_;
	B = b_;
	CurrentID = IDcount;
	id_name[CurrentID] = A->name + "->" + B->name;
	IDcount++;
	dCurrentID = dIDcount;
	dIDcount++;
	A->edges.push_back(this);
	B->edges.push_back(this);

	Edges.push_back(this);
}

void edge::cache_prev_result(double dt, Vector& res, Vector& rest) {
	Current = res(CurrentID) + dt * rest(dCurrentID);
}

node::node(string _name) : name(_name) {

}

int node::constraint_cons_N() {
	throw 1;
}
void node::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {
	throw 1;
}

int node::constraint_td_N() {
	throw 1;
}
void node::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	throw 1;
}

void node::cache_prev_result(double dt, Vector& res, Vector& rest) {
	throw 1;
}

componentnode::componentnode(string _name) : node(_name) {

}

int componentnode::constraint_cons_N() {
	throw 1;
}
void componentnode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {
	throw 1;
}

int componentnode::constraint_td_N() {
	throw 1;
}
void componentnode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	throw 1;
}

/*virtual void prepare_ac_analysis(double omega) {
		throw 1;
}

virtual void prepare_dc_analysis() {
		throw 1;
}*/

void componentnode::cache_prev_result(double dt, Vector& res, Vector& rest) {

}



tyingnode::tyingnode(string _name) : node(_name) {
	VoltageID = IDcount;
	id_name[VoltageID] = name;
	IDcount++;
	dVoltageID = dIDcount;
	dIDcount++;
}

int tyingnode::constraint_cons_N() {
	return 1;
}
void tyingnode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Sum current
	for (auto&& edge : edges) {
		mat(cnt, edge->CurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;
}

int tyingnode::constraint_td_N() {
	return 1;
}
void tyingnode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Sum current derivative
	for (auto&& edge : edges) {
		mat(cnt, edge->dCurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;
}

void tyingnode::cache_prev_result(double dt, Vector& res, Vector& rest) {
	Voltage = res(VoltageID) + dt * rest(dVoltageID);
}



Capacitornode::Capacitornode(string _name, stringstream& ss) : componentnode(_name) {
	string s;
	ss >> s;
	Capacitance = ReadValue(s);
	cout << "Loaded " << _name << " " << Capacitance << endl;
}

int Capacitornode::constraint_cons_N() {
	return 2;
}
void Capacitornode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {

	//Current flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->CurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//Act as voltage source
	mat(cnt, edges[0]->B->VoltageID) = -1;
	mat(cnt, edges[1]->B->VoltageID) = 1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
	vec(cnt) = edges[1]->B->Voltage - edges[0]->B->Voltage;
	++cnt;
}

int Capacitornode::constraint_td_N() {
	return 2;
}
void Capacitornode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current derivative sums to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->dCurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//V' = I/C
	mat(cnt, edges[0]->B->dVoltageID) = -1.0;      // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistancce
	mat(cnt, edges[1]->B->dVoltageID) = 1.0;
	vec(cnt) = edges[0]->Current / Capacitance;
	++cnt;
}



Inductornode::Inductornode(string _name, stringstream& ss) : componentnode(_name) {
	string s;
	ss >> s;
	Inductance = ReadValue(s);
	cout << "Loaded " << _name << " " << Inductance << endl;
}

int Inductornode::constraint_cons_N() {
	return 2;
}
void Inductornode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {

	//Current flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->CurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//Current stays constant
	mat(cnt, edges[0]->CurrentID) = 1;
	vec(cnt) = edges[0]->Current;
	++cnt;
}

int Inductornode::constraint_td_N() {
	return 2;
}
void Inductornode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current derivative flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->dCurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//V1-V0 = L*I0'
	mat(cnt, edges[0]->dCurrentID) = Inductance;
	vec(cnt) = edges[1]->B->Voltage - edges[0]->B->Voltage;
	++cnt;
}



resistornode::resistornode(string _name, stringstream& ss) : componentnode(_name) {
	string s;
	ss >> s;
	Resistance = ReadValue(s);
	cout << "Loaded " << _name << " " << Resistance << endl;
};

int resistornode::constraint_cons_N() {
	return 2;
}
void resistornode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {

	//Current flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->CurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//V = I R
	mat(cnt, edges[0]->B->VoltageID) = 1;
	mat(cnt, edges[1]->B->VoltageID) = -1;
	mat(cnt, edges[0]->CurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
	vec(cnt) = 0;
	++cnt;
}

int resistornode::constraint_td_N() {
	return 2;
}
void resistornode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current derivative flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->dCurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//V' = I' R
	mat(cnt, edges[0]->B->dVoltageID) = -1;
	mat(cnt, edges[1]->B->dVoltageID) = 1;
	mat(cnt, edges[0]->dCurrentID) = Resistance;    // V(0)-----E(0)-----R------E(1)----V(1)    V(0)-V(1)+current*resistance=0
	vec(cnt) = 0;
	++cnt;
}



voltagesourcenode::voltagesourcenode(string _name, stringstream& ss) : componentnode(_name) {
	string s;
	getline(ss, s);
	if (s.rfind(" SINE") == 0) {
		regex rgx("[(](.*)[)]");
		smatch match;
		if (regex_search(s, match, rgx)) {
			ss = stringstream(match[1]);
			string ns;
			ss >> ns;
			Offset = ReadValue(ns);
			ss >> ns;
			Amplitude = ReadValue(ns);
			ss >> ns; 
			Frequency = ReadValue(ns);
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

int voltagesourcenode::constraint_cons_N() {
	return 2;
}
void voltagesourcenode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->CurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//Voltage source has larger voltage on the edge[0] side
	mat(cnt, edges[0]->B->VoltageID) = 1;
	mat(cnt, edges[1]->B->VoltageID) = -1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
	vec(cnt) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude; //V(1)-V(0)=U
	++cnt;
}

int voltagesourcenode::constraint_td_N() {
	return 2;
}
void voltagesourcenode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current derivative flowing in and out has to sum to 0
	for (auto&& edge : edges) {
		mat(cnt, edge->dCurrentID) = 1;
	}
	vec(cnt) = 0;
	++cnt;

	//Voltage has derivative
	mat(cnt, edges[0]->B->dVoltageID) = 1;
	mat(cnt, edges[1]->B->dVoltageID) = -1;          //V(0)----E(0)-----Source(U)------E(1)-------V(1)
	vec(cnt) = time_dep ? (Frequency * 2 * 3.1415926535897932 * Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : 0; //V(1)-V(0)=U
	++cnt;
}



currentsourcenode::currentsourcenode(string _name, stringstream& ss) : componentnode(_name) {
	string s;
	getline(ss, s);
	if (s.rfind(" SINE") == 0) {
		regex rgx("[(](.*)[)]");
		smatch match;
		if (regex_search(s, match, rgx)) {
			ss = stringstream(match[1]);
			string ns;
			ss >> ns;
			Offset = ReadValue(ns);
			ss >> ns;
			Amplitude = ReadValue(ns);
			ss >> ns;
			Frequency = ReadValue(ns);
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

int currentsourcenode::constraint_cons_N() {
	return 2;
}
void currentsourcenode::constraint_cons(Matrix& mat, Vector& vec, double t, int& cnt) {

	//Current source makes current flow from the 0 to the 1
	mat(cnt, edges[0]->CurrentID) = -1;
	vec(cnt) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude;
	++cnt;


	mat(cnt, edges[1]->CurrentID) = 1;
	vec(cnt) = time_dep ? (Offset + Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : Amplitude;
	++cnt;
}

int currentsourcenode::constraint_td_N() {
	return 2;
}
void currentsourcenode::constraint_td(Matrix& mat, Vector& vec, double t, int& cnt) {
	//Current changes sum to 0
	mat(cnt, edges[0]->dCurrentID) = 1;
	mat(cnt, edges[1]->dCurrentID) = -1;
	vec(cnt) = 0;
	++cnt;

	//Current changes
	mat(cnt, edges[1]->dCurrentID) = 1;
	vec(cnt) = time_dep ? (Frequency * 2 * 3.1415926535897932 * Amplitude * sin(Frequency * 2 * 3.1415926535897932 * t)) : 0;
	++cnt;
}

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
		PairT[name] = res = new tyingnode(name); //fancy syntax for a new tyingnode named pair which is equal to the name of Pair
	}
	return res;
}

int main() {
	double timestep = 0.001, stoptime = 1, starttime = 0;

	//Reading in 

	ifstream something("netlist.txt");
	string netlist_line;
	while (getline(something, netlist_line)) {
		if (netlist_line[0] != '*') {
			if (netlist_line[0] == '.') {
				string command;
				stringstream iss(netlist_line);
				iss >> command;
				if (command == ".tran") {
					string tmp;
					iss >> tmp;
					stoptime = ReadValue(tmp);
					iss >> tmp;
					timestep = ReadValue(tmp);
					cout << "Stop: " << stoptime << " Step:" << timestep << endl;
				} else {
					cerr << "Unknown command " << command << endl;
				}
			} else {
				string component, nodefrom, nodeto;
				stringstream iss(netlist_line);
				iss >> component >> nodefrom >> nodeto;
				tyingnode* tyingfrom = Gettyingnode(nodefrom);
				tyingnode* tyingto = Gettyingnode(nodeto);
				componentnode* componentnodepointer;
				switch (component[0]) {
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
		}
	}

	ofstream textfile("out.txt");

	textfile << "time";

	for (auto&& idn : id_name) {
		textfile << "\t" << idn.second;
	}
	textfile << endl;

	for (double time = starttime; time < stoptime; time = time + timestep) {
		int cnt_cons = 0;

		for (auto&& tnode : PairT) {
			cnt_cons += tnode.second->constraint_cons_N();
		}

		for (auto&& cnode : PairC) {
			cnt_cons += cnode.second->constraint_cons_N();
		}

		//Gnd
		cnt_cons++;

		//Step 1: Time-Indep solution
		Matrix mat = Matrix::Zero(cnt_cons, IDcount);
		Vector vec = Vector::Zero(cnt_cons);
		Vector res = Vector::Zero(IDcount);

		//Ground
		mat(0, PairT["0"]->VoltageID) = 1;
		vec(0) = 0;

		int i = 1;

		for (auto&& tnode : PairT) {
			tnode.second->constraint_cons(mat, vec, time, i);
		}

		for (auto&& cnode : PairC) {
			cnode.second->constraint_cons(mat, vec, time, i);
		}

		assert(i == cnt_cons);

		//Solving the matrix of doom
		res = mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vec);

		Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, "", "\t");
		textfile << time << "\t" << res.format(fmt) << endl;

		//Step 2: Time step
		int cnt_td = 0;

		for (auto&& tnode : PairT) {
			cnt_td += tnode.second->constraint_td_N();
		}

		for (auto&& cnode : PairC) {
			cnt_td += cnode.second->constraint_td_N();
		}
		++cnt_td;

		mat = Matrix::Zero(cnt_td, dIDcount);
		vec = Vector::Zero(cnt_td);
		Vector rest = Vector::Zero(dIDcount);

		//Ground
		mat(0, PairT["0"]->dVoltageID) = 1;
		vec(0) = 0;

		i = 1;

		for (auto&& tnode : PairT) {
			tnode.second->constraint_td(mat, vec, time, i);
		}

		for (auto&& cnode : PairC) {
			cnode.second->constraint_td(mat, vec, time, i);
		}

		assert(i == cnt_td);

		//Solving the matrix of doom
		rest = mat.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(vec);

		//Cache results
		for (auto&& tnode : PairT) {
			tnode.second->cache_prev_result(timestep, res, rest);
		}

		for (auto&& cnode : PairC) {
			cnode.second->cache_prev_result(timestep, res, rest);
		}

		for (auto&& edge : Edges) {
			edge->cache_prev_result(timestep, res, rest);
		}
	}

	textfile.close();

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




