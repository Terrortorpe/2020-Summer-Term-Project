#include <sstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main()
{

    vector<vector<string>> data;
    string netlist_line;

    while (getline(cin,netlist_line))
    {
        data.push_back({});
        istringstream iss(netlist_line);
        string netlist_word;
        while(iss>>netlist_word)
        {
            data.back().push_back(netlist_word);
        }
    }

    for (const vector<string> &a: data)
    {
     //   cout<< a[0] << endl; 
        if(a[0].find("V") != string::npos){
            cout<< "Creating Voltage Source" << endl;
        }else if(a[0].find("I") != string::npos){
            cout<< "Creating Current Source" << endl;
        }else if(a[0].find("*") != string::npos){
            cout<< "Comment" << endl;
        }else if(a[0].find("R") != string::npos){
            cout<< "Creating Resistor" << endl;
        }else if(a[0].find("C") != string::npos){
            cout<< "Creating Capacitor" << endl;
        }else if(a[0].find("L") != string::npos){
            cout<< "Creating Inductor" << endl;
        }else if(a[0].find("D") != string::npos){
            cout<< "Creating Diode" << endl;
        }else if(a[0].find("Q") != string::npos){
            cout<< "Creating Transistor" << endl;

        }else if(a[0].find(".tran") != string::npos){
            cout<< "transient simulation" << endl;
        }else if(a[0].find(".end") != string::npos){
            cout<< "end" << endl;

        }else{
            cout<< "Error - unidentified component or simulation: " << a[0] << endl;
        }
    
    }
    /* this code just prints the netlist out, like a test
    for (const auto &a: data)
    {
        for (const auto &b: a)
        {
            cout<<b<<" ";
        }
        cout<<'\n';
    }
     */   
  
}
    