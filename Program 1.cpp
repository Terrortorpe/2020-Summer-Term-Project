// Program 1.cpp : Bu dosya 'main' işlevi içeriyor. Program yürütme orada başlayıp biter.
//

#include <iostream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;
    class Component
    {
    private:
        string typeName;
        string value;
        string firstNode;
        string secondNode;
        string thirdNode;
    public:
        Component();
        ~Component();
        void setValue(string pValue);
        string getValue();
        void setType(string pType);
        string getType();
        void setFirstNode(string pNode);
        void setSecondNode(string pNode);
        void setThirdNode(string pNode);
        string getNode();
    };
    Component::Component()
    {

    }
    Component::~Component()
    {

    }

    void Component::setValue(string pValue)
    {
        value = pValue;
    }

    string Component::getValue()
    {
        return value;
    }
    
    void Component::setType(string pType)
    {
        pType[0] = toupper(pType[0]);
        if (pType[0] == ('V'|'I'|'R'|'C'|'L'|'D'|'Q'))
        {
            typeName = pType;
        }
        else
        {
            cout << "Undefined Type";
        }
    }
    
    string Component::getType()
    {
        return typeName;
    }
    void Component::setFirstNode(string pNode)
    {
        firstNode = pNode;
    }
    void Component::setSecondNode(string pNode)
    {
        secondNode = pNode;
    }
    void Component::setThirdNode(string pNode)
    {
        thirdNode = pNode;
    }

    string Component::getNode()
    {
        int x;
        cout << "Which Node value do you want?";
        cin >> x;
        switch (x)
        {
        case 1:
            return firstNode;
            break;
        case 2:
            return secondNode;
            break;
        case 3:
            return thirdNode;
            break;
        }
    }

    int main()
    {
        vector<Component> v;
        vector<Component> r;
        vector<Component> c;
        vector<Component> l;
        vector<Component> d;
        vector<Component> q;
        vector<Component> i;
        vector<vector<string>> data;
        string netlist_line;

        while (getline(cin, netlist_line))
        {
            data.push_back({});
            istringstream iss(netlist_line);
            string netlist_word;
            while (iss >> netlist_word)
            {
                data.back().push_back(netlist_word);
            }
        }

        for (const vector<string>& a : data)
        {;
               cout<< a[0] << endl; 
            if (a[0].find("V") != string::npos) {
                cout << "Creating Voltage Source" << endl;

            }
            else if (a[0].find("I") != string::npos) {
                cout << "Creating Current Source" << endl;
                Component i;
                i.setType("I");
                i.setValue(a[3]);
                i.setFirstNode(a[1]);
                i.setSecondNode(a[2]);
            }
            else if (a[0].find("*") != string::npos) {
                cout << "Comment" << endl;
            }
            else if (a[0].find("R") != string::npos) {
                cout << "Creating Resistor" << endl;
                Component r;
                r.setType("R");
                r.setValue(a[3]);
                r.setFirstNode(a[1]);
                r.setSecondNode(a[2]);
            }
            else if (a[0].find("C") != string::npos) {
                cout << "Creating Capacitor" << endl;
            }
            else if (a[0].find("L") != string::npos) {
                cout << "Creating Inductor" << endl;
            }
            else if (a[0].find("D") != string::npos) {
                cout << "Creating Diode" << endl;
            }
            else if (a[0].find("Q") != string::npos) {
                cout << "Creating Transistor" << endl;

            }
            else if (a[0].find(".tran") != string::npos) {
                cout << "transient simulation" << endl;
            }
            else if (a[0].find(".end") != string::npos) {
                cout << "end" << endl;

            }
            else {
                cout << "Error - unidentified component or simulation: " << a[0] << endl;
            }

        }
         //this code just prints the netlist out, like a test
        for (const auto &a: data)
        {
            for (const auto &b: a)
            {
                cout<<b<<" ";
            }
            cout<<'\n';
        }
         
        system("PAUSE");

    }
    

