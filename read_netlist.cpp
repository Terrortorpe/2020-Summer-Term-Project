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

    for (const auto &a: data)
    {
        for (const auto &b: a)
        {
            cout<<b<<" ";
        }
        cout<<'\n';
    }
        
}
    