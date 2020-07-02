#include <iostream>
#include <fstream>
#include <unordered_set>

using namespace std;


void read_data(string path, vector< vector<int> >& node2hyperedge, vector< vector<int> >& hyperedge2node, vector< unordered_set<int> >& hyperedge2node_set)
{
	ifstream graphFile(path.c_str());
	string line;
	int num_hyperedges = 0;
	while (getline(graphFile, line))
	{
		size_t pos = 0;
		unordered_set<int> tokens_set;
		vector<int> tokens;
		bool EOL = false;
		int idx;
		while (!EOL)
		{
			int pos = line.find(",");
			if(pos == string::npos){
				pos = line.size();
				EOL = true;
				idx = stoi(line);
			}else{
				idx = stoi(line.substr(0, pos));
				line.erase(0, pos + 1);
			}
			while(idx >= node2hyperedge.size()){
				node2hyperedge.push_back(vector<int>());
			}
			if(node2hyperedge[idx].empty() || node2hyperedge[idx].back() != num_hyperedges){
				node2hyperedge[idx].push_back(num_hyperedges);
				tokens.push_back(idx);
				tokens_set.insert(idx);
			}
		}
		hyperedge2node.push_back(tokens);
		hyperedge2node_set.push_back(tokens_set);
		num_hyperedges++;
	}
	graphFile.close();
}
