#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;

struct hyperwedge{
	int a, b, C_ab;	
};

inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

long long edge_size_limit, edge_capacity_left;
vector<long long> upd_time_proj;
priority_queue< pair<int, int>, vector< pair<int, int> >, greater< pair<int, int> > > pq;


void get_adj(int hyperedge_a, int deg_a, vector< vector< pair<int, int> > >& hyperedge_adj, vector< unordered_map<int, int> >& hyperedge_inter, vector< vector<int> >& node2hyperedge, vector< vector<int> >& hyperedge2node){
	int deg_a_curr = 0;
	hyperedge_adj[hyperedge_a].resize(deg_a);
	for (const int &node: hyperedge2node[hyperedge_a]){
		for (const int &hyperedge_b: node2hyperedge[node]){
			if (hyperedge_b == hyperedge_a) continue;
			if ((upd_time_proj[hyperedge_b] >> 31) ^ hyperedge_a){
				upd_time_proj[hyperedge_b] = ((long long)hyperedge_a << 31) + deg_a_curr;
				hyperedge_adj[hyperedge_a][deg_a_curr++] = {hyperedge_b, 0};
			}else if((int)(upd_time_proj[hyperedge_b] & 0x7FFFFFFFLL) == deg_a_curr){
				upd_time_proj[hyperedge_b] = ((long long)hyperedge_a << 31) + deg_a_curr;
				hyperedge_adj[hyperedge_a][deg_a_curr++] = {hyperedge_b, 0};
			}
			hyperedge_adj[hyperedge_a][(int)(upd_time_proj[hyperedge_b] & 0x7FFFFFFFLL)].second++;
		}
	}
	hyperedge_inter[hyperedge_a].rehash(deg_a);
	for (int i = 0; i < deg_a; i++){
		int hyperedge_b = hyperedge_adj[hyperedge_a][i].first, C_ab = hyperedge_adj[hyperedge_a][i].second;
		hyperedge_inter[hyperedge_a].insert({hyperedge_b, C_ab});
	}
	pq.push({deg_a, hyperedge_a});
}

int main(int argc, char *argv[])
{
	clock_t start;
	clock_t run_start;
	int progress;
	
	long long sampling_size = stoi(argv[1]);
	long double mem_p = stod(argv[2]);

	string graphFile = "dblp_graph.txt";

	cout << "Sampling size: " << sampling_size << endl << endl;

	// Read data
	start = clock();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);

	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
	cout << "Reading data done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

	
	// h_motif counting via hyperwedge smapling
	start = clock(); run_start = clock();
	vector<bool> searched(E, false);
	vector<long long> h_motif(30, 0);
	vector<int> intersection(V, 0);
	vector< vector< pair<int, int> > > hyperedge_adj;
	vector< unordered_map<int, int> > hyperedge_inter;
	hyperedge_adj.resize(E);
	hyperedge_inter.resize(E);
	vector<int> upd_time(E, -1);
	upd_time_proj.resize(E);
	std::fill(upd_time_proj.begin(), upd_time_proj.end(), -1LL);
	
	vector<long long> degs_sum(E + 1, 0);
	for(int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		degs_sum[hyperedge_a + 1] = degs_sum[hyperedge_a];
		long long l_hyperedge_a = (long long)hyperedge_a;
		upd_time[hyperedge_a] = hyperedge_a;
		for (const int &node: hyperedge2node[hyperedge_a]){
			for (const int &hyperedge_b: node2hyperedge[node]){
				if (upd_time[hyperedge_b] ^ hyperedge_a){
					upd_time[hyperedge_b] = hyperedge_a;
					degs_sum[hyperedge_a + 1]++;
				}
			}
		}
	}
	edge_size_limit = (long long)((long double)degs_sum[E] * mem_p); edge_capacity_left = edge_size_limit;
	cout << edge_size_limit << endl;

	mt19937 gen(2020);
	uniform_real_distribution<> urd(0.0, 1.0);
	uniform_int_distribution<long long> dist(0, degs_sum[E] - 1);
	
	//sampling_size = degs_sum[E];
	long long max_batch_size = sampling_size;
	
	vector<long long> sampled_idx(sampling_size, 0);
	vector< vector<long long> > intermediate_buckets;
	vector< vector<long long> > final_buckets;
		
	intermediate_buckets.resize(E);
	final_buckets.resize((degs_sum[E] + E - 1) / E);
	
	for(long long batch_start_idx = 0; batch_start_idx < sampling_size; batch_start_idx += max_batch_size){
		std::fill(upd_time.begin(), upd_time.end(), -1LL);
		long long batch_size = min(max_batch_size, sampling_size - batch_start_idx);
		for (long long sample = 0; sample < batch_size; sample++){
			//sampled_idx[sample] = sample + batch_start_idx; //dist(gen);
			sampled_idx[sample] = dist(gen);
			intermediate_buckets[sampled_idx[sample] % E].push_back(sampled_idx[sample]);
		}

		for (int i = 0; i < E; i++){
			for(const long long &idx: intermediate_buckets[i]){
				final_buckets[idx / E].push_back(idx);
			}
			intermediate_buckets[i].clear();
		}

		int sampled_idx_pointer = 0;
		for(int i = 0; i < (int)final_buckets.size(); i++){
			for(const long long &idx: final_buckets[i]){
				sampled_idx[sampled_idx_pointer++] = idx;
			}
			final_buckets[i].clear();
		}

		unordered_map<int, int> sample_cnt;
		
		int hyperedge_a = -1;
		for (int sample = 0; sample < batch_size; sample++){
			if (sample % 10000 == 0)
				cout << "Sampling: " << sample << " / " << sampling_size << endl;
			
			while(sampled_idx[sample] >= degs_sum[hyperedge_a + 1]){ hyperedge_a++; }

			int deg_a = (int)(degs_sum[hyperedge_a + 1] - degs_sum[hyperedge_a]), size_a = (int)hyperedge2node[hyperedge_a].size();
			if(hyperedge_adj[hyperedge_a].size() ^ deg_a){
				while(edge_capacity_left < deg_a){
					pair<int, int> target = pq.top(); pq.pop();
					edge_capacity_left += target.first;
					hyperedge_adj[target.second].clear();
					hyperedge_inter[target.second].clear();
				} 
				get_adj(hyperedge_a, deg_a, hyperedge_adj, hyperedge_inter, node2hyperedge, hyperedge2node);
				edge_capacity_left -= deg_a;
			}

			int hyperedge_b_idx = sampled_idx[sample] - degs_sum[hyperedge_a];
			int hyperedge_b = hyperedge_adj[hyperedge_a][hyperedge_b_idx].first, C_ab = hyperedge_adj[hyperedge_a][hyperedge_b_idx].second;
			int deg_b = (int)(degs_sum[hyperedge_b + 1] - degs_sum[hyperedge_b]), size_b = (int)hyperedge2node[hyperedge_b].size();
			upd_time[hyperedge_a] = upd_time[hyperedge_b] = sample;

			if(hyperedge_adj[hyperedge_b].size() ^ deg_b){
				bool a_eliminated = false;
				while(edge_capacity_left < deg_b){
					pair<int, int> target = pq.top(); pq.pop();
					if(target.second ^ hyperedge_a){
						edge_capacity_left += target.first;
						hyperedge_adj[target.second].clear();
						hyperedge_inter[target.second].clear();
					}else{
						a_eliminated = true;
					}
				}
				if(a_eliminated){
					pq.push({deg_a, hyperedge_a});
				}
				get_adj(hyperedge_b, deg_b, hyperedge_adj, hyperedge_inter, node2hyperedge, hyperedge2node);
				edge_capacity_left -= deg_b;
			}

			int min_ab = hyperedge_a, max_ab = hyperedge_b;
			if (size_a > size_b) min_ab = hyperedge_b, max_ab = hyperedge_a;
			const auto &nodes = hyperedge2node_set[max_ab]; auto it_end = nodes.end(); int cnt = 0;
			for (const int &node: hyperedge2node[min_ab]){ if(nodes.find(node) != it_end) intersection[cnt++] = node;}

			for (int i = 0; i < deg_b; i++){
				int hyperedge_c = hyperedge_adj[hyperedge_b][i].first, C_bc = hyperedge_adj[hyperedge_b][i].second;
				if (upd_time[hyperedge_c] ^ sample){
					upd_time[hyperedge_c] = sample;
					int size_c = (int)hyperedge2node[hyperedge_c].size();
					int C_ca = 0, g_abc = 0;
					C_ca = hyperedge_inter[hyperedge_a][hyperedge_c];
					const auto &nodes = hyperedge2node_set[hyperedge_c]; auto it_end = nodes.end();	
					for (int k = 0; k < C_ab; k++){ if(nodes.find(intersection[k]) != it_end) g_abc++; }
					int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					h_motif[motif_index]++;
				}
			}

			for (int i = 0; i < deg_a; i++){
				int hyperedge_c = hyperedge_adj[hyperedge_a][i].first, C_ca = hyperedge_adj[hyperedge_a][i].second;
				if (upd_time[hyperedge_c] ^ sample){
					upd_time[hyperedge_c] = sample;
					int size_c = (int)hyperedge2node[hyperedge_c].size();
					int C_bc = 0, g_abc = 0;
					int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					h_motif[motif_index]++;
				}
			}
		}
	}

	int index = 0;
	vector<long double> h_motif_final(30, 0);
	for (int i = 0; i < 30; i++){
		h_motif_final[i] = (long double)h_motif[i];
		h_motif_final[i] *= (long double)degs_sum[E] / sampling_size;
		if (20 <= i && i <= 25)
			h_motif_final[i] /= 4.0;
		else
			h_motif_final[i] /= 6.0;
		if (i == 0 || i == 1 || i == 4 || i == 6) continue;
		cout << "h_motif " << ++index << ": " << h_motif_final[i] << endl;
	}

	double runtime = (double)(clock() - run_start) / CLOCKS_PER_SEC;

	cout << "\nHypergraph motif counting done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;
	cout << "-----------------------------------------" << endl << endl;

	node2hyperedge.clear();
	hyperedge2node.clear();
	hyperedge2node_set.clear();
	hyperedge_adj.clear();
	hyperedge_inter.clear();
	intersection.clear();
	upd_time.clear();
	sampled_idx.clear();
	intermediate_buckets.clear();
	final_buckets.clear();

	return 0;
}
