#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include "read_data.cpp"
#include "motif_id.cpp"

#include <omp.h>
#include <atomic>
//#include <tbb/concurrent_unordered_map.h>

using namespace std;

struct hyperwedge{
	int a, b, C_ab;	
};

inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
	chrono::system_clock::time_point start;
	chrono::system_clock::time_point run_start;
	chrono::duration<double> dur;
	int progress;

	int sampling_size = stoi(argv[1]);
	int num_threads = stoi(argv[2]);

	omp_set_num_threads(num_threads);

	string graphFile = "dblp_graph.txt";

	cout << "Sampling size: " << sampling_size << endl << endl;

	// Read data
	start = chrono::system_clock::now();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);

	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
	dur = std::chrono::system_clock::now() - start;
	cout << "Reading data done: "
		<< dur.count() << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

	// Adjacency list construction
	chrono::system_clock::time_point ss = chrono::system_clock::now();
	start = chrono::system_clock::now();
	run_start = chrono::system_clock::now();
	hyperedge2node.resize(E); hyperedge2node_set.resize(E);
	vector< vector< pair<int, int> > > hyperedge_adj;
	vector< unordered_map<int, int> > hyperedge_inter;
	hyperedge_adj.resize(E);
	hyperedge_inter.resize(E);
	vector< vector<long long> > upd_time(num_threads);
	vector< hyperwedge > W;
	vector< long long > W_cnt(E + 1, 0LL);

	for(int i=0;i<num_threads;i++){
		upd_time[i].resize(E, -1LL);
	}
	#pragma omp parallel for
	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		int tid = omp_get_thread_num();
		int deg_a = 0;
		long long l_hyperedge_a = (long long)hyperedge_a;
		for (const int &node: hyperedge2node[hyperedge_a]){
			for (const int &hyperedge_b: node2hyperedge[node]){
				if (hyperedge_b == hyperedge_a) continue;
				if ((upd_time[tid][hyperedge_b] >> 31) ^ hyperedge_a){
					upd_time[tid][hyperedge_b] = (l_hyperedge_a << 31) + deg_a; deg_a++;
					hyperedge_adj[hyperedge_a].push_back({hyperedge_b, 0});
					if (hyperedge_a < hyperedge_b) W_cnt[hyperedge_a + 1]++;
				}
				hyperedge_adj[hyperedge_a][(int)(upd_time[tid][hyperedge_b] & 0x7FFFFFFFLL)].second++;
			}
		}
	}
	for (int i = 0; i < E; i++) W_cnt[i + 1] += W_cnt[i];
	W.resize(W_cnt[E]);

	#pragma omp parallel for
	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		int deg_a = hyperedge_adj[hyperedge_a].size();
		hyperedge_inter[hyperedge_a].rehash(deg_a);
		int _cnt = 0;
		for (int i = 0; i < deg_a; i++){
			int hyperedge_b = hyperedge_adj[hyperedge_a][i].first;
			int C_ab = hyperedge_adj[hyperedge_a][i].second;
			hyperedge_inter[hyperedge_a].insert({hyperedge_b, C_ab});
			if (hyperedge_a < hyperedge_b) {
				W[W_cnt[hyperedge_a] + _cnt] = hyperwedge{hyperedge_a, hyperedge_b, C_ab};
				_cnt++;
			}
		}
	}
	
	dur = std::chrono::system_clock::now() - start;
	cout << "# of hyperwedges: " << W.size() << "\n";
	cout << "Adjacency list construction done: "
		<< dur.count() << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;


	// h_motif counting via hyperwedge smapling
	start = chrono::system_clock::now();
	vector< vector<long long> > h_motif(num_threads);
	vector< vector<int> > intersection(num_threads);

	mt19937 gen[num_threads];
	uniform_int_distribution<> dist[num_threads];
	for(int i=0;i<num_threads;i++){
		gen[i] = mt19937(2020 * (i + 1));
		dist[i] = uniform_int_distribution<>(0, ((int)W.size())-1-1);
		std::fill(upd_time[i].begin(), upd_time[i].end(), -1LL);
		intersection[i].resize(V);
		h_motif[i].resize(30, 0);
	}
	vector< std::atomic_flag > mutex(E);
	#pragma omp parallel for
	for (int sample = 0; sample < sampling_size; sample++){
		int tid = omp_get_thread_num();

		int sample_index = dist[tid](gen[tid]);
		int hyperedge_a = W[sample_index].a;
		int hyperedge_b = W[sample_index].b;
		int C_ab = W[sample_index].C_ab;
		
		int size_a = (int)hyperedge2node[hyperedge_a].size();
		int size_b = (int)hyperedge2node[hyperedge_b].size();
		int deg_a = (int)hyperedge_adj[hyperedge_a].size();
		int deg_b = (int)hyperedge_adj[hyperedge_b].size();

		upd_time[tid][hyperedge_a] = upd_time[tid][hyperedge_b] = sample;

		int min_ab = hyperedge_a, max_ab = hyperedge_b;
		if (size_a > size_b) min_ab = hyperedge_b, max_ab = hyperedge_a;
		const auto &nodes = hyperedge2node_set[max_ab]; auto it_end = nodes.end(); int cnt = 0;
		for (const int &node: hyperedge2node[min_ab]){ if(nodes.find(node) != it_end) intersection[tid][cnt++] = node;}

		for (int i = 0; i < deg_b; i++){
			int hyperedge_c = hyperedge_adj[hyperedge_b][i].first, C_bc = hyperedge_adj[hyperedge_b][i].second;
			if (upd_time[tid][hyperedge_c] ^ sample){
				upd_time[tid][hyperedge_c] = sample;

				int size_c = (int)hyperedge2node[hyperedge_c].size();
				int C_ca = 0, g_abc = 0;
				while(mutex[hyperedge_a].test_and_set(std::memory_order_acquire));
				C_ca = hyperedge_inter[hyperedge_a][hyperedge_c];
				mutex[hyperedge_a].clear(std::memory_order_release);
				const auto &nodes = hyperedge2node_set[hyperedge_c]; auto it_end = nodes.end();
				for (int k = 0; k < C_ab; k++){ if(nodes.find(intersection[tid][k]) != it_end) g_abc++; }
				
				int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
				h_motif[tid][motif_index]++;
			}
		}

		for (int i = 0; i < deg_a; i++){
			int hyperedge_c = hyperedge_adj[hyperedge_a][i].first, C_ca = hyperedge_adj[hyperedge_a][i].second;
			if (upd_time[tid][hyperedge_c] ^ sample){
				upd_time[tid][hyperedge_c] = sample;
				
				int size_c = (int)hyperedge2node[hyperedge_c].size();
				int C_bc = 0, g_abc = 0;

				int motif_index = get_motif_index_new(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
				h_motif[tid][motif_index]++;
			}
		}
	}

	int index = 0;
	vector<long double> h_motif_final(30, 0);
	for (int i = 0; i < 30; i++){
		for(int j=0;j<num_threads;j++) h_motif_final[i] += h_motif[j][i];
		h_motif_final[i] *= (long double)W.size() / sampling_size;
		if (20 <= i && i <= 25)
			h_motif_final[i] /= 2.0;
		else
			h_motif_final[i] /= 3.0;
		if (i == 0 || i == 1 || i == 4 || i == 6) continue;
		cout << "h-motif " << ++index << ": " << h_motif_final[i] << endl;
	}

	dur = std::chrono::system_clock::now() - run_start;
	double runtime = (double)dur.count();

	dur = std::chrono::system_clock::now() - start;
	cout << "\nHypergraph motif counting done: "
		<< dur.count() << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;
	cout << "-----------------------------------------" << endl << endl;

	W.clear();
	node2hyperedge.clear();
	hyperedge2node.clear();
	hyperedge2node_set.clear();
	h_motif.clear();
	intersection.clear();
	hyperedge_adj.clear();
	hyperedge_inter.clear();
	upd_time.clear();

	return 0;
}
