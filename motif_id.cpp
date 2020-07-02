#include <iostream>

using namespace std;

int id_to_index[128] = {
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0,
	21, 23, 22, 24, 23, 25, 24, 26,
	0, 0, 0, 0, 0, 0, 0, 0,
	21, 22, 23, 24, 23, 24, 25, 26,
	21, 23, 23, 25, 22, 24, 24, 26,
	27, 28, 28, 29, 28, 29, 29, 30,
	1, 2, 2, 3, 2, 3, 3, 4,
	5, 6, 6, 8, 7, 9, 9, 10,
	5, 7, 6, 9, 6, 9, 8, 10,
	11, 13, 12, 14, 13, 15, 14, 16,
	5, 6, 7, 9, 6, 8, 9, 10,
	11, 12, 13, 14, 13, 14, 15, 16,
	11, 13, 13, 15, 12, 14, 14, 16,
	17, 18, 18, 19, 18, 19, 19, 20
};

int get_motif_index_new(int deg_a, int deg_b, int deg_c, int C_ab, int C_bc, int C_ca, int g_abc){
	int a = deg_a - (C_ab + C_ca) + g_abc;
	int b = deg_b - (C_bc + C_ab) + g_abc;
	int c = deg_c - (C_ca + C_bc) + g_abc;
	int d = C_ab - g_abc;
	int e = C_bc - g_abc;
	int f = C_ca - g_abc;
	int g = g_abc;
	int motif_id = (a > 0) + ((b > 0) << 1) + ((c > 0) << 2) + ((d > 0) << 3) + ((e > 0) << 4) + ((f > 0) << 5) + ((g > 0) << 6);
	return id_to_index[motif_id] - 1;
}
