g++ -O3 -std=c++11 -lgomp -fopenmp main_approx_ver2_par.cpp -o approx_ver2_par;
./approx_ver2_par 10000 4;
rm approx_ver2_par;
