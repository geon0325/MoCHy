g++ -O3 -std=c++11 -lgomp -fopenmp main_exact_par.cpp -o exact_par;
./exact_par 4;
rm exact_par;
