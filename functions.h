//
// Created by Aidan Hsu on 3/29/25.
//

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <array>
#include <cstdint>
#include <vector>
#include <algorithm>

using namespace std;
using U32 = uint32_t;

vector<double> linspace(double start, double end, int N);
struct GridData {
    vector<double> data;
    vector<double*> grid;
    int rows, cols;
    GridData(int r, int c);
    GridData(const GridData& other);
    GridData& operator=(const GridData& other);
    void rebuildGridPointers();
    void swapWith(GridData& other);
};
double maxAbsDifference(const GridData& grid1, const GridData& grid2, vector<double>& diffs);

#endif //FUNCTIONS_H
