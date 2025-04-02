//
// Created by Aidan Hsu on 3/29/25.
//

#include "functions.h"
#include <iostream>
#include <array>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using U32 = uint32_t;

vector<double> linspace(double start, double end, int N) {
    vector<double> result(N);
    if (N==1) {
        result[0] = start;
        return result;
    }
    double step = (end - start) / (N - 1);
    for (size_t i = 0; i < N; i++) {
        result[i] = start + i * step;
    }

    return result;
}


GridData::GridData(int r, int c) : data(r * c, 0.0), rows(r), cols(c) {
    grid.resize(r);
    for (int i = 0; i < r; i++) {
        grid[i] = &data[i * c]; // Set row pointers
    }
}
// Copy Constructor with Deep Copy
GridData::GridData(const GridData &other) : rows(other.rows), cols(other.cols), data(other.data) {
    grid.resize(rows);
    for (int i = 0; i < rows; i++) {
        grid[i] = &data[i * cols];
    }
}
// Copy Assignment Operator with Deep Copy
GridData& GridData::operator=(const GridData &other) {
    if (this == &other) return *this;

    rows = other.rows;
    cols = other.cols;
    data = other.data;

    grid.resize(rows);
    for (int i = 0; i < rows; i++) {
        grid[i] = &data[i * cols];
    }

    return *this;
}

void GridData::rebuildGridPointers() {
    for (int i = 0; i < rows; i++) {
        grid[i] = &data[i * cols];
    }
}

void GridData::swapWith(GridData& other) {
    swap(data, other.data);
    swap(grid, other.grid);
    swap(rows, other.rows);
    swap(cols, other.cols);
    rebuildGridPointers();
    other.rebuildGridPointers();
}



double maxAbsDifference(const GridData& grid1, const GridData& grid2, vector<double>& diffs) {
    transform(  grid1.data.begin(), grid1.data.end(),
                grid2.data.begin(), diffs.begin(),
                [](double a, double b) { return fabs(a-b); });

    //return *std::max_element(diffs.begin(), diffs.end());
    return
}
