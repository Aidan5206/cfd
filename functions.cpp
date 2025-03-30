//
// Created by Aidan Hsu on 3/29/25.
//

#include "functions.h"
#include <iostream>
#include <array>

using namespace std;

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
