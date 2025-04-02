//
// Created by Aidan Hsu on 3/25/25.
//

#ifndef GRID_H
#define GRID_H
#include "functions.h"
#include <vector>

using namespace std;

class Grid {
    // This class creates a new grid object
    public:
        int m_IL {};    // Rows
        int m_JL {};    //   Columns
        double m_domain_H {};   // Domain height
        double m_domain_L {};   // Domain length
        double m_cyl_radius {};
        double m_dx {};
        double m_dy {};

        GridData m_x_grid;
        GridData m_y_grid;

        explicit Grid(int IL = 42, int JL = 22, double H = 2, double L = 1, double r = 0.5);

        void InitializeBoundaries();
        void RunAlgebraicMethod();
        void RunGaussSeidelDiffMethod();
        void printGrid() const;
        void displayGrid() const;

};

#endif //GRID_H
