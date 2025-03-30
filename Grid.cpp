//
// Created by Aidan Hsu on 3/25/25.
//

#include "Grid.h"
#include "functions.h"
#include <iostream>
#include <array>
#include <cmath>
#include <matplot/freestanding/plot.h>
#include <matplot/util/common.h>
using namespace std;
using namespace matplot;


Grid::Grid(const int IL, const int JL, const double H, const double L, const double r)
    : m_IL {IL}
    , m_JL {JL}
    , m_domain_H {H}
    , m_domain_L {L}
    , m_cyl_radius {r}  {

    // Allocate a single contiguous block of memory
    m_x_data = new double[m_IL * m_JL]{0.0};  // Initialize to 0
    m_y_data = new double[m_IL * m_JL]{0.0};

    // Allocate the pointers for [i][j] indexing
    m_x_grid = new double *[m_IL];
    m_y_grid = new double *[m_IL];
    for (int i = 0; i < m_IL; i++) {
        m_x_grid[i] = &m_x_data[i * m_JL];  // Set row pointers
        m_y_grid[i] = &m_y_data[i * m_JL];
    }

    cout << "Grid Created! IL: " << m_IL << " JL: " << m_JL << " H: " << m_domain_H << " L: " << m_domain_L << " r: " << m_cyl_radius << endl;
}

Grid::~Grid() {
    delete[] m_x_data;
    delete m_x_grid;
    delete[] m_y_data;
    delete m_y_grid;
}

void Grid::Initialize() {
    m_dx = (m_domain_L - m_cyl_radius) / (m_IL - 1); // Calc the grid spacing
    m_dy = (m_domain_H - m_cyl_radius) / (m_IL - 1);

    // Find theta values for each cylinder grid point
    vector<double> cyl_grid_theta = linspace(M_PI, M_PI/2, m_JL);

    // Init cylinder surface coords to boundary edge (i = m_IL - 1)
    for (size_t j = 0; j < cyl_grid_theta.size(); j++) {
        m_x_grid[m_IL - 1][j] = m_cyl_radius * (2 + cos(cyl_grid_theta[j]));
        m_y_grid[m_IL - 1][j] = m_cyl_radius * sin(cyl_grid_theta[j]);
    }

    // Init ellipse arc coords to boundary edge (i=0)
    vector<double> ellipse_bc_grid_theta = linspace(M_PI, M_PI/2, m_JL); // pi (left) to pi/2 (top)
    for (int j = 0; j < m_JL; j++) {
        m_x_grid[0][j] = m_domain_L * (1 + cos(ellipse_bc_grid_theta[j]));
        m_y_grid[0][j] = m_domain_H * sin(ellipse_bc_grid_theta[j]);
    }

    // Init Outlet (Vertical Wall j = JL - 1)
    vector<double> outlet_grid = linspace(m_domain_H, m_cyl_radius, m_IL);
    for (int i = 0; i < m_IL; i++) {
        m_x_grid[i][m_JL-1] = m_domain_L;
        m_y_grid[i][m_JL-1] = outlet_grid[i];

    }

    // Init Symmetry (Horizontal Wall j = 0)
    vector<double> symmetry_grid = linspace(0, m_domain_L - m_cyl_radius, m_IL);
    for (int i = 0; i < m_IL; i++) {
        m_x_grid[i][0] = symmetry_grid[i];
        m_y_grid[i][0] = 0;
    }

}

void Grid::printGrid() const {
    cout << "X Grid:" << endl;
    for (int i = 0; i < m_IL; i++) {
        for (int j = 0; j < m_JL; j++) {
            cout << m_x_grid[i][j] << " ";
        }
        cout << endl;
    }
}

void Grid::displayGrid() const {
    // Convert data to 1D vectors
    vector<double> x_flat(m_IL * m_JL);
    vector<double> y_flat(m_IL * m_JL);

    // Copy data
    for (int i = 0; i < m_IL; i++) {
        for (int j = 0; j < m_JL; j++) {
            x_flat[i * m_JL + j] = m_x_grid[i][j];
            y_flat[i * m_JL + j] = m_y_grid[i][j];
        }
    }
    scatter(x_flat, y_flat);
    axis(matplot::equal);
    show();
}
