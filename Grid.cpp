//
// Created by Aidan Hsu on 3/25/25.
//

#include "Grid.h"
#include "functions.h"
#include <iostream>
#include <chrono>
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
    , m_cyl_radius {r}
    , m_x_grid(IL, JL)
    , m_y_grid(IL, JL) {

    m_dx = (m_domain_L - m_cyl_radius) / (m_IL - 1); // Calc the grid spacing
    m_dy = (m_domain_H - m_cyl_radius) / (m_IL - 1);

    cout << "Grid Created! IL: " << m_IL << " JL: " << m_JL << " H: " << m_domain_H << " L: " << m_domain_L << " r: " << m_cyl_radius << endl;
}

void Grid::InitializeBoundaries() {
    // Find theta values for each cylinder grid point
    vector<double> grid_theta = linspace(M_PI, M_PI/2, m_JL);

    // Init ellipse arc coords to boundary edge (first row: i = 0)
    for (int j = 0; j < m_JL; j++) {
        m_x_grid.grid[0][j] = m_domain_L * (1 + cos(grid_theta[j])); // x = a + r_x * cos(theta)
        m_y_grid.grid[0][j] = m_domain_H * sin(grid_theta[j]);       // y = r_y * sin(theta)
    }

    // Init cylinder surface coords to boundary edge (last row: i = m_IL-1)
    for (size_t j = 0; j < grid_theta.size(); j++) {
        m_x_grid.grid[m_IL - 1][j] = m_cyl_radius * (2 + cos(grid_theta[j]));
        m_y_grid.grid[m_IL - 1][j] = m_cyl_radius * sin(grid_theta[j]);
    }

    // Init Symmetry (Horizontal Wall, first column: j = 0)
    vector<double> symmetry_grid = linspace(0, m_domain_L - m_cyl_radius, m_IL);
    for (int i = 0; i < m_IL; i++) {
        m_x_grid.grid[i][0] = symmetry_grid[i];
        m_y_grid.grid[i][0] = 0;
    }

    // Init Outlet (Vertical Wall, last column: j = m_JL-1)
    vector<double> outlet_grid = linspace(m_domain_H, m_cyl_radius, m_IL);
    for (int i = 0; i < m_IL; i++) {
        m_x_grid.grid[i][m_JL-1] = m_domain_L;
        m_y_grid.grid[i][m_JL-1] = outlet_grid[i];

    }
}

void Grid::RunAlgebraicMethod() {
    vector<double> grid_theta = linspace(M_PI, M_PI/2, m_JL);
    // Solve interior points (rows i=1 to i=m_IL-2)
    for (int i = 1; i < m_IL-1; i++) {
        for (int j = 1; j < m_JL-1; j++) {
            m_x_grid.grid[i][j] = m_domain_L + (m_domain_L - i * m_dx) * cos(grid_theta[j]);
            m_y_grid.grid[i][j] = (m_domain_H - i * m_dy) * sin(grid_theta[j]);
        }
    }
}

void Grid::RunGaussSeidelDiffMethod() {
    double norm1 = 1; double norm2 = 1;
    const double epsilon = 1e-5;  // Convergence criterion
    const int iter_max = 1864;         // Max iteration count
    int iter = 0;
    GridData x_grid_temp = m_x_grid;
    GridData y_grid_temp = m_y_grid;
    vector<double> diffs(x_grid_temp.data.size());

    auto start = chrono::high_resolution_clock::now();
    while (iter < iter_max && (norm1 > epsilon || norm2 > epsilon)) {
        iter ++;

        // 🔄 Copy current state into temp buffers
        x_grid_temp.data = m_x_grid.data;
        x_grid_temp.rebuildGridPointers();

        y_grid_temp.data = m_y_grid.data;
        y_grid_temp.rebuildGridPointers();

        for (int i = 1; i < m_IL-1; i++) {
            for (int j = 1; j < m_JL-1; j++) {
                double alpha =  (0.5*(x_grid_temp.grid[i][j+1] - x_grid_temp.grid[i][j-1])) * (0.5*(x_grid_temp.grid[i][j+1] - x_grid_temp.grid[i][j-1])) +
                                (0.5*(y_grid_temp.grid[i][j+1] - y_grid_temp.grid[i][j-1])) * (0.5*(y_grid_temp.grid[i][j+1] - y_grid_temp.grid[i][j-1]));

                double beta =   0.5*(x_grid_temp.grid[i+1][j] - x_grid_temp.grid[i-1][j]) *
                                0.5*(x_grid_temp.grid[i][j+1] - x_grid_temp.grid[i][j-1]) +
                                0.5*(y_grid_temp.grid[i+1][j] - y_grid_temp.grid[i-1][j]) *
                                0.5*(y_grid_temp.grid[i][j+1] - y_grid_temp.grid[i][j-1]);

                double gamma =  (0.5*(x_grid_temp.grid[i+1][j] - x_grid_temp.grid[i-1][j])) * (0.5*(x_grid_temp.grid[i+1][j] - x_grid_temp.grid[i-1][j])) +
                                (0.5*(y_grid_temp.grid[i+1][j] - y_grid_temp.grid[i-1][j])) * (0.5*(y_grid_temp.grid[i+1][j] - y_grid_temp.grid[i-1][j]));

                double RHS1 =   -alpha * (x_grid_temp.grid[i+1][j] + x_grid_temp.grid[i-1][j]) +
                                0.5*beta * (x_grid_temp.grid[i+1][j+1] - x_grid_temp.grid[i-1][j+1] - x_grid_temp.grid[i+1][j-1] + x_grid_temp.grid[i-1][j-1]) -
                                gamma * (x_grid_temp.grid[i][j+1] + x_grid_temp.grid[i][j-1]);

                double RHS2 =   -alpha * (y_grid_temp.grid[i+1][j] + y_grid_temp.grid[i-1][j]) +
                                0.5*beta * (y_grid_temp.grid[i+1][j+1] - y_grid_temp.grid[i-1][j+1] - y_grid_temp.grid[i+1][j-1] + y_grid_temp.grid[i-1][j-1]) -
                                gamma * (y_grid_temp.grid[i][j+1] + y_grid_temp.grid[i][j-1]);

                x_grid_temp.grid[i][j] = -RHS1 / (2 * (alpha+gamma));
                y_grid_temp.grid[i][j] = -RHS2 / (2 * (alpha+gamma));
            }
        }
        // Update the norms
        norm1 = maxAbsDifference(x_grid_temp, m_x_grid, diffs);
        norm2 = maxAbsDifference(y_grid_temp, m_y_grid, diffs);
        // Update the matrices
        m_x_grid.swapWith(x_grid_temp);
        m_y_grid.swapWith(y_grid_temp);
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    cout << "Gauss Seidel finished in " << iter << " iterations " << elapsed.count() << " sec" << endl;
    cout << "Norm 1: " << norm1 << " Norm 2: " << norm2 << endl;

}

void Grid::printGrid() const {
    cout << "X Grid:" << endl;
    for (int i = 0; i < m_IL; i++) {
        for (int j = 0; j < m_JL; j++) {
            cout << m_x_grid.grid[i][j] << " ";
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
            x_flat[i * m_JL + j] = m_x_grid.grid[i][j];
            y_flat[i * m_JL + j] = m_y_grid.grid[i][j];
        }
    }

    auto f = figure();
    scatter(x_flat, y_flat, 1);
    axis(matplot::equal);
    show();
}
