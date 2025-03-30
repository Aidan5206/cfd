//
// Created by Aidan Hsu on 3/25/25.
//

#ifndef GRID_H
#define GRID_H

class Grid {
    // This class creates a new grid object
    public:
        int m_IL {};
        int m_JL {};
        double m_domain_H {};   // Domain height
        double m_domain_L {};   // Domain length
        double m_cyl_radius {};
        double m_dx {};
        double m_dy {};
        double **m_x_grid;  // Pointer-to-pointer for [i][j] indexing
        double *m_x_data;   // Contiguous 1D block of memory
        double **m_y_grid;
        double *m_y_data;


        explicit Grid(int IL = 42, int JL = 22, double H = 2, double L = 1, double r = 0.5);
        ~Grid();

        void Initialize();
        void printGrid() const;
        void displayGrid() const;

};

#endif //GRID_H
