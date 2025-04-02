//
// Created by Aidan Hsu on 3/24/25.
//

#include "main.h"
#include "Grid.h"
#include "matplot/matplot.h"
#include "cmath"
#include <iostream>
using namespace std;

int main() {
    Grid CylGrid{82, 42, 2, 1};
    CylGrid.InitializeBoundaries();
    CylGrid.RunAlgebraicMethod();
    //CylGrid.printGrid();
    CylGrid.RunGaussSeidelDiffMethod();
    //CylGrid.displayGrid();
}
