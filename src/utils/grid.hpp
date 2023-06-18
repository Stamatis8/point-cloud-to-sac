#ifndef MY_GRID_HPP
#define MY_GRID_HPP

#include <vector>

std::vector<std::vector<double>> grid(std::vector<std::vector<double>> domain, int Nx, int Ny){
    /*
        Description: create a uniform grid in [domain.at(0).at(0),domain.at(0).at(1)]x[domain.at(1).at(0),domain.at(1).at(1)]
            with Nx points in the horizontal direction and Ny points in the vertical

        Output:
            std::vector<std::vector<double>> pnts
                - pnts.at(i*Nx + j), i = 0,...,Ny-1, j = 0,...,Nx-1 is the point in the ith row
                    and jth column
    */
    std::vector<std::vector<double>> pnts = std::vector<std::vector<double>>(Nx*Ny, std::vector<double>(2,0.0));

    for(int i = 0; i < Ny; i++){
        for(int j = 0; j < Nx; j++){
            pnts.at(i*Nx+j) = {
                domain.at(0).at(0) + (domain.at(0).at(1)- domain.at(0).at(0))*double(j)/(Nx-1),
                domain.at(1).at(0) + (domain.at(1).at(1)- domain.at(1).at(0))*double(i)/(Ny-1),
            };
        }
    }

    return pnts;
}

#endif //MY_GRID_HPP