#ifndef SIMPLEMONTECARLO_HPP
#define SIMPLEMONTECARLO_HPP

#include <vector>
#include <cmath>

std::vector<std::vector<double>> SimpleMonteCarlo(
    std::vector<std::vector<double>> X,
    int N
);

std::vector<std::vector<double>> SimpleMonteCarlo(
    std::vector<std::vector<double>> X,
    int N
){
    /*
        Description: Randomly generates N samples in the design space X

        Input:

            - std::vector<std::vector<double>> X
                number of parameters = X.size()
                ith parameter \in [X.at(i).at(0), X.at(i).at(1)]

            - int N
                number of samples to generate
    */

    std::vector<std::vector<double>> S(N,std::vector<double>(X.size(),0));// Initializing final samples to all zeros

    for (int i = 0; i < N; i++){// iterating N times
        for (int j = 0; j < X.size(); j++){// iterating through sample coordinates
            S.at(i).at(j) = X.at(j).at(0) + (X.at(j).at(1)-X.at(j).at(0))*rand()/RAND_MAX;
        }
    }

    return S;
}


#endif //SIMPLEMONTECARLO_HPP