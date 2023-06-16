#ifndef TRANSPOSEVECTOR_HPP
#define TRANSPOSEVECTOR_HPP

#include <vector>

template<typename data>
std::vector<std::vector<data>> TransposeVector(std::vector<std::vector<data>> M){
    /*
        Description: Returns the transpose of M. M.at(i).size() must be 
            equal to M.at(j).size() for all 0 <= i,j <= M.size()-1
        
        data requirements:
            data(0) must be defined
    */

    int n1 = M.size();
    int n2 = M.at(0).size();

    //initialize result
    std::vector<std::vector<data>> Mtranspose =
        std::vector<std::vector<data>>(n2,std::vector<data>(n1,data(0)));

    for(int i = 0; i < n1; i++){
        for(int j = 0; j < n2; j++){
            Mtranspose.at(j).at(i) = M.at(i).at(j);
        }
    }

    return Mtranspose;

}


#endif // TRANSPOSEVECTOR_HPP