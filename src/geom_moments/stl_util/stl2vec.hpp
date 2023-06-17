#ifndef STL2VEC_HPP
#define STL2VEC_HPP

#include <vector>
#include <string>

#include "stl_reader/stl_reader.h"

std::vector<std::vector<std::vector<double>>> stl2vec(std::string filename){
	/*
		Description: Converts .stl file format to a vector of triangles
		
		Input:
			- std::string filename
				must contain .stl file ending (ie "myfile.stl")
		Output:
			- std::vector<std::vector<std::vector<double>>> T
				T.at(i) contains the vertices T.at(i).at(0), T.at(i).at(1) and T.at(i).at(2)
	*/
	
	stl_reader::StlMesh <float, unsigned int> mesh (filename);

	//Initialize T
	std::vector<std::vector<std::vector<double>>> T(mesh.num_tris(),std::vector<std::vector<double>>(3,std::vector<double>(3,0)));

  	for(size_t itri = 0; itri < mesh.num_tris(); ++itri) {// triangles
      	for(size_t icorner = 0; icorner < 3; ++icorner) {// vertices
          	
          	const float* c = mesh.tri_corner_coords (itri, icorner);
          	
          	T.at(itri).at(icorner).at(0) = c[0];
          	T.at(itri).at(icorner).at(1) = c[1];
          	T.at(itri).at(icorner).at(2) = c[2];
      	}
  	}
  	
  	return T;
} 

#endif //STL2VEC_HPP
