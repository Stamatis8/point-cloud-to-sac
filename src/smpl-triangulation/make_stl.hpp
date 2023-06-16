#ifndef MAKESTL_HPP
#define MAKESTL_HPP

#include <vector>
#include <fstream>

#include "TriangulationClass.hpp"

namespace smpl_triangulation{

	void make_stl(Triangulation<std::vector<double>> t,std::string filename, std::string solidname);
	
	void make_stl(Triangulation<std::vector<double>> t, std::string filename, std::string solidname = ""){
				/*
					Description: Constructs .stl ASCII format in filename from triangulation. If triangulation is 2D, an extra
						coordinate equal to zero is added to the vertices
					
					Input:
						- Triangulation<std::vector<double>> t
							t.faces.size() > 0
					
						- std::string filename
							data is saved in "filename", so .stl or other specifier must be included
						
						- std::string solidname
							name of solid at the top of the stl file
				
					Note: the triangulation is not necessarily oriented, so the unit normals might not belong to the
						same orientation. Make sure to call the orient() method before the make_stl() method if this behaviour
						is not desired
				*/
				
				auto exterior = [] (std::vector<double> a, std::vector<double> b){// exterior product lamda expression
				
					std::vector<double> n(3,0);
					n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1)); 
					n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2));
					n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0));
					double n_mag = std::sqrt(n.at(0)*n.at(0) + n.at(1)*n.at(1) + n.at(2)*n.at(2));
					n.at(0) /= n_mag;
					n.at(1) /= n_mag;
					n.at(2) /= n_mag;
					
					return n;
				};
				
				auto subtract = [] (std::vector<double> a, std::vector<double> b){// lamda subtraction expression
					std::vector<double> result(a.size(),0);
					for (int i = 0; i < a.size(); i++){
						result.at(i) = a.at(i) - b.at(i);
					}
					return result;
				};
				
				/* If 2D converting to 3D */
				
				if (t.nodes.at(0).size() == 2){
					std::vector<std::vector<double>> new_nodes(t.nodes.size(),std::vector<double>(3,0));
					for (int i = 0; i< t.nodes.size(); i++){
						new_nodes.at(i).at(0) = t.nodes.at(i).at(0);
						new_nodes.at(i).at(1) = t.nodes.at(i).at(1);
						new_nodes.at(i).at(2) = 0;
					}
					t.nodes = new_nodes;
				}
				
				/* Writting to file */
				
				std::fstream file;
				file.open(filename,std::ios::out);
				
				std::vector<double> normal;
				std::vector<std::vector<double>> vertices;
				
				file << "solid " << solidname << std::endl;
				for (int i = 0; i < t.faces.size(); i++){
				
					vertices = t.get_vertices(i);
					
					normal = exterior(subtract(vertices.at(2),vertices.at(0)),subtract(vertices.at(1),vertices.at(0)));
					
					file << "	facet normal " << normal.at(0) << " " << normal.at(1) << " " << normal.at(2) << std::endl;
					file << "		outer loop" << std::endl;
						for (int j = 0; j < 3; j++){
							file << "			vertex " << vertices.at(j).at(0) << " " << vertices.at(j).at(1) << " " << vertices.at(j).at(2) << std::endl;
						}
					file << "		endloop" << std::endl;
					file << "	endfacet" << std::endl;
				}
				file << "endsolid";
				file.close();	
		}// make_stl()
}// namespace

#endif // MAKESTL_HPP
