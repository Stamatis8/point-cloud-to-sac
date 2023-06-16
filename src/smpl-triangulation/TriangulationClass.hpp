#ifndef TRIANGULATION_CLASS_HPP
#define TRIANGULATION_CLASS_HPP

#include <vector>
#include <cmath>

namespace smpl_triangulation{
	
	template<class vertex>
	class Triangulation
	{
	public:
		/*
			vertex class requirements:
				- must have default constructor vertex()
		*/
		
		std::vector<vertex> nodes;
		//  a list of all the vertices
		
		std::vector<std::vector<int>> faces;
		/* 
			- A list with all the faces of the triangulation
			- The ith face is structured as follows:
				- faces.at(i).at(0) = index of 1st vertex, from nodes
				- faces.at(i).at(1) = index of 2nd vertex, from nodes
				- faces.at(i).at(2) = index of 3rd vertex, from nodes
				- faces.at(i).at(3) = index of neighboring face opposite 1st node, from faces 
				- faces.at(i).at(4) = index of neighboring face opposite 2nd node, from faces 
				- faces.at(i).at(5) = index of neighboring face opposite 3rd node, from faces 
			- If there is no neighbour opposite some vertex (ie at a corner of the triangulation) then the 
				index of the current face is placed there
				
			example:
				faces.at(0) = {5,4,3,2,1,0}
					- the first face
					- its first vertex is nodes.at(5)
					- its second vertex is nodes.at(4)
					- its third vertex is nodes.at(3)
					- its neighbour opposite its first vertex is faces.at(2)
					- its neighbour opposite its second vertex is faces.at(1)
					- its neighbour opposite its third vertex is faces.at(0) (itself, so it is at the edge of the
																			 convex hull)
		*/
		
		bool is_oriented = false;
		
		/* Constructors */
		
		Triangulation(){}
		
		Triangulation(std::vector<vertex> nodes_in, std::vector<std::vector<int>> faces_in){
			
			/* Todo: Check for validity of triangulation input */
			
			this->nodes = nodes_in;
			this->faces = faces_in;
		}
		
		/* Methods */
		
		operator std::vector<std::vector<vertex>>() {
			/*
				Description: Converts this triangulation data structure into a simpler form, removing connectivity information.
				
				- Output:
					- std::vector<std::vector<vertex>> T
						- T.at(i) represents the ith face
						- T.at(i).at(j) for j = 0,1,2 are the 1st 2nd and 3rd vertices of the ith face
							respectively
			*/
			std::vector<std::vector<vertex>> T(this->faces.size(),std::vector<vertex>(3,vertex()));
			
			for (int i = 0; i < this->faces.size(); i++){
				T.at(i) = this->get_vertices(i);
			}
			
			return T;
		}
		
		std::vector<vertex> get_vertices(int i){
			/*
				Description: Returns the vertices of the ith face
				
				Input:
					- int i
						0 <= i <= (this->faces.size() - 1)
				Output:
					- std::vector<vertex>
						a vector of length 3, each element a vertex
			*/
			
			return std::vector<vertex> {
				this->nodes.at(this->faces.at(i).at(0)),
				this->nodes.at(this->faces.at(i).at(1)),
				this->nodes.at(this->faces.at(i).at(2))
				};
		}
		
		bool orient(int i = 0){
			/*
				Description: Orients the triangulation so that it is in alignment with the orientation
					of ith triangle and sets is_oriented attribute to true. If the triangluation is 
					non-orientable false is returned and is_oriented attribute is set to false
					
					To orient, we start with the ith face and check if it agrees in orientation with its neighbours. If it does
					not, we use flip() to make it so. We note down the ith face. We go to each of its neighbours. For each, we note it
					down. We check if it agrees with its neighbours. If it does not we only use flip() if that neighbour is not noted down. If 
					it is, we have a non-orientable triangulation and return false. We continue until we have gone through every face
				
				Input:
					- int i
						triangulation is attempted to be oriented according to faces.at(i) orientation
						0 <= i <= (faces.size()-1)
				
				Output:
					- bool this->is_oriented
						if triangulation is non-orientable, false is returned. Else, true
			*/
			
			int n = this->faces.size();// num of faces
			int current_face;// index of current face
			std::vector<int> list = {i};// list of faces to check neighbours in current loop
			std::vector<int> next_list;// list of faces to check neighbours in next loop 
			std::vector<int> checked = {i};// list of faces that have correct orientation
			std::vector<int> processed;// list of faces that have had their neighbours checked
			
			while (checked.size() < n){// stop when we have checked entire triangulation
			
				next_list.clear();
				
				for (int f = 0; f < list.size(); f++){// check all faces in list
				
					current_face = list.at(f);// index of current face in this->faces
					processed.push_back(current_face);
					
					for (int neighbour = 3; neighbour < 6; neighbour++){
					
						/* Checking neighbour */
					
						if (this->is_in(this->faces.at(current_face).at(neighbour),checked)){
							// neighbour of current_face has been checked
						
							if (!(this->agree(current_face,neighbour))){// current_face and neighbour have the correct orientation
																		// so triangulation is non-orientable
								return false;
							}	
						}
						else{
							// neighbour of current_face has not been checked
						
							if (!(this->agree(current_face,neighbour))){// current_face has the correct orientation
																		// if its neighbour has the wrong one, correct it
								this->flip(this->faces.at(current_face).at(neighbour));
							}
							
							checked.push_back(this->faces.at(current_face).at(neighbour));
						}
						
						
						/* Checking wether to send neighbour in next list for processing */
						
						if (!this->is_in(this->faces.at(current_face).at(neighbour),processed)){
							//all neighbours of faces in list that have not been processed go to next list
							
							next_list.push_back(this->faces.at(current_face).at(neighbour));
						}
						
					}// neighbours
					
				}// checking list
				
				list = next_list;
			}// while
			
			return true;
			
		}
	private:
		bool agree(int i, int j){
			/*
				Description: Checks if ith triangles orientation agrees with its jth **neighbour**
				
				Input:
					- int i
						0 <= i <= (this->faces.size()-1)
					- int j
						3 <= j <= 5 (indicates the a neighbour of j. See definition of this->faces)
				Output:
					- bool
						true if the orientations of faces.at(i) and faces.at(faces.at(i).at(j)) agree
			*/
			
			/* Check if neighbour is itself (see this->faces definition) */
			
			if (i == this->faces.at(i).at(j)) return true;// face has the same orientation with itself
			
			int a = this->faces.at(i).at((j-2)%3);// index of first common vertex, from nodes
			int b = this->faces.at(i).at((j-1)%3);// index of second common vertex, from nodes
			
			int ia = (j-2)%3;// index of a at ith face
			int ib = (j-1)%3;// index of b at ith face
			
			/* Find ja,jb */
			
			int ja;// index of a at jth face
			int jb;// index of b at jth face
			
			for (int k = 0; k < 3; k++){
				
				if (a == faces.at(faces.at(i).at(j)).at(k)) ja = k;
				
				if (b == faces.at(faces.at(i).at(j)).at(k)) jb = k;
				
			}
			
			/* Compare ia,ib,ja,jb */
			
			if (((ia + 1)%3 == ib && (ja + 1)%3 == jb)
				||((ia - 1)%3 == ib && (ja - 1)%3 == jb)){// Two neighboring triangles have the same orientation iff their common side
									 					  // is oriented differently between them
				// abc is the order of the vertices in both triangles. Therefore they have different orientation
				
				return false;
			}
			else {
			
				// abc and acb is the order of the two triangles. Therefore they have the same orientation
			
				return true;
			}
		}// agree()
		
		void flip(int k, int i = 0, int j = 1){
			/*
				Description: flips the ith vertex with the jth vertex of the kth face
				
				Input:
					- int k
						0 <= k <= (this->faces.size()-1)
					- int i
						0 <= i <= 2
					- int j
						0 <= j <= 2
			*/
			
			if (i != j){// if we are trying to flip the same vertex, do nothing
				
				// Switching neighboring faces
				int previous_jth_face = this->faces.at(k).at(j+3);				
				this->faces.at(k).at(j+3) = this->faces.at(k).at(i+3);
				this->faces.at(k).at(i+3) = previous_jth_face;
				
				// Switching vertices
				int previous_jth_vertex = this->faces.at(k).at(j);
				this->faces.at(k).at(j) = this->faces.at(k).at(i);
				this->faces.at(k).at(i) = previous_jth_vertex;
			
			}
		}// flip()
		
		bool is_in(int i, std::vector<int> v){
			for (int j = 0;j < v.size();j++){
				if (v.at(j) == i) return true;
			}
			return false;
		}
	};// class
	
}// namespace

#endif// TRIANGULATION_CLASS_HPP
