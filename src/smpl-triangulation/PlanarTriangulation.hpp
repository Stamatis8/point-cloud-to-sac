#ifndef PLANARTRIANGULATION_HPP
#define PLANARTRIANGULATION_HPP

#include <vector>
#include <cmath>

#include "TriangulationClass.hpp"

namespace smpl_triangulation{

Triangulation<std::vector<double>> PlanarTriangulation(
		std::vector<double> X,
		std::vector<double> Y,
		int N);

Triangulation<std::vector<double>> PlanarTriangulationStitchTopBottom(
		std::vector<double> X,
		std::vector<double> Y,
		int N);

Triangulation<std::vector<double>> PlanarTriangulationStitchLeftRight(
		std::vector<double> X,
		std::vector<double> Y,
		int N);

Triangulation<std::vector<double>> PlanarTriangulation(
		std::vector<double> X,
		std::vector<double> Y,
		int N)
	{
		/*
			Description: Given the rectangle X \times Y \in R^2, it is subdivided into approximately N triangles
				which span the entire rectangle. First, each interval X,Y is subdivided into M = ceiling(sqrt(N/2)) 
				subintervals. These subintervals partition the rectangle into M^2 subrectangles. Each subrectangle 
				is divided into two triangles by joining opposite vertices. 
				
				Conclusively, given N, a triangulation of X \times Y into 2*ceil(sqrt(N/2))^2 triangles is produced
			
			Note: The output triangulation is oriented
			
			Input:
			
				- std::vector<double> X
					First side of rectangle: [X.at(0), X.at(1)]
				- std::vector<double> Y
					Second side of rectangle: [Y.at(0), Y.at(1)]
				- int N
					N > 0
					
			
			Output:
			
				- std::vector<std::vector<std::vector<double>>> triangles
					triangles.at(i) is the ith triangle consisting of vertices triangles.at(i).at(0), triangles.at(i).at(1),
					triangles.at(i).at(2)
					
					A triangulation consisting of 2*ceil(sqrt(N/2))^2 rectangles is produced
			
		*/
		
		/* Subdividing X,Y */
		
		int M = std::ceil(std::sqrt(N/2)); //number of subintervals in X,Y
		
		std::vector<std::vector<double>> nodes((M+1)*(M+1), std::vector<double>(2,0));// nodes for Triangulation class
		// nodes.at(k*(M+1)) through nodes.at(k*(M+1) + M) for k = 0,...,M are the points in the kth row of the planar discretization
		
		int count = 0;// counts nodes
		
		for (int i = 0; i < (M+1); i++){
			for (int j = 0; j < (M+1); j++){
				nodes.at(count).at(0) = X.at(0) + (X.at(1)-X.at(0))*j/M;
				nodes.at(count).at(1) = Y.at(0) + (Y.at(1)-Y.at(0))*i/M;
				count++;
			}
		}
		
		/* Constructing Triangles */
		
		std::vector<std::vector<int>> faces(2*M*M,std::vector<int>(6,0));// faces for Triangulation class
		// faces.at(k*2*M) through faces.at(k*2*M + 2*M - 1) for k = 0,...,(M-1) are the triangles in the kth row of subrectangles
		// faces.at(k*2*M + 2*j) for j = 0,...,M is the lower right triangle of the jth subrectangle in the kth row of subrectangles
		// faces.at(k*2*M + 2*j + 1) for j = 0,...,M is the upper left triangle of the jth subrectangle in the kth row of subrectangles
				
		//Indexes of triangles in faces, relative to  current triangle
		int down_upper_left = 0;// upper left triangle in subrectangle down
		int top_lower_right = 0;// etc.
		int right_upper_left = 0;
		int left_lower_right = 0;
		
		count = 0;// counts faces
		
		for (int i = 0; i < M; i++){
			for (int j = 0; j < M; j++){
				// (i,j) denotes the (i,j)th discrete subrectangle of X \times Y
				
				down_upper_left = (i-1)*2*M + 2*j + 1;
				top_lower_right = (i+1)*2*M + 2*j;
				right_upper_left = i*2*M + 2*(j+1) + 1;
				left_lower_right = i*2*M + 2*(j-1);
				
				/* lower right triangle */
				
				faces.at(count).at(2) = i*(M+1) + j;// lower left vertex
				faces.at(count).at(1) = i*(M+1) + j + 1;// lower right vertex
				faces.at(count).at(0) = (i+1)*(M+1) + j + 1;// upper right vertex
				
				faces.at(count).at(4) = count + 1;// oppposite triangle to low right corner is the other half of the subrectangle
				
				// Determining neighbours (see Triangulation class definition)
				if (i == 0 && j != (M-1)) {// bottom side but not right corner of XY
					faces.at(count).at(5) = right_upper_left;
					faces.at(count).at(3) = count;// Triang. class convention (points to itself)
				}
				else if (i == 0 && j == (M-1)){// bottom right corner of XY
					faces.at(count).at(5) = count;
					faces.at(count).at(3) = count;
				}
				else if (i == (M-1) && j == (M-1)){// upper right corner of XY
					faces.at(count).at(5) = count;
					faces.at(count).at(3) = down_upper_left;
				}
				else{// has subrectangles to the right and down
					faces.at(count).at(5) = right_upper_left;
					faces.at(count).at(3) = down_upper_left;
				}
				
				count++;
				
				/* upper left triangle */
				
				faces.at(count).at(0) = i*(M+1) + j;// lower left vertex
				faces.at(count).at(1) = (i+1)*(M+1) + j;// upper left vertex
				faces.at(count).at(2) = (i+1)*(M+1) + j + 1;// upper right vertex
				
				faces.at(count).at(4) = count - 1;// oppposite triangle to upper left corner is the other half of the subrectangle
						
				// Determining neighbours (see Triangulation class definition)
				if (i == 0 && j ==0) {// bottom left corner of XY
					faces.at(count).at(3) = top_lower_right;
					faces.at(count).at(5) = count;// Triang. class convention (points to itself)
				}
				else if (i == (M-1) && j == 0){// upper left corner of XY
					faces.at(count).at(3) = count;
					faces.at(count).at(5) = count;
				}
				else if (i == (M-1)){// upper side but not left corner of XY
					faces.at(count).at(3) = count;
					faces.at(count).at(5) = left_lower_right;
				}
				else{// has subrectangles to the left and up
					faces.at(count).at(3) = top_lower_right;
					faces.at(count).at(5) = left_lower_right;
				}
								
				count++;
			}
		}
		
		return(Triangulation<std::vector<double>>(nodes,faces));
		
	}// PlanarTriangulation()


Triangulation<std::vector<double>> PlanarTriangulationStitchTopBottom(
		std::vector<double> X,
		std::vector<double> Y,
		int N)
	{
		/*
			Description: Same as PlanarTriangulation() but also connects upper part of XY (i.e. vertices of the form
				(x,ymax)) to the lower part of XY (i.e. vertices of the form (x,ymin)). To do this, another row of subrectangles
				is appended at the top, with its top two nodes identified with the respective bottom two of the entire XY
				discretization
		*/

		/* Subdividing X,Y */
		
		int M = std::ceil(std::sqrt(N/2)); //number of subintervals in X,Y
		
		std::vector<std::vector<double>> nodes((M+1)*(M+1), std::vector<double>(2,0));// nodes for Triangulation class
		// nodes.at(k*(M+1)) through nodes.at(k*(M+1) + M) for k = 0,...,M are the points in the kth row of the planar discretization
		
		int count = 0;// counts nodes
		
		for (int i = 0; i < (M+1); i++){
			for (int j = 0; j < (M+1); j++){
				nodes.at(count).at(0) = X.at(0) + (X.at(1)-X.at(0))*j/M;
				nodes.at(count).at(1) = Y.at(0) + (Y.at(1)-Y.at(0))*i/M;
				count++;
			}
		}
		
		/* Constructing Triangles */
		
		std::vector<std::vector<int>> faces(2*M*M + 2*M,std::vector<int>(6,0));// faces for Triangulation class
		// faces.at(k*2*M) through faces.at(k*2*M + 2*M - 1) for k = 0,...,M are the triangles in the kth row of subrectangles (kth row
		// 		is below (k+1)th row)
		// faces.at(k*2*M + 2*j) for j = 0,...,M is the lower right triangle of the jth subrectangle in the kth row of subrectangles
		// faces.at(k*2*M + 2*j + 1) for j = 0,...,M is the upper left triangle of the jth subrectangle in the kth row of subrectangles
				
		//Indexes of triangles in faces, relative to current subrectangle
		int down_upper_left = 0;// index (count) of the triangle which is located at the upper left of the 
								// subrectangle below (down) the current subrectangle
		int top_lower_right = 0;// etc.
		int right_upper_left = 0;
		int left_lower_right = 0;
		
		count = 0;// counts faces
		
		for (int i = 0; i < (M+1); i++){
			for (int j = 0; j < M; j++){
				// (i,j) denotes the subrectangle of X \times Y that is on the ith row (ith row is below (i+1)th row)
				//		and jth column (jth columns is to the left of (j+1)th column)
				
				down_upper_left = (i-1)*2*M + 2*j + 1;
				top_lower_right = (i+1)*2*M + 2*j;
				right_upper_left = i*2*M + 2*(j+1) + 1;
				left_lower_right = i*2*M + 2*(j-1);
				
				/* lower right triangle */
				
				faces.at(count).at(2) = i*(M+1) + j;// lower left vertex
				faces.at(count).at(1) = i*(M+1) + j + 1;// lower right vertex
				faces.at(count).at(0) = (i+1)*(M+1) + j + 1;// upper right vertex
				
				faces.at(count).at(4) = count + 1;// opposite triangle to low right corner is the other half of the subrectangle

				//Determining neighbors (see Triangulation class definition)
				if (i > 0 && i < M && j < (M-1)){// triangle has neighbors below and to the right and is not part of 
												 // appended top row to be stitched
					faces.at(count).at(5) = right_upper_left;
					faces.at(count).at(3) = down_upper_left;
				}
				else if (i == 0){//triangle is at bottom side of XY
					faces.at(count).at(3) = M*2*M + 2*j +1;// connect lower right triangle of bottom side of XY to
														   // upper left triangle of top side of XY
					if (j == (M-1)){//triangle is at right corner of bottom side of XY
						faces.at(count).at(5) = right_upper_left;
					}
					else{
						faces.at(count).at(5) = count;
					}
				}
				else if (i == M){//triangle is at appended top row of XY to be stitched with bottom
					faces.at(count).at(0) = (0)*(M+1) + j + 1;// identify upper right vertex with the respective one
														   		// at the bottom of XY
					
					faces.at(count).at(3) = down_upper_left;
					if (j == (M-1)){//triangle is at right corner of bottom side of XY
						faces.at(count).at(5) = right_upper_left;
					}
					else{
						faces.at(count).at(5) = count;
					}
				}
				else if (j == (M-1)){//triangle is at right side of XY but not at bottom neither at appended top row
					faces.at(count).at(3) = down_upper_left;
					faces.at(count).at(5) = count;
				}
				
				count++;
				
				/* upper left triangle */
				
				faces.at(count).at(0) = i*(M+1) + j;// lower left vertex
				faces.at(count).at(1) = (i+1)*(M+1) + j;// upper left vertex
				faces.at(count).at(2) = (i+1)*(M+1) + j + 1;// upper right vertex
				
				faces.at(count).at(4) = count - 1;// opposite triangle to upper left corner is the other half of the subrectangle

				// Determining neighbor (see Triangulation class definition)
				if (i < M && j > 0){// triangle has neighbors to the left and up
					faces.at(count).at(3) = top_lower_right;
					faces.at(count).at(5) = left_lower_right;
				}
				else if (i == M){// triangle is at top row of XY (the one that was appended to be stitched with bottom)
					faces.at(count).at(2) = (0)*(M+1) + j + 1;// identify upper right vertex with the respective one
														   	    // at the bottom of XY
					faces.at(count).at(1) = (0)*(M+1) + j;// identify upper left vertex with the respective one
														   	    // at the bottom of XY
					
					faces.at(count).at(3) = 0*2*M + 2*j;// connect upper left triangle of appended top side of XY to
														// lower right triangle on bottom side of XY
					if (j == 0){//left side of XY
						faces.at(count).at(5) = count;
					}
					else{
						faces.at(count).at(5) = left_lower_right;
					}
				}
				else if (j == 0){//triangle is at left side of XY but not at top appended row
					faces.at(count).at(3) = top_lower_right;
					faces.at(count).at(5) = count;// Triang. class convention (points to itself)
				}
								
				count++;
			}
		}
		
		return(Triangulation<std::vector<double>>(nodes,faces));

	}// PlanarTriangulationStitchTopBottom()

	Triangulation<std::vector<double>> PlanarTriangulationStitchLeftRight(
		std::vector<double> X,
		std::vector<double> Y,
		int N)
	{
		/*
			Description: Same as PlanarTriangulationStitchTopBottom() but instead connects left side of XY (i.e. vertices of the form
				(xmin,y)) to the right side of XY (i.e. vertices of the form (xmax,y)). To do this, the PlanarTriangulationStitchLeftRight()
				is called on Y X instead of X Y, and then the nodes are adjusted back to the original coordinate system
		*/

		Triangulation<std::vector<double>> T = PlanarTriangulationStitchTopBottom(Y,X,N);

		//flip nodes
		double dummy = 0.0;
		for(int i = 0; i < T.nodes.size(); i++){
			dummy = T.nodes.at(i).at(1);
			T.nodes.at(i).at(1) = T.nodes.at(i).at(0);
			T.nodes.at(i).at(0) = dummy;
		}

		return T;
	}

}// namespace



#endif// PLANARTRIANGULATION_HPP
