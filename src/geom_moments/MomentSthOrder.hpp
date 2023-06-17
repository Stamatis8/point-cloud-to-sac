#ifndef MOMENTSTHORDER_HPP
#define MOMENTSTHORDER_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <vector>

#include "NchooseK_cache.hpp"
#include "J_cache.hpp"
#include "stl_util/stl2vec.hpp"

double MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int i,
	int j,
	int k,
	int degree,
	bool is_translation_invariant,
	bool is_scaling_invariant
	);

double MomentSthOrder(
	std::string filename,
	int i,
	int j,
	int k,
	int degree,
	bool is_translation_invariant,
	bool is_scaling_invariant
	);

double MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int i,
	int j,
	int k,
	int degree,
	bool is_translation_invariant = false,
	bool is_scaling_invariant = false
	)
{
	/*
		Description: The input variable triangles describes the triangulation of a surface S, which encloses a volume G. This function calculates the
			s = (i+j+k) order geometric moment of G, according to Pozo's et.al. - *Approximate Series Algorithm For Geometric Moments* which can be found in
			section 4 of [1]. Further, the translation and scaling invariant moments of G can also be calculated as in section 3.2.1 of [2] 
					
		Note: exact result for degree == i + j + k
		
		References:
		
			[1] Efficient 3D Geometric and Zernike moments computation from unstructured surface meshes J. M. Pozo, M. C. Villa-Uriol,
				A. F. Frangi, Senior Member, IEEE
				
			[2] Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation, S. Khan, 
				P. Kaklis, 2022
		
		Input:
		
			- std::vector<std::vector<std::vector<double>>> triangles
				triangles.at(p) is the pth triangle with vertices triangles.at(p).at(q), q=0,1,2. Each vertex must have three elements
			- int i
				i >= 0
			- int j
				j >= 0
			- int k
				k >= 0
			- int degree
				0 <= degree <= i+j+k
				exact result if degree == i+j+k
			- bool is_translation_invariant
				if true, calculates the translation invariant moment form of G (ie centered at it's center of volume)
				By definition all first order moments are equal to zero, so if i+j+k = 1, set M = 0
				default to false
			- bool is_scaling_invariant
				if true, calculates the scaling invariant moment form of G 
				By definition volume is equal to 1, so if i+j+k = 0, set M = 1
				default to false
		Output:
		
			- double M
				the s = (i+j+k) order moment of G, according to options is_translation_invariant and is_scaling_invariant
	*/

	/* Resolving some cases by definition */

	if ((i+j+k)==0 && is_scaling_invariant){
		// if scaling invariant, then by definition volume is equal to 1
	
		return 1;
	}
	else if ((i+j+k)==1 && is_translation_invariant){
		// if centered at volume centroid, then by definition 1st order moments are equal to zero
		
		//return 0;
	}
	
	/* Initializing variables */
	
	int T = triangles.size();// number of triangles
	
	std::vector<std::vector<double>> p_bar = // p_bar.at(i) = centroid of ith triangle
		std::vector<std::vector<double>>(T,std::vector<double>(3,0)); 
	
	std::vector<std::vector<double>> u = // Sec.4 of [1]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	std::vector<std::vector<double>> v = // Sec.4 of [1]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	std::vector<std::vector<double>> u_hat = // Sec.4 of [1]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	std::vector<std::vector<double>> v_hat = // Sec.4 of [1]
		std::vector<std::vector<double>>(T,std::vector<double>(3,0));
	
	auto mag = [](std::vector<double> a){// magnitude lamda function
		return (std::sqrt(a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2)));
	};
	
	auto AreaC = [](std::vector<double> u, std::vector<double> v){// triangle area lamda function
		double n1, n2, n3;
		n1 = (u.at(1)*v.at(2) - u.at(2)*v.at(1)); 
		n2 = (u.at(2)*v.at(0) - u.at(0)*v.at(2));
		n3 = (u.at(0)*v.at(1) - u.at(1)*v.at(0));
		return (std::sqrt(n1*n1 + n2*n2 + n3*n3)/2);
	};
	
	auto VolC = [](std::vector<double> p1, std::vector<double> p2, std::vector<double> p3){// tetrahedra oriented volumes
		return ((
				  p1.at(0)*(p2.at(1)*p3.at(2)-p3.at(1)*p2.at(2))
				- p1.at(1)*(p3.at(2)*p2.at(0)-p2.at(2)*p3.at(0))
				+ p1.at(2)*(p2.at(0)*p3.at(1)-p3.at(0)*p2.at(1))
				)/6);
	};
	
	/* If is_translation_invariant, translate mesh */
	
	double V = 0; //volume
	
	if (is_translation_invariant){
	
		V = MomentSthOrder(triangles,0,0,0,0,false,false);
		
		double Cx = MomentSthOrder(triangles,1,0,0,1,false,false)/V;// centroid x
		double Cy = MomentSthOrder(triangles,0,1,0,1,false,false)/V;// centroid y
		double Cz = MomentSthOrder(triangles,0,0,1,1,false,false)/V;// centroid z
		
		for (int i = 0; i < T; i++){
			//First Vertex
			triangles.at(i).at(0).at(0) -= Cx;
			triangles.at(i).at(0).at(1) -= Cy;
			triangles.at(i).at(0).at(2) -= Cz;
			
			//Second Vertex
			triangles.at(i).at(1).at(0) -= Cx;
			triangles.at(i).at(1).at(1) -= Cy;
			triangles.at(i).at(1).at(2) -= Cz;
			
			//Third Vertex
			triangles.at(i).at(2).at(0) -= Cx;
			triangles.at(i).at(2).at(1) -= Cy;
			triangles.at(i).at(2).at(2) -= Cz;
		}
	}
	
	/* Calculating p_bar, u, v */
	
	for (int i = 0; i < T; i++){
		p_bar.at(i).at(0) = (triangles.at(i).at(0).at(0) + triangles.at(i).at(1).at(0) + triangles.at(i).at(2).at(0))/3;
		p_bar.at(i).at(1) = (triangles.at(i).at(0).at(1) + triangles.at(i).at(1).at(1) + triangles.at(i).at(2).at(1))/3;
		p_bar.at(i).at(2) = (triangles.at(i).at(0).at(2) + triangles.at(i).at(1).at(2) + triangles.at(i).at(2).at(2))/3;
		
		u.at(i).at(0) = triangles.at(i).at(0).at(0) - triangles.at(i).at(2).at(0);
		u.at(i).at(1) = triangles.at(i).at(0).at(1) - triangles.at(i).at(2).at(1);
		u.at(i).at(2) = triangles.at(i).at(0).at(2) - triangles.at(i).at(2).at(2);
		
		v.at(i).at(0) = triangles.at(i).at(1).at(0) - triangles.at(i).at(2).at(0);
		v.at(i).at(1) = triangles.at(i).at(1).at(1) - triangles.at(i).at(2).at(1);
		v.at(i).at(2) = triangles.at(i).at(1).at(2) - triangles.at(i).at(2).at(2);
	}
	
	/* Calculating lamda */
	
	double lamda = 0;// Sec.4 of [2]
	
	double facet_area_sum = 0;// see nominator of lamda definition
	double weight_facet_area_sum = 0;// see denominator of lamda definition
	double A = 0;// area of current rectangle
	double p_bar_mag = 0;// magnitude of current p_bar
	
	for (int i = 0; i < T; i++){
		A = AreaC(u.at(i),v.at(i));
		facet_area_sum += A;
		weight_facet_area_sum += mag(p_bar.at(i))*A;
	}
	
	lamda = std::sqrt(4*std::pow(facet_area_sum,3)/(T*std::sqrt(3)))/weight_facet_area_sum;
	
	/* Calculating u_hat, v_hat */
	
	for (int i = 0; i < T; i++){
		u_hat.at(i).at(0) = u.at(i).at(0)/lamda;
		u_hat.at(i).at(1) = u.at(i).at(1)/lamda;
		u_hat.at(i).at(2) = u.at(i).at(2)/lamda;
		
		v_hat.at(i).at(0) = v.at(i).at(0)/lamda;
		v_hat.at(i).at(1) = v.at(i).at(1)/lamda;
		v_hat.at(i).at(2) = v.at(i).at(2)/lamda;
	}
	
	/* Initialize M calculation */
	
	double M = 0;// final result
	int n = i + j + k;
	double S_ijk;// Sec.4 of [1]
	double S_ijk_r;// Sec.4 of [1]
	
	int i1;// Sec.4 of [1]
	int j1;// etc.
	int k1;
	int k3;
	int kplus;
	int i3;
	int j3;
	
	NchooseK_cache<double> NK;// N choose K cache
	NK.get(std::max(i,std::max(j,k)),0);// Precomputing N choose m pairs for all m and N <= max(i,j,k)
	
	J_cache J;// initializing Jab integral cache. See appendix 1 of [1]
	
	for (int c = 0; c < T; c++){// facets
		
		/* S_ijk calculation */
		
		S_ijk = 0;
		for (int r = 0; r <= degree; r++){
		
			/* S_ijk^r calculation */
			
			S_ijk_r = 0;
			for (int iplus = std::max(0,r-j-k); iplus <= std::min(i,r); iplus++){
				for (int jplus = std::max(0,r-k-iplus);jplus <= std::min(j,r-iplus); jplus++){
					for (int i2 = 0; i2 <= iplus; i2++){
						for (int j2 = 0; j2 <= jplus; j2++){
							kplus = r - iplus - jplus;
							for (int k2 = 0; k2 <= kplus; k2++){
							
								i1 = i - iplus;
								j1 = j - jplus;
								k1 = k - kplus;
								i3 = iplus - i2;
								j3 = jplus - j2;
								k3 = kplus - k2;
							
								S_ijk_r += NK.get(i,iplus)*NK.get(j,jplus)*NK.get(k,kplus)
										 * NK.get(iplus,i2)*NK.get(jplus,j2)*NK.get(kplus,k2)
										 * J.get(i2+j2+k2,r-i2-j2-k2)
										 * std::pow(p_bar.at(c).at(0),i1)*std::pow(p_bar.at(c).at(1),j1)*std::pow(p_bar.at(c).at(2),k1)
										 * std::pow(u_hat.at(c).at(0),i2)*std::pow(u_hat.at(c).at(1),j2)*std::pow(u_hat.at(c).at(2),k2)
										 * std::pow(v_hat.at(c).at(0),i3)*std::pow(v_hat.at(c).at(1),j3)*std::pow(v_hat.at(c).at(2),k3);
									
							}//k2
						}//j2
					}//i2
				}// jplus
			}// iplus
			
			S_ijk += std::pow(lamda,r)*S_ijk_r;
		}
		
		M += 6*VolC(triangles.at(c).at(0),triangles.at(c).at(1),triangles.at(c).at(2))*S_ijk;
	}
	
	M /= (n+3);
	
	/* If is_scaling_invariant, divide final moment M by volume^(1+n/3) */
	
	if (is_scaling_invariant){
		if(!is_translation_invariant){
			// volume has not been calculated
			
			V = MomentSthOrder(triangles,0,0,0,false,false);
		}
		
		M /= std::pow(V,(1+(i+j+k)/3));
	}
	
	return M;
}

double MomentSthOrder(
	std::string filename,
	int i,
	int j,
	int k,
	int degree,
	bool is_translation_invariant = false,
	bool is_scaling_invariant = false
	){
	/*
		Description: Same as MomentSthOrder, but triangles are imported from stl file
	*/
	
	return MomentSthOrder(stl2vec(filename),i,j,k,degree,is_translation_invariant,is_scaling_invariant);
}

#endif // MOMENTSTHORDER_HPP

/*

double DEBUG_TEMP(
		double A, 
		double B, 
		double C, 
		double D,
		double E,
		double F,
		double G,
		double H,
		double I,
		int a,
		int b,
		int c){
		double result = 0;
		int id = 0;
		
		NchooseK_cache NK;// N choose K cache
		NK.get(std::max(a,std::max(b,c)),0);// Precomputing N choose k pairs for all k and N <= max(p,q,r)
		
		for (int ia = 0; ia <= a; ia++){
			for (int ib = 0; ib <= b; ib++){
				for (int ic = 0; ic <= c; ic++){
				
					id = a + b + c + 1 - ia - ib - ic;
				
					for (int ja = 0; ja <= ia; ja++){
						for (int jb = 0; jb <= ib; jb++){
							for (int jc = 0; jc <= ic; jc++){
								for (int jd = 0; jd <= id; jd++){
									result += NK.get(a,ia)*NK.get(b,ib)*NK.get(c,ic)*NK.get(ia,ja)*NK.get(ib,jb)*NK.get(ic,jc)*NK.get(id,jd)
											* std::pow(A,a-ia)*std::pow(D,b-ib)*std::pow(G,c-ic)
											* std::pow(B,ia-ja)*std::pow(C,ja)*std::pow(E,ib-jb)
											* std::pow(F,jb)*std::pow(H,ic-jc)*std::pow(I,jc)
											* std::pow(-1,jd)/(id*(a+b+c+2-ja-jb-jc-jd));
								}//jd
							}//jc
						}//jb
					}//ja
					
				}//ic
			}//ib
		}// ia
		
		return result;
}// DEBUG_TEMP()


double old_MomentSthOrder(
	std::vector<std::vector<std::vector<double>>> triangles,
	int p,
	int q,
	int r,
	bool is_invariant
	)
{
	/*
		Description: The input variable triangles describes the triangulation of some surface S, which encloses some volume G. This function calculates the
			s = (p+q+r) order moment of G, as described in the appendix of 
			
			[1] Geometric Moment-Dependent Global Sensitivity Analysis without Simulation Data: Application to Ship Hull Form Optimisation - S. Khan - P. Kaklis -
			2022
		
		Input:
		
			- std::vector<std::vector<std::vector<double>>> triangles
				triangles.at(i) is the ith triangle with vertices triangles.at(i).at(j), j=0,1,2. Each vertex must have three elements
			- int p
				p >= 0
			- int q
				q >= 0
			- int r
				r >= 0
			- bool is_invariant
				calculates the invariant moment form of G (invariant to scaling and translation transformations) as described in eq.8 of [1] 
	
		Output:
		
			- double M
				the s = (p+q+r) moment of G
				
		Note: The method used in this implementation is the analytical evaluation of the integral of f (see appendix of [1]) over each triangle.
	*\/
	
	double M = 0;
	double I_t = 0;// The integral of f (see appendix of [1]) over the triangle t
	
	double A,B,C,D,E,F,G,H,I;// Factors in analytitcal expression of I_t
	double G1,G1_bar,G1_tilde;// Summation terms in I_t expression
	double power_row;//common factor in G1,G1_bar,G1_tilde	
	
	std::vector<double> a(3,0), b(3,0), c(3,0);// triangle is parametrized as (u,v) = au+bv+c
	double a_mag,b_mag;// Factors in analytical expression of I_t
	std::vector<double> n(3,0);// a,b exterior product
	int id; // Term in analytical expression of I_t
	
	NchooseK_cache NK;// N choose K cache
	NK.get(std::max(p,std::max(q,r)),0);// Precomputing N choose k pairs for all k and N <= max(p,q,r)
	
	double volume = 0;// In case invariant M is requested
	std::vector<double> centroid(3,0);// In case invariant M is requested
	
	for (int t = 0; t < triangles.size(); t++){
		
		/* Evaluating A through I, a, b, c, n *\/
		
		a.at(0) = triangles.at(t).at(1).at(0) - triangles.at(t).at(0).at(0);
		a.at(1) = triangles.at(t).at(1).at(1) - triangles.at(t).at(0).at(1);
		a.at(2) = triangles.at(t).at(1).at(2) - triangles.at(t).at(0).at(2);
		
		b.at(0) = triangles.at(t).at(2).at(0) - triangles.at(t).at(0).at(0);
		b.at(1) = triangles.at(t).at(2).at(1) - triangles.at(t).at(0).at(1);
		b.at(2) = triangles.at(t).at(2).at(2) - triangles.at(t).at(0).at(2);
		
		c = triangles.at(t).at(0); 
		
		a_mag = std::sqrt(a.at(0)*a.at(0) + a.at(1)*a.at(1) + a.at(2)*a.at(2));
		b_mag = std::sqrt(b.at(0)*b.at(0) + b.at(1)*b.at(1) + b.at(2)*b.at(2));
		
		//n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1))/(a_mag*b_mag); 
		//n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2))/(a_mag*b_mag);
		//n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0))/(a_mag*b_mag);
		
		n.at(0) = (a.at(1)*b.at(2) - a.at(2)*b.at(1)); 
		n.at(1) = (a.at(2)*b.at(0) - a.at(0)*b.at(2));
		n.at(2) = (a.at(0)*b.at(1) - a.at(1)*b.at(0));
		double n_mag = std::sqrt(n.at(0)*n.at(0) + n.at(1)*n.at(1) + n.at(2)*n.at(2));
		n.at(0) /= n_mag;
		n.at(1) /= n_mag;
		n.at(2) /= n_mag;

		
		A = a.at(0);
		B = b.at(0);
		C = c.at(0);
		D = a.at(1);
		E = b.at(1);
		F = c.at(1);
		G = a.at(2);
		H = b.at(2);
		I = c.at(2);
	
		if (is_invariant){//amending A through I (translational invariance)
			volume = old_MomentSthOrder(triangles,0,0,0,false);
			
			/***** Todo: Handle volume == 0 *****\/
			
			centroid.at(0) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(1) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			centroid.at(2) = old_MomentSthOrder(triangles,1,0,0,false)/volume;
			
			C -= centroid.at(0);
			F -= centroid.at(1);
			I -= centroid.at(2);
		}		
		
		/* Evaluating I_t *\/
		
		/* DEBUG: Evaluating I_t another way *\/
		
		double term1 = 0;
		double term2 = 0;
		double term3 = 0;
		
		term1 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p+1,q,r);
		term2 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p,q+1,r);
		term3 = DEBUG_TEMP(A,B,C,D,E,F,G,H,I,p,q,r+1);
		
		term1 *= n.at(0)/(p+1);
		term2 *= n.at(1)/(q+1);
		term3 *= n.at(2)/(r+1);
		
		//I_t = a_mag*b_mag*(term1 + term2 + term3)/3;
		
		double ab = a.at(0)*b.at(0) + a.at(1)*b.at(1) + a.at(2)*b.at(2);
		I_t = std::sqrt(a_mag*a_mag*b_mag*b_mag - ab*ab)*(term1+term2+term3)/3;
				
		M += I_t;
	}//triangles
	
	if (is_invariant){//amending M (scaling invariance)
		M /= std::pow(volume,1+(p+q+r)/3);
	}
	
	return M;

}// MomentSthOrder()
*/
