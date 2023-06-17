#ifndef J_CACHE_HPP
#define J_CACHE_HPP

#include <vector>
#include <cmath>

struct J_cache{
	/*
		Description:
			This struct is intended as a lookup table for the Jab integrals, as calculated in the appendix of [1]. The Jab integral is computed only once,
			at the first time it is requested. 
			
			To access Jab, use the method .get(a,b)
			
			Jab if a function of Kab, where Kab can be computed recursively
			
			The J and K tables will have dubs with bool values to indicate which slots in them have actually been calculated.
		
			The interface only gives access to the .get method (ie only the J table) 
			
			Note: Jab = Jba and Kab = Kba, so J, K tables only have lower triangular values
		
		References:
		
			[1] Efficient 3D Geometric and Zernike moments computation from unstructured surface meshes J. M. Pozo, M. C. Villa-Uriol,
				A. F. Frangi, Senior Member, IEEE
	*/
public:
	double get(int a, int b){
		/*
			Description:
				returns Jab. If necessary calculates it and saves it at this->J
		*/
		if (b > a) {// Jab = Jba. Bypassing Kab for b > a
			return this->get(b,a);
		}
		else{
			if ((a + 1) <= this->J.size() && (b + 1) <= this->J.at(a).size() && this->J_dub.at(a).at(b)){
				// Jab has been calculated
				return this->J.at(a).at(b);
			}
			else{
			
				//calculate Jab
				
				double Jab = (this->get_K(a,b+1)-std::pow(-1,b+1)*this->get_K(a,0))/(std::pow(3,a+b+2)*(b+1));
				
				this->add_J(a, b, Jab);//save Jab
				
				return Jab;
			}
		}
	}// get()
		
private:
	std::vector<std::vector<double>> J;
	// J.at(a).at(b) equals Jab
	
	std::vector<std::vector<bool>> J_dub;
	// size of J_dub is at all times equal to J
	// if J_dub.at(a).at(b) = true, Jab has been calculated, else not
	
	std::vector<std::vector<double>> K;
	// K.at(a).at(b) equals Kab
	
	std::vector<std::vector<bool>> K_dub;
	// size of K_dub is at all times equal to K
	// if K_dub.at(a).at(b) = true, Kab has been calculated, else not

	double get_K(int a, int b){
		/*
			Description:
				returns Kab. If necessary calculates it recursively and saves it at this->K. All
				recursive instances Kij are saved to K
		*/
		if (b > a) {// Kab = Kba. Bypassing Kab for b > a
			return this->get_K(b,a);
		}
		else{
			if ((a + 1) <= this->K.size() && (b + 1) <= this->K.at(a).size() && this->K_dub.at(a).at(b)){
				// Kab has been calculated
				return this->K.at(a).at(b);
			}
			else{
			
				//calculate Kab
				
				double Kab = 0;
				if (b == 0){
					//base case
					Kab = (std::pow(2,a+1) - std::pow(-1,a+1))/(a+1);
				}
				else{
					//recursion
					Kab = this->get_K(a,b-1) - this->get_K(a+1,b-1);
				} 
				
				this->add_K(a, b, Kab);//save Kab
				
				return Kab;
			}
		}
	}// get_K()
	
	void add_J(int a, int b, double Jab){
		/*
			Description:
				adds Jab to J. If necessary size of Jab is increased. This method is useful because
				compatibility of J_dub with J must be satisfied
		*/
		
		if((a + 1) > this->J.size()){
			this->J.resize(a+1,std::vector<double>(1,0));
			this->J_dub.resize(a+1,std::vector<bool>(1,false));
		}
		
		if((b + 1) > this->J.at(a).size()){
			this->J.at(a).resize(b+1,0);
			this->J_dub.at(a).resize(b+1,false);
		}
		
		this->J.at(a).at(b) = Jab;
		this->J_dub.at(a).at(b) = true;
	}// add_J()
	
	void add_K(int a, int b, double Kab){
		/*
			Description:
				adds Kab to K. If necessary size of Kab is increased. This method is useful because
				compatibility of K_dub with K must be satisfied
		*/
		
		if((a + 1) > this->K.size()){
			this->K.resize(a+1,std::vector<double>(1,0));
			this->K_dub.resize(a+1,std::vector<bool>(1,false));
		}
		
		if((b + 1) > this->K.at(a).size()){
			this->K.at(a).resize(b+1,0);
			this->K_dub.at(a).resize(b+1,false);
		}
		
		this->K.at(a).at(b) = Kab;
		this->K_dub.at(a).at(b) = true;
	}// add_K()
	
}; // J_cache

#endif // J_CACHE_HPP
