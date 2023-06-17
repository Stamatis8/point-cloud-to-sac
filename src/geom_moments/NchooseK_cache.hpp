#ifndef NCHOOSEK_CACHE_HPP
#define NCHOOSEK_CACHE_HPP

#include <vector>

template<typename scalar = double>
struct NchooseK_cache{
	std::vector<std::vector<scalar>> nk;
	/*
		- the nk table contains (n,k) pairs. nk.at(i).at(j) is equal to (i choose j)
		- throughout the program execution **in the scope that NchooseK is initialized at** , more and more values are calculated and included in nk
		- let M = nk.size. All elements nk.at(i).at(j), for i = 0,...,M-1 are calculated
	*/

	scalar get(int n, int k){
		/*
			Description: checks if (n,k) is calculated in nk. If it is calculated in nk, then nk.at(n).at(k) is returned.
				If it is not calculated then nk.size() - 1 <  n. Then, all values nk.at(i) for i = nk.size() to n are
				calculated
			
			Input:
				
				- int n
					0 <= n
				- int k
					0 <= k <= n
			Output:
			
				- int out
					out = (n choose k)
					
			Note: Since at any point, all values nk.at(i) for i = 1,...,nk.size() are calculated, (n choose k) is calculated using its recursive definition
		*/
		
		/* Check if (n choose k) has been calculated */
	
		if ((nk.size()-1) >= n && nk.size() != 0) {//If nk.size() == 0, then nk.size() - 1 is not negative one in this expression but max(int) - 1. This is 
												   // because nk.size() == 0 is unsigned type which if negative loops around
			return nk.at(n).at(k);
		}
		
		/* If it has not been calculated, calculate it recursively */
		
		int prev = nk.size() - 1;// previously, maximum nk index was nk.at(prev)
	
		nk.resize(n+1);
		
		for (int i = prev+1; i <= n; i++){// nk.at(prev+1) through nk.at(n) must be filled
			
			nk.at(i) = std::vector<scalar>(i+1,0);
			
			nk.at(i).at(0) = 1;// (i choose 0) = 1
			nk.at(i).at(i) = 1;// (i choose i) = 1
			
			if (i != 0){// if not base of recursion
				for (int j = 1; j < i; j++){
					nk.at(i).at(j) = nk.at(i-1).at(j-1) + nk.at(i-1).at(j);
				}
			}
		}
		
		return nk.at(n).at(k);
		
	}// get()
}; // NchooseK_cache



#endif// NCHOOSEK_HPP
