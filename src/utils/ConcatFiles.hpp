#ifndef CONCATFILES_HPP
#define CONCATFILES_HPP

#include <string>
#include <vector>

#include "ReadFile.hpp"
#include "WriteToFile.hpp"

void ConcatFilesNumeric(std::string File1, std::string File2, std::string FileOut, bool vertical = true){
	/*
		Description: Copy Contents of File1 and File2 to FileOut in that order. Contents must be numeric
			(i.e. File3 starts with File1 and continues with File2)

			if vertical = true then FileOut starts with all lines of File1 and then below them continues with
				all lines of File 2

			if vertical = false then each line of FileOut starts with the respective line of File1 and continues
				with the respective line of File 2. If the sizes of files 1 and 2 differ, then the extra 
				"necessary" lines are treated as empty strings.
	*/
	
	//Read Files
	
	std::vector<std::vector<double>> F = ReadFile(File1);
	std::vector<std::vector<double>> F_ = ReadFile(File2);
	
	if(vertical){

		F.insert(F.end(),F_.begin(),F_.end());
		
		WriteToFile(F,FileOut);

	}
	else{
		
		std::vector<std::vector<double>> F__;

		for(int i = 0; i<std::min(F.size(),F_.size()); i++){
			//while there is data in both files, add both of them line by line normally
			F__.push_back(F.at(i));
			F__.at(i).insert(F__.at(i).end(),F_.at(i).begin(),F_.at(i).end());
		}

		if(F.size() > F_.size()){//if F > F_, insert the remaining part of F at F__
			for(int i = F_.size(); i < F.size(); i++){
				F__.push_back(F.at(i));
			}
		}
		else{//if F < F_, insert the remaining part of F_ at F__
			for(int i = F.size(); i < F_.size(); i++){
				F__.push_back(F_.at(i));
			}
		}

		WriteToFile(F__, FileOut);
	}
		
}

#endif //CONCATFILES_HPP