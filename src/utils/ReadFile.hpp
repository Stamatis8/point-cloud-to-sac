#ifndef READFILE_HPP
#define READFILE_HPP

#include <string>
#include <vector>
#include <fstream>

std::vector<std::vector<double>> ReadFile(std::string filename, std::string separator = " ");

std::vector<std::string> separate_string(std::string string, std::string separator);

std::vector<std::vector<double>> ReadFile(std::string filename, std::string separator)
{
	/*
		Description: Saves contents of filename into a vec<vec<item>> structure. Each line is separated into 
            disctinct items via the 'separator'.
			
		Input:
			- std::string filename
				file to write X at. Must include extension (ie "data.txt" or "data.dat")
            - std::string separator
                separator in each line of filename
	
		Output:
		
			- std::vector<std::vector<scalar>> out
				out.at(i) = ith line of 'filename'
				out.at(i).at(j) = jth item in ith line of filename
				
		Notes: If filename exists, then it is overwritten
	*/
	
	std::fstream file;
    
	std::vector<std::vector<double>> out;
    
    std::string line;
    
    std::vector<std::string> words;// each item in 'line' separated by 'separator'
	
    std::vector<double> numbers;// 'words' converted to 'double'
    
	file.open(filename,std::ios::in);
	
    while(std::getline(file,line)){
        words = separate_string(line,separator);
        
        numbers = std::vector<double>(words.size(),0);
        
        for(int i = 0; i < words.size();i++){
            numbers.at(i) = std::stod(words.at(i));
        }
        
        out.push_back(numbers);
        
    }
        
	file.close();	
	
	return(out);
	
}// WriteToFile()

std::vector<std::string> separate_string(std::string string, std::string separator){
    /*
        Description: Separates 'string' according to 'separator' and returns it in a vector. Separator must be a single character.
        
        Example:
            string = "abc def  ghi    j klm nop", separator = " " -> {"abc","def","ghi","j","klm","nop"}
            
    */
    
    std::vector<std::string> out;

    std::string substring = "";
    
    bool was_completing_substring = false;// last iteration appended to substring
    
    for(int i = 0; i < string.length(); i++){
        if(std::string(1,string.at(i)) == separator){
            if(was_completing_substring){// substring has a complete word
                was_completing_substring = false;// no longer completing a substring
                out.push_back(substring);
                substring = "";
            }
            continue;
        }
        else{
            substring = substring + std::string(1,string.at(i));
            was_completing_substring = true;
            
            if(i==(string.length()-1) && was_completing_substring){//finished 'string' and was completing a substring
                out.push_back(substring);
            }
        }
    }
    
    return out;
}

#endif// READFILE_HPP
