#ifndef WRITETOFILE_HPP
#define WRITETOFILE_HPP

#include <string>
#include <vector>
#include <fstream>

template<typename scalar>
bool WriteToFile(
	scalar X,
	std::string filename,
	std::string message = ""
);

template<typename scalar>
bool WriteToFile(
	std::vector<scalar> X,
	std::string filename,
	bool write_as_row,
	std::string message = ""
);

template<typename scalar>
bool WriteToFile(
	std::vector<std::vector<scalar>> X,
	std::string filename,
	std::string message = ""
	);

template<typename scalar>
bool WriteToFile(
	std::vector<std::vector<scalar>> X,
	std::string filename,
	std::string message
	)
{
	/*
		Description: Writes each element of X to a row in filename. Each element of X.at(i) is separated by
			a space in the respective line in filename
			
		Input:
			
			- std::vector<std::vector<double>> X
				X is written in filename. X.at(i) is written at ith row of filename. Each element of X.at(i) is
					separated by a space
			- std::string filename
				file to write X at. Must include extension (ie "data.txt" or "data.dat")
	
		Output:
		
			- bool out
				Signifies completion
				
		Notes: If filename exists, then it is overwritten
	*/
	
	std::fstream file;
	bool out;
	
	file.open(filename,std::ios::out);
	if (file.is_open()){
	
		out = true;
	
		if (message != "") {//print message
			file << "#" << message << std::endl;
		}

		for (int i = 0;i < X.size();i++){
			for (int j = 0; j < X.at(i).size(); j++){
				file << X.at(i).at(j);
				if (j != (X.at(i).size() - 1)) file << " ";// Do not add space after last element of each row in file
			}
			if (i != (X.size() - 1)) file << std::endl;// Do not add newline after last row in file
		}
	}
	else{
		out = false;	
	}
	file.close();	
	
	return(out);
	
}// WriteToFile()

template<typename scalar>
bool WriteToFile(
	scalar X,
	std::string filename,
	std::string message
) {
	return WriteToFile(std::vector<std::vector<scalar>> {{X}}, filename, message);
}

bool WriteToFile(std::string message, std::string filename) {
	
	/*
		Description: create filename and write message in first line

		If filename already exists, it is deleted
	*/
	std::fstream file;
	bool out;

	file.open(filename, std::ios::out);
	if (file.is_open()) {
		out = true;
		file << message;
	}
	else {
		out = false;
	}
	file.close();

	return out;
}

template<typename scalar>
bool WriteToFile(
	std::vector<scalar> X,
	std::string filename,
	bool write_as_row,
	std::string message
) {
	/*
		if print_as_row = true, writes X on the first line of filename. Otherwise,
		each element of X occupies one line of filename
	*/

	std::vector<std::vector<scalar>> Y;
	if(write_as_row){//write X as a row
		Y.push_back(X);
	}
	else{//write X as a column
		Y = std::vector<std::vector<scalar>>(X.size(),std::vector<scalar>(1,scalar(0)));
		for(int i = 0; i < X.size(); i++){
			Y.at(i) = {X.at(i)};
		}
	}

	return WriteToFile(Y, filename, message);
}

#endif// WRITETOFILE_HPP
