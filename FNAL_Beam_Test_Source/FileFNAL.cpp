/*
 * FileFNAL.cpp
 *
 *  Created on: Apr 1, 2019
 *      Author: dylan
 */

#include "../FNAL_Beam_Test_Source/FileFNAL.h"

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;


vector<vector<int>> read_2d_file(string file_path) {
	ifstream file(file_path);
	if(file.fail()) {
		cout << "Could not find file." << endl;
	}
	vector<vector<int>> data;
	int temp_int;
	string temp_string;
	while(getline(file, temp_string)) {
		vector<int> temp_row;
		istringstream temp_stream(temp_string);
		while(temp_stream >> temp_int) {
			temp_row.push_back(temp_int);
		}
		data.push_back(temp_row);
	}

	return(data);
}

