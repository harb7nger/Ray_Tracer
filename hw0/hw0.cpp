#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
using namespace std;

// struct to store RGB values
struct RGB {
  	int red;
  	int green;
  	int blue;
}; 

// parses the file and gets the image size
pair<int, int> getDimensions(string fileName) {
    ifstream infile(fileName);
    int a, b;
	string imsize;
	infile >> imsize >> a >> b;

	//checks if the input is valid
	if ((imsize != "imsize") || a <= 0 || b <= 0) return make_pair(-1, -1);
	return make_pair(a, b);
}

// creates a new file which ends with .ppm from the input file name and puts in the header information
string createFileAndWriteHeader(string fileName, pair<int, int>& dimensions) {

	// parse the file name to create output file
	// https://stackoverflow.com/questions/757933/how-do-you-change-the-filename-extension-stored-in-a-string-in-c
	string outputFile = fileName.substr(0, fileName.find_last_of('.'))+".ppm";
	
	ofstream file;
	file.open(outputFile);
	file << "P3" << endl;
	file << "#simple color patterns" << endl;
	file << dimensions.first << " " << dimensions.second << endl; 
	file << "255" << endl;
	file.close();

	return outputFile;
}

// helper to create patterns, uses weighted values for RGB, 1 for red, 2 for green & 3 for blue 
vector<vector<RGB>> getPatterns(pair<int, int>& dimensions) {
	vector<vector<RGB>> colors;
	for (int i = 0; i < dimensions.second; i++) {	
	    vector<RGB> row;
            for (int j = 0; j < dimensions.first; j++) {
		     	struct RGB val;
		     	val.red = (i*1+ j%100)%255; val.green = (i*2+j%100)%255; val.blue = (i*3+j%100)%255;
		     	row.push_back(val);
            }
	    colors.push_back(row);
	}
	return colors;
}

// calls on the above helper function to get the pattern in a vector form and prints it to the filename specified
void drawPatterns(string fileName, pair<int, int>& dimensions) {

	ofstream file;
	file.open(fileName, ios_base::app);
	
	// gets the patterns in vector form
	vector<vector<RGB>> pattern = getPatterns(dimensions);

	// write pattern to file
	for (int i = 0; i < dimensions.second;i++) {
	    for (int j = 0; j < dimensions.first; j++) {
			file << pattern[i][j].red << " " << pattern[i][j].green << " " << pattern[i][j].blue << endl; 
	    }
	}
	file.close();
	return;
}    


int main (int argc, char** argv) {
	// check for arguments
	if (argc < 2) {
		cout << "Too few arguments please retry =D" << endl;
		return 0;
    }

    // if there are too many arguments
    if (argc > 2) {
    	cout << "Too many arguments please retry" << endl;
    	return 0;
    }
    
	// get dimensions
	pair<int, int> dimensions = getDimensions(argv[1]);
	
	if (dimensions.first < 0 ) {
		cout << "Please check the input" << endl;
		return 0;
	}

	// get output filename
	string outputFile = createFileAndWriteHeader(argv[1], dimensions);
	// draw patterns in output file
	drawPatterns(outputFile, dimensions);
	return 0;
} 
