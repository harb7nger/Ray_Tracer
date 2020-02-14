#include "hw1.h"
#include <ctype.h>
#include <cmath>
#include <iostream>
#include <float.h>
#include <fstream>
#include <sstream>
#include <utility>
#include <unordered_map>
#include <string>

using namespace std;

// creates a new file which ends with .ppm from the input file name and puts in the header information

bool checkFloat (string word) {
	int decimalCount = 0;
	if (word[0] == '-' || word[0] == '+') word = word.substr(1, word.size()-1);
	for (char& c: word) {
		if (c == '.') {
			decimalCount++;
		} else if (!isdigit(c)) {
			return false;
		}
	}
	return decimalCount	<= 1 ? true:false;
} 

bool checkInt (string word) {
	for (char& c: word) if (!isdigit(c)) {return false;}
	return true;
}

vector<string> tokenizeLine(string line, char delim) {
	stringstream ss(line);
	vector<string> words;
	string word;
	while (getline(ss, word, delim)) {
  		words.push_back(word);
	}
	return words;
}

// returns the sum of two VectorTypes 
VectorType sum(VectorType a, VectorType b) {
	VectorType result;
	result.dx = a.dx + b.dx;
	result.dy = a.dy + b.dy;
	result.dz = a.dz + b.dz;
	return result;
}

// returns the sum of a point and a VectorType
VectorType sum(PointType a, VectorType b) {
	VectorType result;
	result.dx = a.x + b.dx;
	result.dy = a.y + b.dy;
	result.dz = a.z + b.dz;
	return result;
}


// mulitplies a vector with a given scalar
VectorType multiplyScalar(VectorType a, float b) {
	VectorType result = {b*a.dx, b*a.dy, b*a.dz};
	return result;
}

// outputs a given vector
void displayVector(VectorType a) {
	cout << a.dx << "," << a.dy << "," << a.dz << endl;
}

// display ray
void displayRay(RayType a) {
	cout << a.x << "," << a.y << "," << a.z << ":"
		<< a.dx << "," << a.dy << "," << a.dz << endl;
}


// returns a normalized VectorType
VectorType getUnitVector(VectorType vector) {
	float a2 = pow(vector.dx, 2), b2 = pow(vector.dy, 2),
	      c2 = pow(vector.dz, 2), mag = 0;
	mag = sqrt(a2 + b2 + c2);
	VectorType result = {vector.dx/mag, vector.dy/mag, vector.dz/mag};
	return result;
}

// returns the negative of  a vector
VectorType negativeOfVector(VectorType v) {
	VectorType	res = {-v.dx, -v.dy, -v.dz};
	return res;
}


// returns cross product of two VectorTypes
VectorType getCrossProduct(VectorType a, VectorType b) {
	VectorType res;
	res.dx = a.dy*b.dz - a.dz*b.dy;
	res.dy = a.dz*b.dx - a.dx*b.dz; 
	res.dz = a.dx*b.dy - b.dx*a.dy;
	return res;
}

RayType createRay(PointType p, VectorType a) {
	VectorType dir = sum(p, negativeOfVector(a));
	dir = getUnitVector(dir);
	RayType res = {p.x, p.y, p.z, -dir.dx, -dir.dy, -dir.dz};
	return res;
}

vector<vector<RayType>> getRays(Image im) {
	int w = im.width, h = im.height;

	vector<vector<RayType>> res(h, (vector<RayType>){});

	VectorType	delCH =
			multiplyScalar(sum(im.window.ur, negativeOfVector(im.window.ul)),
			    1.0/(2.0*(float)w));
	VectorType	delCV =
	        multiplyScalar(sum(im.window.ll, negativeOfVector(im.window.ul)),
	            1.0/(2.0*(float)h));

	VectorType	delCenter = sum(im.window.ul, sum(delCH, delCV));

	VectorType	delH =
	        multiplyScalar(sum(im.window.ur, negativeOfVector(im.window.ul)),
	            1.0/((float)w));
	VectorType	delV = 
	        multiplyScalar(sum(im.window.ll, negativeOfVector(im.window.ul)),
	            1.0/((float)h));

	for (int i=0; i<h; i++) {
		for (int j=0; j<w; j++) {
			VectorType v = sum(delCenter, sum(multiplyScalar(delV, i),
			                    multiplyScalar(delH, j)));
			RayType	ray = createRay(im.eye, v);
			res[i].push_back(ray);
		}
	}

	return res;
}

float getSphereIntersectionDistance(RayType	ray, SphereType sphere) {
	float A = 1.0;
	float B = 2.0 *(ray.dx*(ray.x-sphere.x) +  ray.dy*(ray.y-sphere.y)
			+ ray.dz*(ray.z-sphere.z));
	float C = (ray.x-sphere.x)*(ray.x-sphere.x) 
			+ (ray.y-sphere.y)*(ray.y-sphere.y)
	    	+ (ray.z-sphere.z)*(ray.z-sphere.z) 
	        - sphere.r*sphere.r;

	float modulo = B*B - 4.0*A*C;
	if (modulo<0.0) return FLT_MAX;
	else if (modulo==0.0) return -B/(2.0*A);
	else {
		float first, second;
		first = (-B+sqrt(modulo))/(2.0*A);
		second = (-B-sqrt(modulo))/(2.0*A);

		//if (modulo>0.0) cout<<t1<<" "<<t2<<endl;

		if (first>0.0 && second>0.0) return min(first, second);
		else if (first>0.0) return first;
		else if(second>0.0) return second;
		else return FLT_MAX;
	}
}

ColorType shadeRay(Image im, int objId) {
	return im.spheres[objId].m.c;
}

ColorType traceRay(RayType ray, Image im){
	int objId = -1;
	float minDist = FLT_MAX;
	
	for (int i=0; i<im.spheres.size(); i++){
		float dist = getSphereIntersectionDistance(ray, im.spheres[i]);
		// there is not intersection, try the next sphere
		if (dist == FLT_MAX) continue;
		else if (minDist > dist){
			minDist = dist;
			objId = i;
			cout<<"Found Intersection"<<endl;
		}
	}
	// return background color if there is no intersection
	return minDist != FLT_MAX ? shadeRay(im, objId) : im.backgroundColor;
}

void initializeImagePlane(Image& im) {
    im.U = getCrossProduct(im.viewDirection, im.upDirection);
    im.U = getUnitVector(im.U);
	im.V = getCrossProduct(im.U, im.viewDirection);
	im.V = getUnitVector(im.V);
	im.W = {-im.viewDirection.dx, -im.viewDirection.dy, -im.viewDirection.dz};
	im.W = getUnitVector(im.W);
	return;
}

void initializeViewingWindow(Image& im) {
	float width = 2*DOW*tan((im.hfov/2.0)*(PI/180.0));
	float height = width / ((float)im.width/(float)im.height);

	ViewingWindow vw;

	VectorType viewOrig = sum(im.eye, multiplyScalar(getUnitVector(im.viewDirection), DOW));

	vw.ul = sum(viewOrig, sum(multiplyScalar(im.U, -width/2.0), multiplyScalar(im.V, height/2.0))); 
	vw.ur = sum(viewOrig, sum(multiplyScalar(im.U, width/2.0), multiplyScalar(im.V, height/2.0))); 
	vw.ll = sum(viewOrig, sum(multiplyScalar(im.U, -width/2.0), multiplyScalar(im.V, -height/2.0))); 
	vw.lr = sum(viewOrig, sum(multiplyScalar(im.U, width/2.0), multiplyScalar(im.V, -height/2.0))); 

	im.window = vw;
	return;
}


string createFileAndWriteHeader(string fileName, int width, int height) {

	// parse the file name to create output file
	// https://stackoverflow.com/questions/757933/how-do-you-change-the-filename-extension-stored-in-a-string-in-c
	string outputFile = fileName.substr(0, fileName.find_last_of('.'))+".ppm";
	
	ofstream file;
	file.open(outputFile);
	file << "P3" << endl;
	file << "#Ray tracing hw1" << endl;
	file << width << " " << height << endl; 
	file << "255" << endl;
	file.close();

	return outputFile;
}

void drawImage(string fileName, int width, int height, vector<vector<ColorType>>& image) {

	ofstream file;
	file.open(fileName, ios_base::app);
	
	// write pattern to file
	for (int i = 0; i < height; i++) {
	    for (int j = 0; j < width; j++) {
			file <<(int)(255*image[i][j].red) << " "
			  << (int)(255*image[i][j].green) << " "
			  << (int) (255*image[i][j].blue) << endl; 
	    }
	}
	file.close();
	return;
}   

vector<vector<ColorType>> initilialzieImage(Image im) {
	vector<vector<ColorType>> image(im.height, vector<ColorType>(im.width, {0.0, 0.0, 0.0}));
	return image;
}

// parses the file and gets the image size
Image readInput(string fileName) {
	Image image;
	MaterialType material;
  	ifstream infile; infile.open(fileName);
	string LINE;
	unordered_map<string, int> cases = {
		{"bkgcolor", 0}, {"cylinder", 1}, {"eye", 2}, {"hfov", 3},
		{"imsize", 4}, {"mtlcolor", 5}, {"sphere", 6}, {"viewdir", 7},
		{"updir", 8}
	};

	vector<bool> valid(9, false);
	vector<string> tokens;
    while (!infile.eof()) {
    	getline(infile, LINE);
    	tokens = tokenizeLine(LINE, ' ');
    	if (!tokens.size()) continue;
    	switch (cases[tokens[0]]) {
    		case 0: {
				if (tokens.size()!=4 || !checkFloat(tokens[1])
				    || !checkFloat(tokens[2])
				    || !checkFloat(tokens[3])) throw -1;
				float red = stof(tokens[1]), green = stof(tokens[2]),
					blue = stof(tokens[3]);

				if (red<0.0 || red>1.0 || green<0.0 || green>1.0
				    || blue<0.0 || blue>1.0) throw -1;
				valid[cases["bkgcolor"]] = true;
				ColorType color = {red, green, blue};
				image.backgroundColor = color;
				cout << "bkgcolor" << endl;
				break;
			}

    		case 1: {

    		}

    		case 2: {
    			if (tokens.size()!=4 || !checkFloat(tokens[1])
    				|| !checkFloat(tokens[2])
    				|| !checkFloat(tokens[3])) throw -1;
				
				float x = stof(tokens[1]), y = stof(tokens[2]),
				    z = stof(tokens[3]);

				PointType eye = {(float)x, (float)y, (float)z};
				image.eye = eye;
				valid[cases["eye"]] = true;
				cout << "eye" << endl;
				break;
    		}

    		case 3: {
    			if (tokens.size()!=2 || !checkFloat(tokens[1])) throw -1;
				float angle = stof(tokens[1]);
				image.hfov = angle;
				valid[cases["hfov"]] = true;
				cout << "hfov" << endl;
				break;
    		}

    		case 4: {
    			if (tokens.size()!=3 || !checkInt(tokens[1])
    			    || !checkInt(tokens[2])) throw -1;
				int width = stoi(tokens[1]), height = stoi(tokens[2]);
				if (width<=0 || height<=0) throw NULL;
					
				image.width	 = width; image.height = height;
				valid[cases["imsize"]] = true;
				cout << "imsize" << endl;
				break;
    		}

    		case 5: {
    			if (tokens.size()!=4 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) 
    			    || !checkFloat(tokens[3])) throw -1;
				float red = stof(tokens[1]), green = stof(tokens[2]),
					blue = stof(tokens[3]);

				if (red<0.0 || red>1.0 || green<0.0 || green>1.0
				    || blue<0.0 || blue>1.0) throw -1;
				valid[cases["mtlcolor"]] = true;
				ColorType mtlcolor = {(float)red, (float)green, (float)blue};
				material = {mtlcolor};
				cout << "mtlcolor" << endl;
				break;	
    		}

    		case 6: {
    			if (tokens.size()!=5 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4])) throw -1;
				
				float x = stof(tokens[1]), y = stof(tokens[2]),
					z = stof(tokens[3]);
				
				float radius = stof(tokens[4]);

				if (radius<=0 || !valid[cases["mtlcolor"]]) throw -1;
				valid[cases["sphere"]] = true;

				SphereType sphere = {(float)x, (float)y, (float)z,
				    (float)radius, material};
				image.spheres.push_back(sphere);
				cout << "sphere" << endl;
				break;
    		}

    		case 7: {
    			if (tokens.size()!=4 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2])
    			    || !checkFloat(tokens[3])) throw -1;
				
				float x = stof(tokens[1]), y = stof(tokens[2]),
				    z = stof(tokens[3]);
				    // convert to unit lenght ? 
				VectorType viewdir = {(float)x, (float)y, (float)z};
				image.viewDirection = viewdir;
				valid[cases["viewdir"]] = true;
				cout << "viewdir" << endl;
				break;
    		}

    		case 8: {
    			if (tokens.size()!=4 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2])
    			    || !checkFloat(tokens[3])) throw -1;
					float x = stof(tokens[1]), y = stof(tokens[2]),
					    z = stof(tokens[3]);

					VectorType up = {(float)x, (float)y, (float)z};
					image.upDirection = up;
					valid[cases["updir"]] = true;
					cout << "updir" << endl;
					break;	
    		}

    		default:
    			throw -1;
    	}
    }

    // TODO check if all parameters are set or not, if not throw exception

    // return the initialized image parameters
   	return image;
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
    
    try {
    	Image im = readInput(argv[1]);
    	vector<vector<ColorType>> image = initilialzieImage(im);
    	initializeImagePlane(im);
    	initializeViewingWindow(im);
    	vector<vector<RayType>> rays = getRays(im);

    	for (int i=0; i<im.height; i++) {
    		for (int j=0; j<im.width; j++) {
    			image[i][j] = traceRay(rays[i][j], im);
    		}
    	}

    	string file = createFileAndWriteHeader(argv[1], im.width, im.height);
    	drawImage(file, im.width, im.height, image);
    } catch (int e){
    	cout << "something is wrong" << endl;
    }
	return 0;
}



