#include "hw2.h"
#include <algorithm>
#include <ctype.h>
#include <cmath>
#include <float.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>

using namespace std;

// checks if a given string represents a float value or not
bool checkFloat (string word) {
	int decimalCount = 0;
	if (word[0] == '-' || word[0] == '+') {
	    word = word.substr(1, word.size()-1);
	}

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

// tokenize a line of text, used for parsing the input file
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

// returns the vector a->b given two points a, b 
VectorType getVector(PointType a, PointType b) {
	VectorType result;
	result.dx = b.x - a.x;
	result.dy = b.y - a.y;
	result.dz = b.z - a.z;
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

// display point
void displayPoint(PointType a) {
	cout << a.x << "," << a.y << "," << a.z << endl;
}
	

// returns the magnitude of a vector
float getMagnitude(VectorType v) {
    float a2 = pow(v.dx, 2), b2 = pow(v.dy, 2),
          c2 = pow(v.dz, 2), mag = 0;
    mag = sqrt(a2 + b2 + c2);
    return mag;
}

// returns a normalized VectorType
VectorType getUnitVector(VectorType vector) {
	float a2 = pow(vector.dx, 2), b2 = pow(vector.dy, 2),
	      c2 = pow(vector.dz, 2), mag = 0;
	mag = sqrt(a2 + b2 + c2);
	VectorType result = {vector.dx/mag, vector.dy/mag, vector.dz/mag};
	return result;
}

// returns the negative of a vector
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

// returns a scalar dot product between two VectorType
float getDotProduct(VectorType a, VectorType b) {
	float res = a.dx*b.dx + a.dy*b.dy + a.dz*b.dz;
	return res;
}

// returns a scaled color
ColorType scaleColor(ColorType color, float scale) {
	float red = scale*color.red, green = scale*color.green, blue = scale*color.blue;
    ColorType res = {red, green, blue};
	return res;
}

// adds the resultant of two ColorType
ColorType addColors(ColorType a, ColorType b) {
    float red = a.red+b.red, green = a.green+b.green, blue = a.blue+b.blue;
	ColorType res = {red, green, blue};
	return res;
}

// returns the dot product of two colors, usually intensities and color
ColorType dotProduct(ColorType a, ColorType b) {
    float red = a.red*b.red, green = a.green*b.green, blue = a.blue*b.blue;
    ColorType res = {red, green, blue};
	return res;
}

// returns rounded values of color intensities
ColorType clampColor(ColorType a) {
    float red = a.red > 1.0 ? 1.0 : a.red,
		  green = a.green > 1.0 ? 1.0 : a.green,
	      blue = a.blue > 1.0 ? 1.0 : a.blue;
	ColorType res = {red, green, blue};
	return res;
}


// creates a ray given a point and vector
RayType createRay(PointType p, VectorType a) {
	VectorType dir = sum(p, negativeOfVector(a));
	dir = getUnitVector(dir);
	RayType res = {p.x, p.y, p.z, -dir.dx, -dir.dy, -dir.dz};
	return res;
}

// 
RayType getRay(PointType a, VectorType v) {
	RayType res = {a.x, a.y, a.z, v.dx, v.dy, v.dz};
	return res;
}

// given ray and distance returns the point
PointType getPoint(RayType ray, float dist) {
	float x = ray.x + dist*ray.dx, y = ray.y + dist*ray.dy, z = ray.z + dist*ray.dz;
	PointType res = {x, y, z};
	return res;
}

// returns the attenuation factor for a given light source
float getAttnFactor(LightType light, float dist) {
	return (float)1.0/(light.c1+light.c2*dist+light.c3*dist*dist);
}

// creates all the rays that will be used to probe through 
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

// to get intersection distance for each ray and sphere pair
float getSphereIntersectionDistance(RayType	ray, SphereType sphere) {
	float A = 1.0;
	float B = 2.0*(ray.dx*(ray.x-sphere.x) + ray.dy*(ray.y-sphere.y)
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

// returns the shade of a particular ray
ColorType shadeRay(Image im, int objId, RayType ray, float dist) {
    SphereType sphere = im.spheres[objId];
    ColorType amb = scaleColor(sphere.m.alb, sphere.m.ka), diff = {0.0, 0.0, 0.0},
			  spec = {0.0, 0.0, 0.0};	
	PointType intPt = getPoint(ray, dist), center = {sphere.x, sphere.y, sphere.z}; 
	VectorType surfNorm = getVector(center, intPt), zero = {0.0, 0.0, 0.0}; 
	surfNorm = getUnitVector(surfNorm);
	VectorType L;
	// iterate through all light sources
	for (LightType& light: im.lights) {
		float shadowFlag = 1.0;
	    if (light.pointLight) {// point light case
			PointType l = {light.x, light.y, light.z};
			L = getVector(intPt, l);
		} else { // directional light case
      	 	L = {light.x, light.y, light.z};	
		}
		// diffusion related terms of phong equation
		float mag = getMagnitude(L);
    	L = getUnitVector(L);
		float ln = getDotProduct(L, surfNorm);
		ln = ln < 0.0 ? 0.0:ln;
		ColorType lightDiff = scaleColor(sphere.m.alb, ln*sphere.m.kd);
		lightDiff = dotProduct(lightDiff, light.c); // multiplying intensities
		
        // specular related terms of phong equation	
		VectorType V = getVector(intPt, im.eye);
	    V = getUnitVector(V);
		VectorType H = sum(V, L);
    	H = getUnitVector(H);	
		float nh = getDotProduct(H, surfNorm);
		nh = nh < 0.0 ? 0.0:nh;
		nh = pow(nh, sphere.m.n); 
		ColorType lightSpec = scaleColor(sphere.m.spec, nh*sphere.m.ks);
		lightSpec = dotProduct(lightSpec, light.c); // multiplying intensities

		// to calculate the shadow component
		for (int t=0; t<SHADOWTESTCOUNT; t++) {
            float j_x = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
	   		float j_y = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
			float j_z = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
			j_x = j_x*JITTER; j_y = j_y*JITTER; j_z = j_z*JITTER;
			//displayPoint({j_x, j_y, j_z});
    
		   RayType ray = createRay(intPt, {light.x+j_x, light.y+j_y, light.z+j_z});
		   if (!light.pointLight) { // if the light is directional
		       ray = getRay(intPt, L);
		   }
		   float minDist = FLT_MAX;
		   for (int i=0; i<im.spheres.size(); i++) { // check for intersection with all spheres
		   		if (i == objId) continue; // not to check for sphere where the intersection lies on
		   		float dist = getSphereIntersectionDistance(ray, im.spheres[i]);
		   		if (dist == FLT_MAX) continue;	
		   		if (minDist > dist && dist > EPI) minDist = dist;
		   }

		   if (!light.pointLight && minDist!=FLT_MAX) shadowFlag += 0.0;
		   else if (light.pointLight && minDist!=FLT_MAX) {
		   		if (minDist < mag) {shadowFlag += 0.0;}
		   		else {shadowFlag += 1.0;}
		   } else if (minDist==FLT_MAX) {shadowFlag += 1.0;}
		}

	   	if (!light.pointLight && shadowFlag<1.0) {shadowFlag = 0.0;}
		else {shadowFlag = shadowFlag/(float)SHADOWTESTCOUNT;}
		
		float atF = 1; // attenuation factor
		if (light.atten) atF = getAttnFactor(light, mag); 

		diff = addColors(diff, scaleColor(lightDiff, atF*shadowFlag)); // add the effect of this light
 		spec = addColors(spec, scaleColor(lightSpec, atF*shadowFlag));// add the effect of this light 
	}
	
    // add up all the diffusion and spec terms for all light sources
	ColorType res = addColors(diff, amb); res = addColors(res, spec);
	res = clampColor(res); // clamp intensities to prevent overflow
	return res;
}

// to trace the given ray
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
			//cout<<"Itersection"<<endl;
		}
	}
	// return background color if there is no intersection
	return minDist != FLT_MAX ? shadeRay(im, objId, ray, minDist) : im.backgroundColor;
}

// initialize the image plane, all the vectors such as U, V & W
void initializeImagePlane(Image& im) {
    im.U = getCrossProduct(im.viewDirection, im.upDirection);
	if (getMagnitude(im.U) == 0.0) throw 1; // if u is of zero notify user
    im.U = getUnitVector(im.U);
	im.V = getCrossProduct(im.U, im.viewDirection);
	im.V = getUnitVector(im.V);
	im.W = {-im.viewDirection.dx, -im.viewDirection.dy, -im.viewDirection.dz};
	im.W = getUnitVector(im.W);
	return;
}

// initialize the viewing window
void initializeViewingWindow(Image& im) {
	float width = 2*DOW*tan((im.hfov/2.0)*(PI/180.0));
	float height = width / ((float)im.width/(float)im.height);

	ViewingWindow vw;

	VectorType viewOrig = sum(im.eye, multiplyScalar(getUnitVector(im.viewDirection), DOW));

	// defining the corners of the viewing window
	vw.ul = sum(viewOrig, sum(multiplyScalar(im.U, -width/2.0), multiplyScalar(im.V, height/2.0))); 
	vw.ur = sum(viewOrig, sum(multiplyScalar(im.U, width/2.0), multiplyScalar(im.V, height/2.0))); 
	vw.ll = sum(viewOrig, sum(multiplyScalar(im.U, -width/2.0), multiplyScalar(im.V, -height/2.0))); 
	vw.lr = sum(viewOrig, sum(multiplyScalar(im.U, width/2.0), multiplyScalar(im.V, -height/2.0))); 

	im.window = vw;
	return;
}

// creates a new file which ends with .ppm from the input file name and puts in the header information
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

// draws the image into the ppm file
void drawImage(string fileName, int width, int height, vector<vector<ColorType>>& image) {

	ofstream file;
	file.open(fileName, ios_base::app);
	// write image to file
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

// creates a vector which holds all the colors to be written the image
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
		{"bkgcolor", 0}, {"eye", 1}, {"hfov", 2},
		{"imsize", 3}, {"light", 4}, {"mtlcolor", 5}, {"sphere", 6}, {"viewdir", 7},
		{"updir", 8}, {"attlight", 9} 
	};

	vector<bool> validArgs(9, false);
	vector<string> tokens;
    while (!infile.eof()) {
    	getline(infile, LINE);
    	cout << LINE << endl;
    	tokens = tokenizeLine(LINE, ' ');
    	if (!tokens.size()) continue;
		// initialize inputs for the image through a switch case
    	switch (cases[tokens[0]]) {
    		case 0: {
				if (tokens.size()!=4 || !checkFloat(tokens[1])
				    || !checkFloat(tokens[2])
				    || !checkFloat(tokens[3])) throw -1;
				float red = stof(tokens[1]), green = stof(tokens[2]),
					blue = stof(tokens[3]);

				if (red<0.0 || red>1.0 || green<0.0 || green>1.0
				    || blue<0.0 || blue>1.0) throw -1;

				validArgs[cases["bkgcolor"]] = true;
				ColorType color = {red, green, blue};
				image.backgroundColor = color;
				cout << "bkgcolor" << endl;
				break;
			}

    		case 1: {
    			if (tokens.size()!=4 || !checkFloat(tokens[1])
    				|| !checkFloat(tokens[2])
    				|| !checkFloat(tokens[3])) throw -1;
				
				float x = stof(tokens[1]), y = stof(tokens[2]),
				    z = stof(tokens[3]);

				PointType eye = {(float)x, (float)y, (float)z};
				image.eye = eye;
				validArgs[cases["eye"]] = true;
				cout << "eye" << endl;
				break;
    		}

    		case 2: {
    			if (tokens.size()!=2 || !checkFloat(tokens[1])) throw -1;
				float angle = stof(tokens[1]);
				image.hfov = angle;
				validArgs[cases["hfov"]] = true;
				cout << "hfov" << endl;
				break;
    		}

    		case 3: {
    			if (tokens.size()!=3 || !checkInt(tokens[1])
    			    || !checkInt(tokens[2])) throw -1;
				int width = stoi(tokens[1]), height = stoi(tokens[2]);
				if (width<=0 || height<=0) throw -1;
					
				image.width	 = width; image.height = height;
				validArgs[cases["imsize"]] = true;
				cout << "imsize" << endl;
				break;
    		}

			case 4: {
				if (tokens.size()!=8 || !checkFloat(tokens[1])
					|| !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4]) || !checkFloat(tokens[5])
    			    || !checkFloat(tokens[6]) || !checkFloat(tokens[7])) throw -1;

				float x = stof(tokens[1]), y = stof(tokens[2]), z = stof(tokens[3]),
					  w = stof(tokens[4]), r = stof(tokens[5]), g= stof(tokens[6]),
                      b = stof(tokens[7]);
				if ((w!=1.0 && w!=0) || r<0.0 || r>1.0 || g<0.0 || g>1.0 
					|| b<0.0 || b>1.0) throw -1; 
				ColorType color = {(float)r, (float)g, (float)b};
                if (w == 1) { // for point lights
			     	LightType light 	
						= {(float)x, (float)y, (float)z, true, false, color, 0, 0, 0};	
					image.lights.push_back(light);
				} else { // for directional lights
				   	VectorType v = {(float)x, (float)y, (float)z};
				    v = negativeOfVector(v);
					v = getUnitVector(v);
				    LightType light 
						= {v.dx, v.dy, v.dz, false, false, color, 0, 0, 0};	
					image.lights.push_back(light);
				}
				validArgs[cases["light"]] = true;
				cout << "lights" << endl;
				break;
			}

    		case 5: {
    			if (tokens.size()!=11 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4]) || !checkFloat(tokens[5])
    			    || !checkFloat(tokens[6]) || !checkFloat(tokens[7])
    			    || !checkFloat(tokens[8]) || !checkFloat(tokens[9])
    			    || !checkFloat(tokens[10])) throw -1;

				float oDr = stof(tokens[1]), oDg = stof(tokens[2]),
					  oDb = stof(tokens[3]), oSr = stof(tokens[4]),
					  oSg = stof(tokens[5]), oSb = stof(tokens[6]),
					  ka = stof(tokens[7]), kd = stof(tokens[8]),
					  ks = stof(tokens[9]), n = stof(tokens[10]);
                               

				if (oDr<0.0 || oDr>1.0 || oDg<0.0 || oDg>1.0 || oDb<0.0 || oDb>1.0
				    || oSr<0.0 || oSr>1.0 || oSg<0.0 || oSg>1.0 || oSb<0.0 || oSb>1.0
				    || ka<0.0 || ka>1.0 || kd<0.0 || kd>1.0 || ks<0.0 || ks>1.0) throw -1;

				validArgs[cases["mtlcolor"]] = true;
				ColorType oD = {(float)oDr, (float)oDg, (float)oDb},
				          oS = {(float)oSr, (float)oSg, (float)oSb};
				
				material = {oD, oS, (float)ka, (float)kd, (float)ks, (float)n};
				// cout << "mtlcolor" << endl;
				break;	
    		}

    		case 6: {
    			if (tokens.size()!=5 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4])) throw -1;
			        cout << "here1" << endl;	
				float x = stof(tokens[1]), y = stof(tokens[2]),
					z = stof(tokens[3]);
				
				float radius = stof(tokens[4]);

				if (radius<=0 || !validArgs[cases["mtlcolor"]]) throw -1;
				validArgs[cases["sphere"]] = true;

			        cout << "here2" << endl;	
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
				VectorType viewdir = {(float)x, (float)y, (float)z};
				image.viewDirection = viewdir;
				validArgs[cases["viewdir"]] = true;
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
					validArgs[cases["updir"]] = true;
					cout << "updir" << endl;
					break;	
    		}

			case 9: {
				if (tokens.size()!=11 || !checkFloat(tokens[1])
					|| !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4]) || !checkFloat(tokens[5])
    			    || !checkFloat(tokens[6]) || !checkFloat(tokens[7])
    			    || !checkFloat(tokens[8]) || !checkFloat(tokens[9])
    			    || !checkFloat(tokens[10])) throw -1;

				float x = stof(tokens[1]), y = stof(tokens[2]), z = stof(tokens[3]),
					  w = stof(tokens[4]), r = stof(tokens[5]), g= stof(tokens[6]),
                      b = stof(tokens[7]), c1 = stof(tokens[8]),
					  c2 = stof(tokens[9]), c3 = stof(tokens[10]);
				if (w!=1.0 || r<0.0 || r>1.0 || g<0.0 || g>1.0 
					|| b<0.0 || b>1.0) throw -1; 
				ColorType color = {(float)r, (float)g, (float)b};
			    LightType light 
						= {(float)x, (float)y, (float)z, true, true, color, c1, c2, c3};	
				image.lights.push_back(light);
				validArgs[cases["light"]] = true;
				// cout << "lights" << endl;
				break;
			}

    		default: {
    			throw -1;
    		}
   		}
    }

	// check if all the parameters are set or throw error
	// TODO remove the exception for cylinder once implemented
	for (const bool& valid : validArgs) if (!valid) throw 0;	
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

		//use the traced ray to get the image
    	for (int i=0; i<im.height; i++) {
    		for (int j=0; j<im.width; j++) {
    			image[i][j] = traceRay(rays[i][j], im);
    		}
    	}

    	string file = createFileAndWriteHeader(argv[1], im.width, im.height);
    	drawImage(file, im.width, im.height, image);
    } catch (int e) {
		switch(e) {
			case -1: 
				cout << "please verify input" << endl;
				break;
			case 0:
				cout << "all parameters not specified, please check" << endl;
				break;
			case 1:
				cout << "U vector is zero, please check view and up direction" << endl; 
				break;
		}
    }

	return 0;
}



