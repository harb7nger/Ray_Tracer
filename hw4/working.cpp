#include "working.h"
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

bool checkInt (string word) { // checking only for +ve numbers
	for (char& c: word) if (!isdigit(c)) return false;
	return true;
}

// tokenize a line of text, used for parsing the input file
vector<string> tokenizeLine(string line, char delim) {
	stringstream ss(line);
	vector<string> words;
	string word;
	while (getline(ss, word, delim)) {
  		if (word != "") words.push_back(word);
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

// display texture coordinates
void displayTextCoor(TexturePoint t) {
	cout << t.u << "," << t.v << endl;
}

// display face
void displayFace(FaceType f) {
	cout << "vertices" << endl;
	displayPoint(f.v1);
	displayPoint(f.v2);
	displayPoint(f.v3);
	cout << "textures" << endl;
	displayTextCoor(f.vt1);
	displayTextCoor(f.vt2);
	displayTextCoor(f.vt3);
	cout << "normals" << endl;
	displayVector(f.vn1);
	displayVector(f.vn2);
	displayVector(f.vn3);
	cout << "type:" << f.type << endl;
	return;
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

// returns a scalar dot product between VectorType and a PointType
float getDotProduct(VectorType a, PointType b) {
	float res = a.dx*b.x + a.dy*b.y + a.dz*b.z;
	return res;
}

// get reflected ray given incident ray
RayType getReflectedRay(RayType incidentRay, VectorType norm, PointType interPt) {
	VectorType ray = {-incidentRay.dx, -incidentRay.dy, -incidentRay.dz};
	ray = getUnitVector(ray);
	float dotProd = getDotProduct(norm, ray);
	VectorType prod = multiplyScalar(norm, 2*dotProd);
    RayType res = {interPt.x, interPt.y, interPt.z,
		  	   	   prod.dx-ray.dx, prod.dy-ray.dy, prod.dz-ray.dz};	
	return res;
}

// get fresnel refelctance coefficient
float getFresnelRC(RayType incidentRay, VectorType norm, float fo) {
	VectorType ray = {-incidentRay.dx, -incidentRay.dy, -incidentRay.dz};
	ray = getUnitVector(ray);
	float dotProd = getDotProduct(norm, ray);
	float exp = pow((1.0-dotProd),5.0);
	float fr = fo + (1.0-fo)*exp;
	return fr;
}

// returns if colors are equal
bool areEqual(ColorType a, ColorType b) {		
	if (a.red == b.red
		&& a.green == b.green
		&& a.blue == b.blue) return true;
	return false;
}

// returns if two points are equal
bool areEqual(PointType a, PointType b) {
	if (a.x == b.x
	    && a.y == b.y
		&& a.z == b.z) return true;
	return false;
}

bool areEqual(TriIntType a, TriIntType b) {
	if (a.objId == b.objId
		&& a.dist == b.dist
		&& areEqual(a.intPt, b.intPt)
		&& areEqual(a.bcc, b.bcc)) return true;
	return false;
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


// displays color
void displayColor(ColorType c) {
	cout << "red:" << c.red
		 << "green:" << c.green
		 << "blue:" << c.blue << endl;
}

// returns resultant color when considering two colors
ColorType getResultant(ColorType a, ColorType b, Image& im) {
	if (areEqual(a, im.backgroundColor) && areEqual(b, im.backgroundColor)) {
		return im.backgroundColor;
	} else if (!areEqual(a, im.backgroundColor) && areEqual(b, im.backgroundColor)) {
		return a;
	} else if (areEqual(a, im.backgroundColor) && !areEqual(b, im.backgroundColor)) {
		return b;
	}  
    return clampColor(addColors(a, b));
}


// creates a ray given a point and vector
RayType createRay(PointType p, VectorType a) {
	VectorType dir = getUnitVector(a);
	RayType res = {p.x, p.y, p.z, dir.dx, dir.dy, dir.dz};
	return res;
}

// creates a ray given a point and vector
RayType createRay(PointType p, PointType a) {
	VectorType dir = getVector(p, a);
	dir = getUnitVector(dir);
	RayType res = {p.x, p.y, p.z, dir.dx, dir.dy, dir.dz};
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
	PointType res = {(float)x, (float)y, (float)z};
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
			VectorType v = sum(delCenter, 
				sum(multiplyScalar(delV, i),multiplyScalar(delH, j)));
			PointType pt = {v.dx, v.dy, v.dz};
			RayType	ray = createRay(im.eye, pt);
			res[i].push_back(ray);
		}
	}

	return res;
}

// to get intersection distance for each ray and sphere pair
float getSphereIntersectionDistance(RayType	ray, SphereType sphere) {
	float A = 1.0;
	float B = 2.0*(ray.dx*(ray.x-sphere.x) + ray.dy*(ray.y-sphere.y) + ray.dz*(ray.z-sphere.z));
	float C = (ray.x-sphere.x)*(ray.x-sphere.x) + (ray.y-sphere.y)*(ray.y-sphere.y) + (ray.z-sphere.z)*(ray.z-sphere.z) - sphere.r*sphere.r;

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

// get the distance for a ray to intersect a plane given by a face
float getPlaneIntersectionDistance(RayType ray, FaceType face) {
	PointType p0 = face.v1, p1 = face.v2, p2 = face.v3,
			  rayOri = {(float)ray.x, (float)ray.y, (float)ray.z};
	VectorType e1 = getVector(p0, p1), e2 = getVector(p0, p2);
	VectorType n = getCrossProduct(e1, e2),
			   rayDir = {(float)ray.dx, (float)ray.dy, (float)ray.dz};
    float den = getDotProduct(n, rayDir);
	if (den == 0.0) return FLT_MAX;
	// optimization, checking denominator before forming the equation
	float D = -getDotProduct(n, p0);
	float num = -(getDotProduct(n, rayOri)+D);
	float t = num/den;
	return t<0.0? FLT_MAX:t; // if t is neg implies the point is behind the ray orig
}

// to get barycentri coordinates
PointType getBarycentricCoord(PointType p, FaceType face) {
	// intersection point wrt plane
    //displayPoint(p);
    PointType p0 = face.v1, p1 = face.v2, p2 = face.v3;
    VectorType e1 = getVector(p0, p1), e2 = getVector(p0, p2),
			   e3 = getVector(p1, p), e4 = getVector(p2, p);
    
    float A = getMagnitude(getCrossProduct(e1, e2)),
		  a = getMagnitude(getCrossProduct(e3, e4)),
		  b = getMagnitude(getCrossProduct(e4, e2)),
		  c = getMagnitude(getCrossProduct(e1, e3));
	float alpha = a/A, beta = b/A, gamma = c/A;
			
    // cout << "alpha" << alpha << "beta" << beta << "gamma" << gamma << endl;
	float sum = alpha+beta+gamma, one = (float)1.0;
    if (alpha < 0.0 || alpha > 1.0 // to check if bcc are valid
        || beta < 0.0 || beta > 1.0
        || gamma < 0.0 || gamma > 1.0 || abs(sum - 1) > BTOL) return PD;
    return {alpha, beta, gamma};
}


TriIntType getFaceIntersection(RayType ray, Image& im) {
	float minDist = FLT_MAX;
	int objId = -1;
/*	cout << "printing ray again" << endl;
	displayRay(ray);*/
	PointType bcc = PD, intPt;
	if (im.faces.size() == 1) return DI;
	for (int i=1; i<im.faces.size(); i++) {
		float dist = getPlaneIntersectionDistance(ray, im.faces[i]); 
		if (dist == FLT_MAX) continue;
		PointType pt = getPoint(ray, dist);
		PointType b = getBarycentricCoord(pt, im.faces[i]);
		if (minDist > dist && dist > EPI && !areEqual(b, PD)) {
			minDist = dist;
			objId = i;
			bcc = b;
			intPt = pt;
		}
	}
	return {objId, minDist, intPt, bcc};
}

// returns the shade of a particular ray
ColorType shadeRay(int objType, Image& im, int objId, PointType intPt, PointType bcc, RayType ray, int recDepth) {
	// initialize 
	ColorType amb, odlam, oslam,  
			  diff = {0.0, 0.0, 0.0}, spec = {0.0, 0.0, 0.0};	
	float kd = 0, ks = 0, n = 0, fo = 0;
    VectorType surfNorm;

	if (objType == 1) { // handle everything as a face 
		FaceType face = im.faces[objId];
		kd = face.m.kd; ks = face.m.ks; n = face.m.n;
		amb = scaleColor(face.m.alb, face.m.ka);
		odlam = face.m.alb; oslam = face.m.spec;
		fo = face.m.fo;
		PointType p0 = face.v1, p1 = face.v2, p2 = face.v3;	
		VectorType e1 = getVector(p0, p1), e2 = getVector(p0, p2);
		if (face.type == 2 || face.type == 3) {
			surfNorm = sum(sum(multiplyScalar(face.vn1, bcc.x), multiplyScalar(face.vn2, bcc.y)),
				multiplyScalar(face.vn3, bcc.z));
		} else {
			surfNorm = getCrossProduct(e1, e2);
		}

		if (face.t != -1) {
			// cout << "BANNER" << endl;
			TextureType& texture =  im.textures[face.t];
			float hor = bcc.x*face.vt1.u + bcc.y*face.vt2.u + bcc.z*face.vt3.u;
			float vert = bcc.x*face.vt1.v + bcc.y*face.vt2.v + bcc.z*face.vt3.v;

			int i = (int)(vert*(float)(texture.height-1));
			int j = (int)(hor*(float)(texture.width-1));

			odlam = texture.list[i][j];
		}
			
	} else { // handle everything as a sphere
		SphereType sphere = im.spheres[objId];
		kd = sphere.m.kd; ks = sphere.m.ks; n = sphere.m.n;
    	amb = scaleColor(sphere.m.alb, sphere.m.ka);
		odlam = sphere.m.alb; oslam = sphere.m.spec;
		fo = sphere.m.fo;
	 	PointType center = {sphere.x, sphere.y, sphere.z}; 
		surfNorm = getVector(center, intPt); 

		if (sphere.t != -1) {
			TextureType& texture = im.textures[sphere.t];
			VectorType norm = getUnitVector(surfNorm);
			float vert = acos(norm.dz)/PI;
			float theta = atan2(norm.dy, norm.dz);
			theta = theta < 0.0 ? (float) (theta + 2*PI):(float)theta;
			float hor = theta/(2*PI);

			int i = (int)(vert*(float)(texture.height-1));
			int j = (int)(hor*(float)(texture.width-1));

			odlam = texture.list[i][j];
		}
	}
	/*displayVector(surfNorm);
    displayColor(amb);*/

    surfNorm = getUnitVector(surfNorm);
	VectorType L;
	VectorType V = getVector(intPt, im.eye);
	float distEyeIntPt = getMagnitude(V);
    V = getUnitVector(V);

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
		//cout << "ln:" << ln << endl;
		ln = ln < 0.0 ? 0.0:ln;
		ColorType lightDiff = scaleColor(odlam, ln*kd);
		lightDiff = dotProduct(lightDiff, light.c); // multiplying intensities
		
        // specular related terms of phong equation	
		VectorType H = sum(V, L);
    	H = getUnitVector(H);	
		float nh = getDotProduct(H, surfNorm);
		//cout << "nh:" << nh << endl;
		nh = nh < 0.0 ? 0.0:nh;
		nh = pow(nh, n); 
		ColorType lightSpec = scaleColor(oslam, nh*ks);
		lightSpec = dotProduct(lightSpec, light.c); // multiplying intensities
		//cout << "before" << endl;
	    //displayColor(lightDiff);
		//displayColor(lightSpec);
		//cout << "after" << endl;
		// to calculate the shadow component
		for (int t=0; t<SHADOWTESTCOUNT; t++) {
            float j_x = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
	   		float j_y = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
			float j_z = static_cast<float>(rand())/static_cast<float>(RAND_MAX);
			j_x = j_x*JITTER; j_y = j_y*JITTER; j_z = j_z*JITTER;
			//displayPoint({j_x, j_y, j_z});
    	   PointType endPt = {light.x+j_x, light.y+j_y, light.z+j_z};
		   RayType ray = createRay(intPt, endPt);
		   if (!light.pointLight) { // if the light is directional
		       ray = getRay(intPt, L);
		   }

		   float minDist = FLT_MAX;
		   for (int i=0; i<im.spheres.size(); i++) { // check for intersection with all spheres
		   		if (i == objId && objType == 0) continue; // not to check for sphere where the intersection lies on
		   		float dist = getSphereIntersectionDistance(ray, im.spheres[i]);
		   		if (dist == FLT_MAX) continue;	
		   		if (minDist > dist && dist > EPI) minDist = dist;
		   }
		   
		   TriIntType intr = getFaceIntersection(ray, im);
		   if (intr.dist < minDist && (objType == 0 || objId != intr.objId)) minDist = intr.dist;

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
		//displayColor(diff);
		//displayColor(spec);
	}

	ColorType refColor = DC; 
	if (recDepth+1 < MAXRECUR) {
        float frc = getFresnelRC(ray, surfNorm, fo);
       	RayType reflRay = getReflectedRay(ray, surfNorm, intPt);    
	    ColorType refColor = traceRay(ray, im, recDepth+1);
	    //displayColor(refColor);	
        refColor = scaleColor(refColor, frc);
	}
	
    // add up all the diffusion and spec terms for all light sources
	ColorType res = addColors(diff, amb); res = addColors(res, spec);
	res = addColors(refColor, res);
	res = clampColor(res); // clamp intensities to prevent overflow
	//displayColor(res);
	return res;
}

pair<float, int> getSphereIntersection(RayType ray, Image& im) {
	float minDist = FLT_MAX;
	int objId = -1;
	for (int i=0; i<im.spheres.size(); i++){
		float dist = getSphereIntersectionDistance(ray, im.spheres[i]);
		// there is not intersection, try the next sphere
		if (dist == FLT_MAX) {
			continue;
		}
		else if (minDist > dist && dist > EPI) {
			minDist = dist;
			objId = i;
			//cout<<"Itersection"<<endl;
		}
	}
	return make_pair(minDist, objId);
}


/*
 * traces a ray through all objects described in the scene
 * iterates through all spheres firstly
 * iterates through all triangles secondly
 * if intersection distance of sphere is lesser
 	* return shading as per the sphere
 * if intersection distance of triangle is lesser
	* getBarycentricCoordinates
	* if barycentricCoordinates are invalid, repeat the process for circle 
	* else shade ray as per the triangle
*/
ColorType traceRay(RayType ray, Image& im, int recDepth) {	
	float minDist = FLT_MAX;
    pair<float, int> sphereInt = getSphereIntersection(ray, im);
    // cout<<result.first<<endl;
    TriIntType triInt = getFaceIntersection(ray, im);
	if (sphereInt.first < triInt.dist) {
	     //cout << "after:" << result.first << endl;
	     //cout << "afterS:" << result.second << endl;
	     //cout << (result.first == FLT_MAX) << endl;
		  PointType pt = getPoint(ray, sphereInt.first);
		  if (recDepth > 0) cout << "type" << 0 << "id" << sphereInt.second << endl;
	      return shadeRay(0, im, sphereInt.second, pt, PD, ray, recDepth);
	} else {
	     //cout << "after:" << triInt.dist << endl;
	     if (triInt.dist == FLT_MAX) return im.backgroundColor;
	     //cout << triInt.objId << endl;
		if (recDepth > 0) cout << "type" << 1 << "id" << triInt.objId << endl;
    	return shadeRay(1, im, triInt.objId, triInt.intPt, triInt.bcc, ray, recDepth);
	}
	return im.backgroundColor;
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


// to parse the texture file
TextureType parseTextureFile(string fileName) {
	//cout << "parsing" << endl;
	ifstream file(fileName);
	if (file.is_open()) {
		string line;
		getline(file, line);

		vector<string> headerInfo = tokenizeLine(line, ' ');

		if (!checkInt(headerInfo[1]) || !checkInt(headerInfo[2]) || !checkInt(headerInfo[3])) throw 7;

		int w = stoi(headerInfo[1]), h = stoi(headerInfo[2]), maxClr = stoi(headerInfo[3]);

		int pixel = 0, col = 0, max = w*h;

		vector<ColorType> parsed;

		//cout<<maxClr<<endl;

		float inv = (float)(1.0/maxClr);

		while (pixel < max && getline(file, line)) {
			vector<string> tks = tokenizeLine(line, ' ');
			ColorType c;
			for (int i = 0; i < tks.size(); i++) {
				switch(col%3) {
					case 0: {
						c.red = (float)stoi(tks[i])*inv;
						col++;
						break;
					}

					case 1: {
						c.green = (float)stoi(tks[i])*inv;
						col++;
						break;
					}

					case 2: {
						c.blue = (float)stoi(tks[i])*inv;
						col++;
						pixel++;
						break;
					}
				}
			}
			parsed.push_back(c);
		}

		file.close();

		vector<vector<ColorType>> result;

		pixel = 0;

		for (int i=0; i<h; i++) {
			vector<ColorType> temp;
			for (int j=0; j<w; j++) {
				temp.push_back(parsed[pixel++]);
			}
			result.push_back(temp);
		}

		TextureType texture;
		texture.list = result;
		texture.width = w;
		texture.height = h;

		return texture;
	} else {throw 8;}
	cout << "parsing failed" << endl;
	file.close();
	return {}; // return empty object
}


// parses the file and gets the image size
Image readInput(string fileName) {

	Image image;
	MaterialType material;
	int texture = 0;
  	ifstream infile; infile.open(fileName);
	string LINE;
	unordered_map<string, int> cases = {
		{"bkgcolor", 0}, {"eye", 1}, {"hfov", 2}, {"imsize", 3},
		{"light", 4}, {"mtlcolor", 5}, {"sphere", 6}, {"viewdir", 7},
		{"updir", 8}, {"attlight", 9},
		{"v", 10}, {"vn", 11}, {"vt", 12}, {"f", 13}, {"texture", 14}
	};

	vector<bool> validArgs(14, false);
	vector<string> tokens;
    while (!infile.eof()) {
    	getline(infile, LINE);
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
				// cout << "bkgcolor" << endl;
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
				// cout << "eye" << endl;
				break;
    		}

    		case 2: {
    			if (tokens.size()!=2 || !checkFloat(tokens[1])) throw -1;
				float angle = stof(tokens[1]);
				image.hfov = angle;
				validArgs[cases["hfov"]] = true;
				// cout << "hfov" << endl;
				break;
    		}

    		case 3: {
    			if (tokens.size()!=3 || !checkInt(tokens[1])
    			    || !checkInt(tokens[2])) throw -1;
				int width = stoi(tokens[1]), height = stoi(tokens[2]);
				if (width<=0 || height<=0) throw -1;
					
				image.width	 = width; image.height = height;
				validArgs[cases["imsize"]] = true;
				// cout << "imsize" << endl;
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
				// cout << "lights" << endl;
				break;
			}

    		case 5: {
    			if (tokens.size()!=13 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4]) || !checkFloat(tokens[5])
    			    || !checkFloat(tokens[6]) || !checkFloat(tokens[7])
    			    || !checkFloat(tokens[8]) || !checkFloat(tokens[9])
    			    || !checkFloat(tokens[10]) || !checkFloat(tokens[11])
					|| !checkFloat(tokens[12])) throw -1;


				float oDr = stof(tokens[1]), oDg = stof(tokens[2]),
					  oDb = stof(tokens[3]), oSr = stof(tokens[4]),
					  oSg = stof(tokens[5]), oSb = stof(tokens[6]),
					  ka = stof(tokens[7]), kd = stof(tokens[8]),
					  ks = stof(tokens[9]), n = stof(tokens[10]),
                      alpha = stof(tokens[11]), eta = stof(tokens[12]);         

				if (oDr<0.0 || oDr>1.0 || oDg<0.0 || oDg>1.0 || oDb<0.0 || oDb>1.0
				    || oSr<0.0 || oSr>1.0 || oSg<0.0 || oSg>1.0 || oSb<0.0 || oSb>1.0
				    || ka<0.0 || ka>1.0 || kd<0.0 || kd>1.0 || ks<0.0 || ks>1.0
					|| alpha<0.0 || alpha>1.0 || eta < 1.0) throw -1;

				float fo = pow(((eta-1)/(eta+1)), 2.0); 

				validArgs[cases["mtlcolor"]] = true;
				ColorType oD = {(float)oDr, (float)oDg, (float)oDb},
				          oS = {(float)oSr, (float)oSg, (float)oSb};
				
				material = {oD, oS, (float)ka, (float)kd, (float)ks, (float)n, alpha, fo};
				// cout << "mtlcolor" << endl;
				break;	
    		}

    		case 6: {
    			if (tokens.size()!=5 || !checkFloat(tokens[1])
    			    || !checkFloat(tokens[2]) || !checkFloat(tokens[3])
    			    || !checkFloat(tokens[4])) throw -1;
				float x = stof(tokens[1]), y = stof(tokens[2]),
					z = stof(tokens[3]);
				
				float radius = stof(tokens[4]);

				if (radius<=0 || !validArgs[cases["mtlcolor"]]) throw -1;
				validArgs[cases["sphere"]] = true;
				int t = texture > 0 ? texture : -1;
				SphereType sphere = {(float)x, (float)y, (float)z,
				    (float)radius, material, t};
				image.spheres.push_back(sphere);
				// cout << "sphere" << endl;
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
				// cout << "viewdir" << endl;
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
					// cout << "updir" << endl;
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

			case 10: {
                if (tokens.size()!=4 || !checkFloat(tokens[1])
					|| !checkFloat(tokens[2]) || !checkFloat(tokens[3])) throw 2;

				float x = stof(tokens[1]), y = stof(tokens[2]), z = stof(tokens[3]);
				PointType vertex = {(float)x, (float)y, (float)z};
				image.vertices.push_back(vertex);
				// cout << "vertex" << endl;
				break;
			}

			case 11: {
				if (tokens.size()!=4 || !checkFloat(tokens[1])
					|| !checkFloat(tokens[2]) || !checkFloat(tokens[3])) throw 3;
				float dx = stof(tokens[1]), dy = stof(tokens[2]), dz = stof(tokens[3]);
			 	VectorType normal = {(float)dx, (float)dy, (float)dz};
				image.norms.push_back(normal);
				// cout << "vertex normal" << endl;
				break;
			}

			case 12: {
				if (tokens.size()!=3 || !checkFloat(tokens[1])
					|| !checkFloat(tokens[2])) throw 4;
				float u = stof(tokens[1]), v = stof(tokens[2]);
			 	TexturePoint texturePt = {(float)u, (float)v};
				image.texts.push_back(texturePt);
				// cout << "texture" << endl;
				break;
			}

			case 13: {
                if (tokens.size()!=4) throw -1;

				vector<string> first = tokenizeLine(tokens[1], '/'),
						second = tokenizeLine(tokens[2], '/'),
						third = tokenizeLine(tokens[3], '/');

				// get lenght of each token set
				int fs = first.size(), ss = second.size(), ts = third.size();
				// get the bounds for all vectors
				int vertexSize = image.vertices.size(), normSize = image.norms.size(),
					textSize = image.texts.size();

				// initialize lookups for face as 0s
				int v1 = 0, v2 = 0, v3 = 0, // indices for all attr
					vt1 = 0, vt2 = 0, vt3 = 0, 
					vn1 = 0, vn2 = 0, vn3 = 0;

				int type = 0;
				// count the number of slashes, thats how we know the input type
				int slash = count(tokens[1].begin(), tokens[1].end(), '/');

				// if lengths don't match
				if ((fs!=ss)||(ss!=ts)||(ts!=fs)) throw 5;

				/* performs indirect hashing of sorts
				 * v -> 1
				 * v/vt -> 3
				 * v//vn -> 4
				 * v/vt/vn -> 5
				 * the corresponding types are 0, 1, 2, 3
				*/
				fs = fs+slash;
				// handles differnt pncs of face inputs
				switch(fs) {
					  case 1: {
 			 			v1 = stoi(first[0]), v2 = stoi(second[0]), v3 = stoi(third[0]);					
						type = 0; // stands for a basic face
						break;			 
					  }

					  case 3: {
 			 			v1 = stoi(first[0]), v2 = stoi(second[0]), v3 = stoi(third[0]);					
						vt1 = stoi(first[1]), vt2 = stoi(second[1]), vt3 = stoi(third[1]);
						type = 1; // stands for texture without smooth shading
						break;			 
					  }
					    	
					  case 4: {
 			 			v1 = stoi(first[0]), v2 = stoi(second[0]), v3 = stoi(third[0]);					
						vn1 = stoi(first[1]), vn2 = stoi(second[1]), vn3 = stoi(third[1]);
						type = 2; // stands for face with smooth shading only
						break;			 

					  }

					  case 5: {
 			 			v1 = stoi(first[0]), v2 = stoi(second[0]), v3 = stoi(third[0]);					
						vt1 = stoi(first[1]), vt2 = stoi(second[1]), vt3 = stoi(third[1]);
						vn1 = stoi(first[2]), vn2 = stoi(second[2]), vn3 = stoi(third[2]);
						type = 3; // stands for texture and smooth shading 
						break;			 
					  }
				}
					 
					   
				if (v1 >= vertexSize || v2 >= vertexSize || v3 >= vertexSize) throw 5;
				if (vn1 >= normSize || vn2 >= normSize || vn3 >= normSize) throw 5;

				PointType P1 = image.vertices[v1], P2 = image.vertices[v2],
					  P3 = image.vertices[v3]; // vertices 

				TexturePoint T1 = image.texts[vt1], T2 = image.texts[vt2],
				            T3 = image.texts[vt3]; // texture points 

				VectorType N1 = image.norms[vn1], N2 = image.norms[vn2], 
					   N3 = image.norms[vn3]; // norms
				//cout << type << endl;
				int t = texture > 0 ? texture : -1;
				FaceType face = {P1, P2, P3, T1, T2, T3,
						N1, N2, N3, material, type, t};

				//displayFace(face);

				image.faces.push_back(face);				
				break;
			}

			case 14: {
				if (tokens.size()!=2) throw 6;
				TextureType tex = parseTextureFile(tokens[1]);
				image.textures.push_back(tex);
				// cout << "height" << tex.height << "width" << tex.width << endl; 
				texture++;
				break;	
			}

    		default: {
    			throw -1;
    		}
   		}
    }

	// check if all the parameters are set or throw error
	// for (const bool& valid : validArgs) if (!valid) throw 0;	
    // return the initialized image parameters
       // exit(5);
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
			    
				image[i][j] = traceRay(rays[i][j], im, 0); 
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
			case 2:
				cout << "Invalid vertex entry encountered, please check input" << endl;
				break;
			case 3:
				cout << "Invalid vertex normal entry encountered, please check input" << endl;
				break;
			case 4:
				cout << "Invalid vertex texture entry encountered, please check input" << endl;
				break;
			case 5:
				cout << "Invalid face entry encountered, please check input" << endl;
				break;
			case 6:
				cout << "Texture file parsing failed please check" << endl;
				break;
			case 7:
				cout << "Error occured while parsing header for texture file" << endl;
				break;
			case 8:
				cout << "Error, texture file did not open" << endl;
		}
    }

	return 0;
}
