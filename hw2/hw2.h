#include <vector>

// distance between view origin in viewing window
#define DOW 5

// defining the constant pie
#define PI 3.14159265

// defining the number of test rays
#define SHADOWTESTCOUNT 121

// jitter
#define JITTER 0.6

// tolerance for shadows
#define EPI 0.061
// struct to store color at a point in the image plane
struct ColorType {
    float red;
    float green;
    float blue;
};

// material type as per phong illumination equation 
struct MaterialType {
	// albedo and specular constants
	ColorType alb, spec;
	// weights 
    float ka, kd, ks; 
	// fall off
	float n;
};

// struct to store information about a ray
struct RayType {
	// origin of the line
	float x, y, z;
	// point minus the origin, the direction
	float dx, dy, dz;
};

struct LightType {
	// the point or direction
	float x, y, z;
	// point or direction light
	bool pointLight, atten;
    // color
	ColorType c;
	// attenuation constants
	float c1, c2, c3;
};

// struct to store information about a sphere
struct SphereType {
	// the center of the sphere
	float x, y, z;
	// radius of the sphere
	float r;
	// for now the shere is of solid color
	MaterialType m;
};

// struct to store information about a point
struct PointType {
	float x, y, z;
};

// struct to store information about a vector along a direction
struct VectorType {
	float dx, dy, dz;
};

// struct to store depth cue related information
struct DepthCue {
    ColorType c;
	float amax, amin, dmax, dmin;
};

/* struct to define the viewing window
 * defines the following parameters
 * ul: upper left corner of the viewing window
 * ur: upper right corner of the viewing window
 * ll: lower left corner of the vieweing window
 * lr: lower right corner of the viewing window
 */
struct ViewingWindow {
	VectorType ul, ur, ll, lr;
};

// struct to define the image at large
struct Image {
	int width, height;
	VectorType viewDirection, upDirection, U, V, W;
	PointType eye;
	ViewingWindow window;
	DepthCue depthQ;
	// horizontal field of view
	float hfov;
	// defining the background color of the image
	ColorType backgroundColor;
	bool depthQFlag;
	// vector to store all the spheres in the image
	std::vector<SphereType> spheres;	
	// vector to store all light sources
	std::vector<LightType> lights;	
};
