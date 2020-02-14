#include <vector>

// distance between view origin in viewing window
#define DOW 5
// defining constant pie
#define PI 3.14159265

// struct to store color at a point in the image plane
struct ColorType {
    float red;
    float green;
    float blue;
};

// defines the material type for the project
struct MaterialType {
	ColorType c;
};

// struct to store information about a ray
struct RayType {
	// origin of the line
	float x, y, z;
	// point minus the origin
	float dx, dy, dz;
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

// struct to store information about a cylinder 
struct CylinderType {
	// the center of the sphere
	float x, y, z;
	// radius of the sphere
	float r;
	// for now the shere is of solid color
	MaterialType m;
};

// struct to define the viewing window
struct ViewingWindow {
	VectorType ul, ur, ll, lr;
};

// struct to define the image at large
struct Image {
	int width, height;
	VectorType viewDirection, upDirection, U, V, W;
	PointType eye;
	ViewingWindow window;
	float hfov;
	ColorType backgroundColor;
	std::vector<SphereType> spheres;	
};