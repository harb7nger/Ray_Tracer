#include <climits>
#include <float.h>
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

// dummy point
#define PD {(float)FLT_MAX, (float)FLT_MAX, (float)FLT_MAX}

// dummy intersection
#define DI {INT_MAX, (float)FLT_MAX, PD, PD}

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

// struct to store information about a point
struct PointType {
	float x, y, z;
};

// struct to store intersection information
struct TriIntType {
	int objId;
	float dist;
	PointType intPt, bcc;
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
	// to store texture information
	int t;

};

// struct to store information about a texture 
struct TexturePoint {
	float u, v;
};

// struct to store texture 
struct TextureType {
	int width, height;
	std::vector<std::vector<ColorType>> list;
};


// struct to store information about a vector along a direction
struct VectorType {
	float dx, dy, dz;
};

// struct to store information related to the face of a triangle
struct FaceType {
	// vertices of a face
    PointType v1, v2, v3;
    // texture coordinates 
	TexturePoint vt1, vt2, vt3;
	// vertex normals
	VectorType vn1, vn2, vn3;
	// material properties
    MaterialType m;
    // point type and texture index
	int type, t;
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

/* struct to define the image at large
 * all 1 indexed vectors are initialized with dummies
 */
struct Image {
	int width, height;
	VectorType viewDirection, upDirection, U, V, W;
	PointType eye;
	ViewingWindow window;
	// horizontal field of view
	float hfov;
	// defining the background color of the image
	ColorType backgroundColor;
	// vector to store all the spheres in the image
	std::vector<SphereType> spheres;	
	// vector to store all light sources
	std::vector<LightType> lights;	
	// vector to store all vertices
	std::vector<PointType> vertices{{}};
	// adding first value to take care of 1 indexing
	std::vector<FaceType> faces{{}};
	// vector to store all faces in the image
	std::vector<TexturePoint> texts{{}};
	// vector to store all texture points
	std::vector<VectorType> norms{{}};
	// vector to store all norms
	std::vector<TextureType> textures{{}};
};
