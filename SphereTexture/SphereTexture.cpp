//#define _CRT_SECURE_NO_WARNINGS 1
//#include <vector>
//#include <algorithm>
//#include <cmath>
//
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
//
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//
//#define PI 3.1415
//
//class Vector3 {
//public:
//    explicit Vector3(double x = 0, double y = 0, double z = 0) {
//        coords[0] = x;
//        coords[1] = y;
//        coords[2] = z;
//    };
//    double operator[](int i) const { return coords[i]; };
//    double& operator[](int i) { return coords[i]; };
//    double normSquared() {
//        return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
//    }
//    Vector3 normalize() {
//        double norm = sqrt(normSquared());
//        return Vector3(coords[0] / norm, coords[1] / norm, coords[2] / norm);
//    }
//private:
//    double coords[3];
//};
//
//Vector3 operator+(const Vector3& a, const Vector3& b) {
//    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
//}
//
//Vector3 operator-(const Vector3& a, const Vector3& b) {
//    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
//}
//
//Vector3 operator*(double a, const Vector3& b) {
//    return Vector3(a * b[0], a * b[1], a * b[2]);
//}
//
//Vector3 operator*(const Vector3& a, double b) {
//    return Vector3(a[0] * b, a[1] * b, a[2] * b);
//}
//
//Vector3 operator/(const Vector3& a, double b) {
//    return Vector3(a[0] / b, a[1] / b, a[2] / b);
//}
//
//Vector3 operator*(const Vector3& a, const Vector3& b) {
//    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
//}
//
//double dot(const Vector3& a, const Vector3& b) {
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//class Ray {
//public:
//    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
//    }
//    Vector3 C, u;
//};
//
//class Sphere {
//public:
//    Sphere(const Vector3& O, double R, const Vector3& albedo, bool isMirror = false) : O(O), R(R), albedo(albedo), isMirror(isMirror) {
//    }
//    bool intersection(const Ray& r, Vector3& P, Vector3& N, double& t) {
//        double a = 1;
//        double b = 2 * dot(r.u, r.C - O);
//        double c = (r.C - O).normSquared() - R * R;
//        double delta = b * b - 4 * a * c;
//        if (delta < 0) return false;
//
//        double sqrtDelta = sqrt(delta);
//        double t2 = (-b + sqrtDelta) / (2 * a);
//        if (t2 < 0) return false;
//
//        double t1 = (-b - sqrtDelta) / (2 * a);
//        if (t1 > 0) t = t1;
//        else t = t2;
//
//        P = r.C + t * r.u;
//        N = (P - O).normalize();
//
//        return true;
//    };
//    Vector3 O;
//    double R;
//    Vector3 albedo;
//    bool isMirror;
//};
//class Scene {
//public:
//    Scene() {};
//    bool intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& isMirror) {
//        t = 1E9;
//        bool intersected = false;
//        for (int i = 0; i < objects.size(); i++) {
//            Vector3 Pobj, Nobj;
//            double tobj;
//            if (objects[i].intersection(r, Pobj, Nobj, tobj) && tobj < t) {
//                intersected = true;
//                t = tobj;
//                P = Pobj;
//                N = Nobj;
//                albedo = objects[i].albedo;
//                isMirror = objects[i].isMirror;
//            }
//        }
//        return intersected;
//    }
//    Vector3 getColour(const Ray& r, int depth) {
//        if (depth > 5) return Vector3(0., 0., 0.);
//        else {
//            double t;
//            bool isMirror;
//            Vector3 P, N, albedo;
//            Vector3 color(0, 0, 0);
//            if (intersection(r, P, N, albedo, t, isMirror)) {
//                if (isMirror) {
//                    Vector3 reflectedDirection = r.u - 2 * dot(r.u, N) * N;
//                    Ray reflectedRay(P + 0.00001 * N, reflectedDirection);
//                    return getColour(reflectedRay, depth + 1);
//                }
//                else {
//                    double PLNorm = sqrt((lightPosition - P).normSquared());
//                    Vector3 Pshadow, Nshadow, albedoshadow;
//                    double tshadow;
//                    bool isMirrorshadow;
//                    Ray shadowRay(P + 0.001 * N, (lightPosition - P) / PLNorm);
//                    if (intersection(shadowRay, Pshadow, Nshadow, albedoshadow, tshadow, isMirrorshadow) && tshadow < PLNorm) {
//                        color = Vector3(0., 0., 0.);
//                    }
//                    else {
//                        double maxDotProduct = std::max(0., dot(N, (lightPosition - P) / PLNorm));
//                        color = lightIntensity / (4 * PI * PLNorm * PLNorm) * albedo / PI * maxDotProduct;
//                    }
//                }
//            }
//            return color;
//        }
//    }
//
//    std::vector<Sphere> objects;
//    Vector3 lightPosition;
//    Vector3 lightIntensity;
//};
//
//int main() {
//    int W = 512;
//    int H = 512;
//    Vector3 cameraPosition(0, 0, 55);
//    Scene scene;
//
//    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0), true);
//    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(1.0, 0.0, 0.0)); // red
//    Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
//    Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
//    Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
//    Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
//    Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown
//
//    scene.objects.push_back(S1);
//    scene.objects.push_back(Sleft);
//    scene.objects.push_back(Sright);
//    scene.objects.push_back(Stop);
//    scene.objects.push_back(Sbottom);
//    scene.objects.push_back(Sback);
//    scene.objects.push_back(Sfront);
//
//    double fov = 60 * PI / 180;
//    scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
//    scene.lightPosition = Vector3(-10, 20, 40);
//
//    std::vector<unsigned char> image(W * H * 3, 0);
//    for (int i = 0; i < H; i++) {
//        for (int j = 0; j < W; j++) {
//
//            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
//            u = u.normalize();
//
//            Ray r(cameraPosition, u);
//            Vector3 color = scene.getColour(r, 0);
//
//            color[0] = pow(color[0], 0.45);  // gamma correction
//            color[1] = pow(color[1], 0.45);
//            color[2] = pow(color[2], 0.45);
//
//            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
//            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
//            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
//        }
//    }
//    stbi_write_png("MirrorSphere.png", W, H, 3, &image[0], 0);
//
//    return 0;
//
//}



#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415


class Vector3 { //Definition of a class for 3D vectors
public:
	explicit Vector3(double x = 0, double y = 0, double z = 0) {
		coordinates[0] = x;
		coordinates[1] = y;
		coordinates[2] = z;
	};
	double operator[](int i) const { return coordinates[i]; };
	double& operator[](int i) { return coordinates[i]; };
	double NormSquared() { //Squared norm of a vector
		return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
	}
	Vector3 Normalize() { //Normalize a vector
		double norm = sqrt(NormSquared());
		return Vector3(coordinates[0] / norm, coordinates[1] / norm, coordinates[2] / norm);
	}
private:
	double coordinates[3];
};

Vector3 operator+(const Vector3& a, const Vector3& b) { //Vector addition
	return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector3 operator-(const Vector3& a, const Vector3& b) { //Vector subtraction
	return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector3 operator-(const Vector3& a) { //Inverse of a vector
	return Vector3(-a[0], -a[1], -a[2]);
}

Vector3 operator*(double a, const Vector3& b) { //Multiplication of a vector by a scalar
	return Vector3(a * b[0], a * b[1], a * b[2]);
}

Vector3 operator*(const Vector3& a, double b) { //Idem
	return Vector3(a[0] * b, a[1] * b, a[2] * b);
}

Vector3 operator/(const Vector3& a, double b) { //Division of a vector by a scalar
	return Vector3(a[0] / b, a[1] / b, a[2] / b);
}

Vector3 operator*(const Vector3& a, const Vector3& b) {
	return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
double Dot(const Vector3& a, const Vector3& b) { //Dot product between two vectors
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Ray { //Class for rays (characterized by their origin and direction)
public:
	Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {
	}
	Vector3 origin, direction;
};

class Sphere { //Class for spheres (adding an attribute to determine if it is a mirror or not, and another for transparency)
public:
	Sphere(const Vector3& center, double radius, const Vector3& albedo, bool isMirror = false, bool isTransparent = false) : center(center), radius(radius), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {
	}
	bool intersection(const Ray& ray, Vector3& point, Vector3& normal, double& t) {
		// Solve at² + bt + c = 0
		double a = 1;
		double b = 2 * Dot(ray.direction, ray.origin - center);
		double c = (ray.origin - center).NormSquared() - radius * radius;
		double delta = b * b - 4 * a * c;
		if (delta < 0) return false; //no solution, no intersection
		double sqrtDelta = sqrt(delta);
		double t2 = (-b + sqrtDelta) / (2 * a);
		if (t2 < 0) return false; //no positive solution, no intersection

		double t1 = (-b - sqrtDelta) / (2 * a);
		if (t1 > 0) //t is the smallest positive value between t1 and t2
			t = t1;

		else
			t = t2;

		point = ray.origin + t * ray.direction; //Intersection point
		normal = (point - center).Normalize(); //Normal at intersection point

		return true;

	};
	Vector3 center;
	double radius;
	Vector3 albedo;
	bool isMirror, isTransparent;
};

class Scene { //Class for the scene, which can contain multiple spheres
public:
	Scene() {};
	bool intersection(const Ray& ray, Vector3& point, Vector3& normal, Vector3& albedo, double& t, bool& isMirror, bool& isTransparent) {
		t = 1E9;
		bool intersected = false;
		for (int i = 0; i < objects.size(); i++) { //for each sphere in the scene
			Vector3 pointObj, normalObj;
			double tObj;
			if (objects[i].intersection(ray, pointObj, normalObj, tObj) && tObj < t) { //if there is a closer intersection than the existing ones, then we consider this one and not the others
				intersected = true;
				t = tObj;
				point = pointObj;
				normal = normalObj;
				albedo = objects[i].albedo;
				isMirror = objects[i].isMirror;
				isTransparent = objects[i].isTransparent;
			}

		}
		return intersected;
	}
	Vector3 getPixelColor(const Ray& ray, int bounce) {
		if (bounce > 5) {
			return Vector3(0., 0., 0.); //if we exceed the maximum number of bounces, return black
		}
		else {
			double t;
			bool isMirror, isTransparent;
			Vector3 point, normal, albedo, color(0., 0., 0.); //initialize color to black
			if (intersection(ray, point, normal, albedo, t, isMirror, isTransparent)) { //if there is an intersection with an object in the scene
				if (isMirror) { //if it's a mirror object
					Vector3 reflectedDirection = ray.direction - 2 * Dot(ray.direction, normal) * normal; //direction of reflection
					Ray reflectedRay(point + 0.0001 * normal, reflectedDirection); //reflected ray, starting from the intersection point and going in the reflected direction
					color = getPixelColor(reflectedRay, bounce + 1);
				}
				else {
					if (isTransparent) { //if it's a transparent object
						double n1 = 1., n2 = 1.4; //refractive indices
						Vector3 n2Normal = normal;
						if (Dot(ray.direction, normal) > 0) { //if we're exiting the sphere, swap n and the normal
							std::swap(n1, n2);
							n2Normal = -normal;
						}
						double angle = 1 - n1 * n1 / (n2 * n2) * (1 - Dot(ray.direction, n2Normal) * Dot(ray.direction, n2Normal));
						if (angle < 0) { //if the angle is less than 0, there is reflection
							Vector3 reflectedDirection = ray.direction - 2 * Dot(ray.direction, normal) * normal; //direction of reflection
							Ray reflectedRay(point + 0.0001 * normal, reflectedDirection); //reflected ray, starting from the intersection point and going in the reflected direction
							color = getPixelColor(reflectedRay, bounce + 1);
						}
						else {
							Vector3 tangentComponent = n1 / n2 * (ray.direction - Dot(ray.direction, n2Normal) * n2Normal); //tangential component of the refracted direction
							Vector3 normalComponent = -sqrt(angle) * n2Normal; //normal component of the refracted direction
							Vector3 refractedDirection = tangentComponent + normalComponent; //refracted direction
							Ray refractedRay(point - 0.0001 * n2Normal, refractedDirection); //refracted ray, starting from the intersection point and going in the refracted direction
							color = getPixelColor(refractedRay, bounce + 1);
						}
					}
					else { //if it's a diffuse object
						double normLP = sqrt((lightPosition - point).NormSquared()); //norm of LP
						Vector3 pointShadow, normalShadow, albedoShadow;
						double tShadow;
						bool isMirrorShadow, isTransparentShadow;
						Ray shadowRay(point + 0.001 * normal, (lightPosition - point) / normLP); //ray starting from the intersection point and directed towards the light source
						if (intersection(shadowRay, pointShadow, normalShadow, albedoShadow, tShadow, isMirrorShadow, isTransparentShadow) && tShadow < normLP) { //if there is an intersection between the intersection point and the light source
							color = Vector3(0., 0., 0.); //not illuminated = shadow
						}
						else {
							double prov = std::max(0., Dot(normal, (lightPosition - point) / normLP));
							color = lightIntensity / (4 * M_PI * normLP * normLP) * albedo / M_PI * prov;
						}
					}
				}
			}
			return color;
		}
	}

	std::vector<Sphere> objects;
	Vector3 lightPosition;
	Vector3 lightIntensity;
};

int main() {
	int W = 512;
	int H = 512;
	Vector3 cameraPosition(0, 0, 55);
	Scene scene;
	Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0), false, true);
	Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(1.0, 0.0, 0.0)); // red
	Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
	Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
	Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
	Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
	Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown

	scene.objects.push_back(S1);
	scene.objects.push_back(Sleft);
	scene.objects.push_back(Sright);
	scene.objects.push_back(Stop);
	scene.objects.push_back(Sbottom);
	scene.objects.push_back(Sback);
	scene.objects.push_back(Sfront);

	double fov = 60 * M_PI / 180;
	scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
	scene.lightPosition = Vector3(-10, 20, 40);

	std::vector<unsigned char> image(W * H * 3, 0);
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
			u = u.Normalize();

			Ray r(cameraPosition, u);
			Vector3 color = scene.getPixelColor(r, 0);

			color[0] = pow(color[0], 0.45);  // gamma correction
			color[1] = pow(color[1], 0.45);
			color[2] = pow(color[2], 0.45);

			image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
			image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
			image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
		}
	}
	stbi_write_png("TransparentSphere.png", W, H, 3, &image[0], 0);

	return 0;

}