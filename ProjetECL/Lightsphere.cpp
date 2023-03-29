//#define _CRT_SECURE_NO_WARNINGS
//#include <iostream>
//#include <vector>
//#include <cmath>
//
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
//
//
//const double PI = 3.14159265;
//
//class Vector3 {
//public:
//    explicit Vector3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
//    double lengthSquared() const { return x * x + y * y + z * z; }
//    Vector3 normalize() const {
//        double len = std::sqrt(lengthSquared());
//        return Vector3(x / len, y / len, z / len);
//    }
//
//    double x, y, z;
//};
//
//Vector3 operator+(const Vector3& a, const Vector3& b) {
//    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
//}
//
//Vector3 operator-(const Vector3& a, const Vector3& b) {
//    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
//}
//
//Vector3 operator*(double scalar, const Vector3& v) {
//    return Vector3(scalar * v.x, scalar * v.y, scalar * v.z);
//}
//
//double dotProduct(const Vector3& a, const Vector3& b) {
//    return a.x * b.x + a.y * b.y + a.z * b.z;
//}
//
//class Ray {
//public:
//    Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {}
//    Vector3 origin, direction;
//};
//
//class Sphere {
//public:
//    Sphere(const Vector3& center, double radius) : center(center), radius(radius) {}
//
//    bool intersects(const Ray& ray) {
//        double a = 1;
//        double b = 2 * dotProduct(ray.direction, ray.origin - center);
//        double c = (ray.origin - center).lengthSquared() - radius * radius;
//        double delta = b * b - 4 * a * c;
//        return delta >= 0;
//    }
//
//    Vector3 center;
//    double radius;
//};
//
//int main() {
//    int width = 512;
//    int height = 512;
//
//    Vector3 cameraPos(0, 0, 55);
//    Vector3 sphereCenter(0, 0, 0);
//    double sphereRadius = 10;
//    Sphere sphere(sphereCenter, sphereRadius);
//    double LightIntensity = 100000000;
//    double fov = 60 * PI / 180;
//
//    std::vector<unsigned char> image(width * height * 3, 0);
//
//    for (int i = 0; i < height; ++i) {
//        for (int j = 0; j < width; ++j) {
//            Vector3 direction(j - width / 2, i - height / 2, -width / (2 * std::tan(fov / 2)));
//            direction = direction.normalize();
//            Ray ray(cameraPos, direction);
//            Vector3 color(0, 0, 0);
//
//            if (sphere.intersects(ray)) {
//                color = Vector3(255, 255, 255);
//            }
//
//            image[(i * width + j) * 3 + 0] = static_cast<unsigned char>(color.x);
//            image[(i * width + j) * 3 + 1] = static_cast<unsigned char>(color.y);
//            image[(i * width + j) * 3 + 2] = static_cast<unsigned char>(color.z);
//        }
//    }
//
//    stbi_write_png("whitesphere.png", width, height, 3, &image[0], 0);
//
//    return 0;
//}

////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////


#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415

class Vector3 { //Definition of a class for vectors
public:
    explicit Vector3(double x = 0, double y = 0, double z = 0) {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    };
    double operator[](int i) const { return coordinates[i]; };
    double& operator[](int i) { return coordinates[i]; };
    double NormeAuCarre() { //square norm
        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
    }
    Vector3 normalize() { //to normalize a vector
        double norm = sqrt(NormeAuCarre());
        return Vector3(coordinates[0] / norm, coordinates[1] / norm, coordinates[2] / norm);
    }
private:
    double coordinates[3];
};

Vector3 operator+(const Vector3& a, const Vector3& b) { //sum of two vectors
    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector3 operator-(const Vector3& a, const Vector3& b) { //difference of two vectors
    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a double
    return Vector3(a * b[0], a * b[1], a * b[2]);
}

Vector3 operator*(const Vector3& a, double b) { //idem
    return Vector3(a[0] * b, a[1] * b, a[2] * b);
}

Vector3 operator/(const Vector3& a, double b) { //division of a vector by a double
    return Vector3(a[0] / b, a[1] / b, a[2] / b);
}


double dot(const Vector3& a, const Vector3& b) { //dot product between two vectors
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Ray { //Class for rays (characterized by their center and their direction)
public:
    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
    }
    Vector3 C, u;
};

class Sphere { //Class for spheres
public:
    Sphere(const Vector3& O, double R) : O(O), R(R) {
    }
    bool intersection(const Ray& r, Vector3& P, Vector3& N) {
        // solves a*t² + b*t + c = 0
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c = (r.C - O).NormeAuCarre() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; //no solution, therefore no intersection

        double sqrt_delta = sqrt(delta);
        double t2 = (-b + sqrt_delta) / (2 * a);
        if (t2 < 0) return false; //no positive solution, therefore no intersection

        double t1 = (-b - sqrt_delta) / (2 * a);
        double t; //smallest positive value between t1 and t2
        if (t1 > 0)
            t = t1;
        else
            t = t2;

        P = r.C + t * r.u; // Intersection point
        N = (P - O).normalize(); //Normal at intersection point

        return true;

    };
    Vector3 O;
    double R;
};

int main() {
    int W = 512;
    int H = 512;

    Vector3 C(0, 0, 55);
    Vector3 O(0, 0, 0);
    double R = 10; //radius of the sphere
    Sphere S(O, R);
    double fov = 60 * M_PI / 180;
    double I = 1E7; //light intensity
    Vector3 SphereColor(1, 1, 0);
    Vector3 L(-10, 20, 40); //light source

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
            u = u.normalize();
            Ray r(C, u);
            Vector3 P, N;
            Vector3 color(0, 0, 0);
            if (S.intersection(r, P, N)) { //in case of intersection with the sphere, we calculate the color of the pixel
                double normPL = sqrt((L - P).NormeAuCarre()); //PL norm
                double prov = std::max(0., dot(N, (L - P) / normPL));
                color = I / (4 * M_PI * normPL * normPL) * SphereColor / M_PI * prov;
            }

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("coloredsphere.png", W, H, 3, &image[0], 0);

    return 0;
}
