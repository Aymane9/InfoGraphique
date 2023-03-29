#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <random>

#define M_PI 3.1415

static std::default_random_engine engine(10);
static std::uniform_real_distribution<double> uniform(0, 1);

class Vector3 { //Defining a class for 3D vectors
public:
    explicit Vector3(double x = 0, double y = 0, double z = 0) {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    };
    double operator[](int i) const { return coordinates[i]; };
    double& operator[](int i) { return coordinates[i]; };
    double NormSquared() { //squared norm
        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
    }
    Vector3 Normalize() { //normalize a vector
        double norm = sqrt(NormSquared());
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

Vector3 operator-(const Vector3& a) { //inverse of a vector
    return Vector3(-a[0], -a[1], -a[2]);
}

Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a scalar
    return Vector3(a * b[0], a * b[1], a * b[2]);
}

Vector3 operator*(const Vector3& a, double b) { //same as above
    return Vector3(a[0] * b, a[1] * b, a[2] * b);
}

Vector3 operator/(const Vector3& a, double b) { //division of a vector by a scalar
    return Vector3(a[0] / b, a[1] / b, a[2] / b);
}

Vector3 operator*(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector3 CrossProduct(const Vector3& a, const Vector3& b) { //cross product of two vectors
    return Vector3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double DotProduct(const Vector3& a, const Vector3& b) { //dot product of two vectors
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector3 TermByTermProduct(const Vector3& a, const Vector3& b) { //term-by-term product of two vectors
    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector3 RandomInUnitSphere(const Vector3& N) { //Returns a vector zN+xT1+y*T2, where (N,T1,T2) is a coordinate system, and x,y,z are random variables following a cosine probability distribution (the radius is more likely to be close to N)
    double u1 = uniform(engine); // random number between 0 and 1
    double u2 = uniform(engine);
    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
    double z = sqrt(u2);
    Vector3 T1;
    if (N[0] < N[1] && N[0] < N[2]) {
        T1 = Vector3(0, N[2], -N[1]);
    }
    else {
        if (N[1] < N[2] && N[1] < N[0]) {
            T1 = Vector3(N[2], 0, -N[0]);
        }
        else {
            T1 = Vector3(N[1], -N[0], 0);
        }
    }
    T1 = T1.Normalize();
    Vector3 T2 = CrossProduct(N, T1);
    return x * T1 + y * T2 + z * N;
}

class Ray { //Class for rays (characterized by their origin and direction)
public:
    Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {
    }
    Vector3 origin, direction;
};

class Sphere { //Class for spheres (added an attribute to indicate if it's a mirror or not, and another for transparency)
public:
    Sphere(const Vector3& center, double radius, const Vector3& albedo, bool isMirror = false, bool isTransparent = false) : center(center), radius(radius), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {
    }
    bool Intersect(const Ray& ray, Vector3& P, Vector3& N, double& t) {
        // solve a*t² + b*t + c = 0
        double a = 1;
        double b = 2 * DotProduct(ray.direction, ray.origin - center);
        double c = (ray.origin - center).NormSquared() - radius * radius;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; //no solution, no intersection

        double sqrtDelta = sqrt(delta);
        double t2 = (-b + sqrtDelta) / (2 * a);
        if (t2 < 0) return false; //no positive solution, no intersection

        double t1 = (-b - sqrtDelta) / (2 * a);
        if (t1 > 0) // t is the smallest positive value between t1 and t2
            t = t1;
        else
            t = t2;

        P = ray.origin + t * ray.direction; //Intersection point
        N = (P - center).Normalize(); //Normal at the intersection point

        return true;
    }
    Vector3 center;
    double radius;
    Vector3 albedo;
    bool isMirror, isTransparent;
};


class Scene {
public:
    Scene() {};
    bool Intersect(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& mirror, bool& transp, int& id) {
        t = 1E9;
        bool intersected = false;
        for (int i = 0; i < objects.size(); i++) {
            Vector3 P_object, N_object;
            double t_object;
            if (objects[i].Intersect(r, P_object, N_object, t_object) && t_object < t) {
                intersected = true;
                t = t_object;
                P = P_object;
                N = N_object;
                albedo = objects[i].albedo;
                mirror = objects[i].isMirror;
                transp = objects[i].isTransparent;
                id = i;
            }

        }
        return intersected;
    }

    Vector3 GetColor(const Ray& r, int bounce, bool lastDiffuse) {
        if (bounce > 5) return Vector3(0., 0., 0.);
        else {
            double t;
            bool mirror, transp;
            Vector3 P, N, albedo;
            Vector3 color(0, 0, 0);
            int id;
            if (Intersect(r, P, N, albedo, t, mirror, transp, id)) {
                if (id == 0) {
                    if (bounce == 0 || !lastDiffuse)
                        return lightIntensity / (4 * M_PI * M_PI * objects[0].radius * objects[0].radius);
                    else
                        return Vector3(0, 0, 0);
                }

                if (mirror) {
                    Vector3 reflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N;
                    Ray reflectedRay(P + 0.00001 * N, reflectedDirection);
                    return GetColor(reflectedRay, bounce + 1, false);
                }
                else {
                    if (transp) {
                        double n1 = 1., n2 = 1.4;
                        Vector3 N2 = N;
                        if (DotProduct(r.direction, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - DotProduct(r.direction, N2) * DotProduct(r.direction, N2));
                        if (angle < 0) {
                            Vector3 reflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N;
                            Ray reflectedRay(P + 0.00001 * N, reflectedDirection);
                            return GetColor(reflectedRay, bounce + 1, false);
                        }
                        Vector3 T_t = n1 / n2 * (r.direction - DotProduct(r.direction, N2) * N2);
                        Vector3 T_n = -sqrt(angle) * N2;
                        Vector3 refractedDirection = T_t + T_n;
                        Ray refractedRay(P - 0.0001 * N2, refractedDirection);
                        return GetColor(refractedRay, bounce + 1, false);
                    }
                    else {
                        Vector3 w = RandomInUnitSphere((P - lightPosition).Normalize());
                        Vector3 xp = w * objects[0].radius + objects[0].center;
                        Vector3 Pxp = xp - P;
                        double normPxp = sqrt(Pxp.NormSquared());
                        Pxp = Pxp.Normalize();
                        Vector3 P_shadow, N_shadow, albedo_shadow;
                        double t_shadow;
                        bool mirror_shadow, transp_shadow;
                        int id_shadow;
                        Ray shadowRay(P + 0.00001 * N, Pxp);
                        if (Intersect(shadowRay, P_shadow, N_shadow, albedo_shadow, t_shadow, mirror_shadow, transp_shadow, id_shadow) && t_shadow < normPxp - 0.0001) {
                            color = Vector3(0., 0., 0.);
                        }
                        else {
                            color = lightIntensity / (4 * M_PI * M_PI * objects[0].radius * objects[0].radius) * albedo / M_PI * std::max(0., DotProduct(N, Pxp)) * std::max(0., -DotProduct(w, Pxp)) / (normPxp * normPxp) / (std::max(0., -DotProduct((lightPosition - P).Normalize(), w)) / (M_PI * objects[0].radius * objects[0].radius));

                            Vector3 wi = RandomInUnitSphere(N);
                            Ray randomRay(P + 0.00001 * N, wi);
                            color = color + TermByTermProduct(albedo, GetColor(randomRay, bounce + 1, true));

                        }
                    }
                    return color;
                }
            }
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

    Sphere Slumiere(scene.lightPosition, 5, Vector3(1., 0.3, 0.2));
    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0));
    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(1.0, 0.0, 0.0)); // red
    Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
    Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
    Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
    Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
    Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown

    scene.objects.push_back(Slumiere);
    scene.objects.push_back(S1);
    scene.objects.push_back(Sleft);
    scene.objects.push_back(Sright);
    scene.objects.push_back(Stop);
    scene.objects.push_back(Sbottom);
    scene.objects.push_back(Sback);
    scene.objects.push_back(Sfront);


    double fov = 70 * M_PI / 180;
    scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
    
    scene.lightPosition = Vector3(0, 20, 40);
    int nmbRays = 3;

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {



            Vector3 color(0, 0, 0);
            for (int k = 0; k < nmbRays; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
                Vector3 u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2)));
                u = u.Normalize();

                Ray r(cameraPosition, u);
                color = color + scene.GetColor(r, 0,false);
        }
            color = color / nmbRays;


            color[0] = std::pow(color[0], 0.45);//application du gamma
            color[1] = std::pow(color[1], 0.45);
            color[2] = std::pow(color[2], 0.45);


            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("SoftShadows.png", W, H, 3, &image[0], 0);

    return 0;
}


