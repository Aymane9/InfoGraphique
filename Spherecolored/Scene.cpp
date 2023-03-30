////////////////////////////////////////////////////////////////////////////

// First try at a scene

////////////////////////////////////////////////////////////////////////////

//#define _CRT_SECURE_NO_WARNINGS 1
//#include <vector>
//#include <algorithm>
//#include <iostream>
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
//
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//
//#define M_PI 3.1415
//
//class Vector3 { // Definition of a class for vectors
//public:
//    explicit Vector3(double x = 0, double y = 0, double z = 0) {
//        coordinates[0] = x;
//        coordinates[1] = y;
//        coordinates[2] = z;
//    };
//    double operator[](int i) const { return coordinates[i]; };
//    double& operator[](int i) { return coordinates[i]; };
//    double NormSquared() { //norm squared
//        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
//    }
//    Vector3 Normalize() { //to normalize a vector
//        double norm = sqrt(NormSquared());
//        return Vector3(coordinates[0] / norm, coordinates[1] / norm, coordinates[2] / norm);
//    }
//private:
//    double coordinates[3];
//};
//
//Vector3 operator+(const Vector3& a, const Vector3& b) { //sum of two vectors
//    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
//}
//
//Vector3 operator-(const Vector3& a, const Vector3& b) { //difference of two vectors
//    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
//}
//
//Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a double
//    return Vector3(a * b[0], a * b[1], a * b[2]);
//}
//
//Vector3 operator*(const Vector3& a, double b) { //idem
//    return Vector3(a[0] * b, a[1] * b, a[2] * b);
//}
//
//Vector3 operator/(const Vector3& a, double b) { //division of a vector by a double
//    return Vector3(a[0] / b, a[1] / b, a[2] / b);
//}
//
//
//double Dot(const Vector3& a, const Vector3& b) { //dot product between 2 vectors
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//
//class Ray { //Class for rays (characterized by their center and their direction)
//public:
//    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
//    }
//    Vector3 C, u;
//};
//
//class Sphere { //class for spheres
//public:
//    Sphere(const Vector3& O, double R, const Vector3& albedo) : O(O), R(R), albedo(albedo) {
//    }
//    bool Intersection(const Ray& r, Vector3& P, Vector3& N, double& t) {
//        // solves a*t² + b*t + c = 0
//        double a = 1;
//        double b = 2 * Dot(r.u, r.C - O);
//        double c = (r.C - O).NormSquared() - R * R;
//        double delta = b * b - 4 * a * c;
//        if (delta < 0) return false; // no solution, thus no intersection
//
//        double deltaSqrt = sqrt(delta);
//        double t2 = (-b + deltaSqrt) / (2 * a);
//        if (t2 < 0) return false; // no positive solution, thus no intersection
//
//        double t1 = (-b - deltaSqrt) / (2 * a);
//        if (t1 > 0) // t is the smallest positive value between t1 and t2
//            t = t1;
//
//        else
//            t = t2;
//
//        P = r.C + t * r.u; // intersection point
//        N = (P - O).Normalize(); //normal at the intersection point
//
//        return true;
//
//    };
//    Vector3 O;
//    double R;
//    Vector3 albedo;
//};
//
//class Scene { // scene class, which can contain multiple spheres
//public:
//    Scene() {};
//    bool Intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo) {
//        double t = 1E9;
//        bool intersected = false;
//        for (int i = 0; i < objects.size(); i++) { // for each sphere in the scene
//            Vector3 PObject, NObject;
//            double tObject;
//            if (objects[i].Intersection(r, PObject, NObject, tObject) && tObject < t) { // if there's an intersection closer than existing ones, then we take this one into account and not the others
//                intersected = true;
//                t = tObject;
//                P = PObject;
//                N = NObject;
//                albedo = objects[i].albedo;
//            }
//
//        }
//        return intersected;
//    }
//    std::vector<Sphere> objects;
//
//};
//
//int main() {
//    int W = 512;
//    int H = 512;
//
//    Scene scene;
//    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0)); // bright red
//    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(0.8, 0.2, 0.9)); // magenta
//    Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
//    Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
//    Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
//    Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
//    Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown
//    Vector3 C(0, 0, 55);
//
//    scene.objects.push_back(S1);
//    scene.objects.push_back(Sfront);
//    scene.objects.push_back(Sback);
//    scene.objects.push_back(Stop);
//    scene.objects.push_back(Sbottom);
//    scene.objects.push_back(Sright);
//    scene.objects.push_back(Sleft);
//
//
//    double fov = 70 * M_PI / 180;
//    double I = 10000000; //light intensity
//    Vector3 L(0, 20, 40); //light source
//
//    std::vector<unsigned char> image(W * H * 3, 0);
//    for (int i = 0; i < H; i++) {
//        for (int j = 0; j < W; j++) {
//
//            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
//            u = u.Normalize();
//
//            Ray r(C, u);
//            Vector3 P, N, albedo;
//            Vector3 color(0, 0, 0);
//            if (scene.Intersection(r, P, N, albedo)) { // in case of intersection with an object in the scene
//                double normPL = sqrt((L - P).NormSquared()); //norm of PL
//                double prov = std::max(0., Dot(N, (L - P) / normPL));
//                color = I / (4 * M_PI * normPL * normPL) * albedo / M_PI * prov;
//            }
//
//            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
//            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
//            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
//        }
//    }
//    stbi_write_png("scene.png", W, H, 3, &image[0], 0);
//
//    return 0;
//}


////////////////////////////////////////////////////////////////////////////

// Scene with gamma correction

////////////////////////////////////////////////////////////////////////////


//#define _CRT_SECURE_NO_WARNINGS 1
//#include <vector>
//#include <algorithm>
//#include <iostream>
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
//
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//
//#define M_PI 3.1415
//
//class Vector3 { // Definition of a class for vectors
//public:
//    explicit Vector3(double x = 0, double y = 0, double z = 0) {
//        coordinates[0] = x;
//        coordinates[1] = y;
//        coordinates[2] = z;
//    };
//    double operator[](int i) const { return coordinates[i]; };
//    double& operator[](int i) { return coordinates[i]; };
//    double NormSquared() { //norm squared
//        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
//    }
//    Vector3 Normalize() { //to normalize a vector
//        double norm = sqrt(NormSquared());
//        return Vector3(coordinates[0] / norm, coordinates[1] / norm, coordinates[2] / norm);
//    }
//private:
//    double coordinates[3];
//};
//
//Vector3 operator+(const Vector3& a, const Vector3& b) { //sum of two vectors
//    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
//}
//
//Vector3 operator-(const Vector3& a, const Vector3& b) { //difference of two vectors
//    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
//}
//
//Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a double
//    return Vector3(a * b[0], a * b[1], a * b[2]);
//}
//
//Vector3 operator*(const Vector3& a, double b) { //idem
//    return Vector3(a[0] * b, a[1] * b, a[2] * b);
//}
//
//Vector3 operator/(const Vector3& a, double b) { //division of a vector by a double
//    return Vector3(a[0] / b, a[1] / b, a[2] / b);
//}
//
//
//double Dot(const Vector3& a, const Vector3& b) { //dot product between 2 vectors
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//
//class Ray { //Class for rays (characterized by their center and their direction)
//public:
//    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
//    }
//    Vector3 C, u;
//};
//
//class Sphere { //class for spheres
//public:
//    Sphere(const Vector3& O, double R, const Vector3& albedo) : O(O), R(R), albedo(albedo) {
//    }
//    bool Intersection(const Ray& r, Vector3& P, Vector3& N, double& t) {
//        // solves a*t² + b*t + c = 0
//        double a = 1;
//        double b = 2 * Dot(r.u, r.C - O);
//        double c = (r.C - O).NormSquared() - R * R;
//        double delta = b * b - 4 * a * c;
//        if (delta < 0) return false; // no solution, thus no intersection
//
//        double deltaSqrt = sqrt(delta);
//        double t2 = (-b + deltaSqrt) / (2 * a);
//        if (t2 < 0) return false; // no positive solution, thus no intersection
//
//        double t1 = (-b - deltaSqrt) / (2 * a);
//        if (t1 > 0) // t is the smallest positive value between t1 and t2
//            t = t1;
//
//        else
//            t = t2;
//
//        P = r.C + t * r.u; // intersection point
//        N = (P - O).Normalize(); //normal at the intersection point
//
//        return true;
//
//    };
//    Vector3 O;
//    double R;
//    Vector3 albedo;
//};
//
//class Scene { // scene class, which can contain multiple spheres
//public:
//    Scene() {};
//    bool Intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo) {
//        double t = 1000000000;
//        bool intersected = false;
//        for (int i = 0; i < objects.size(); i++) { // for each sphere in the scene
//            Vector3 PObject, NObject;
//            double tObject;
//            if (objects[i].Intersection(r, PObject, NObject, tObject) && tObject < t) { // if there's an intersection closer than existing ones, then we take this one into account and not the others
//                intersected = true;
//                t = tObject;
//                P = PObject;
//                N = NObject;
//                albedo = objects[i].albedo;
//            }
//
//        }
//        return intersected;
//    }
//    std::vector<Sphere> objects;
//
//};
//
//int main() {
//    int W = 512;
//    int H = 512;
//
//    Scene scene;
//    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0)); // bright red
//    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(0.8, 0.2, 0.9)); // magenta
//    Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
//    Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
//    Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
//    Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
//    Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown
//    Vector3 C(0, 0, 55);
//
//    scene.objects.push_back(S1);
//    scene.objects.push_back(Sfront);
//    scene.objects.push_back(Sback);
//    scene.objects.push_back(Stop);
//    scene.objects.push_back(Sbottom);
//    scene.objects.push_back(Sright);
//    scene.objects.push_back(Sleft);
//
//
//    double fov = 70 * M_PI / 180;
//    double I = 10000000000; //light intensity
//    Vector3 L(0, 20, 40); //light source
//
//    std::vector<unsigned char> image(W * H * 3, 0);
//    for (int i = 0; i < H; i++) {
//        for (int j = 0; j < W; j++) {
//
//            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
//            u = u.Normalize();
//
//            Ray r(C, u);
//            Vector3 P, N, albedo;
//            Vector3 color(0, 0, 0);
//            if (scene.Intersection(r, P, N, albedo)) { // in case of intersection with an object in the scene
//                double normPL = sqrt((L - P).NormSquared()); //norm of PL
//                double prov = std::max(0., Dot(N, (L - P) / normPL));
//                color = I / (4 * M_PI * normPL * normPL) * albedo / M_PI * prov;
//            }
//
//
//            color[0] = std::pow(color[0], 0.45);//applying gamma correction
//            color[1] = std::pow(color[1], 0.45);
//            color[2] = std::pow(color[2], 0.45);
//
//            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
//            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
//            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
//        }
//    }
//    stbi_write_png("scenegammacorrected.png", W, H, 3, &image[0], 0);
//
//    return 0;
//}


////////////////////////////////////////////////////////////////////////////

// Scene with gamma correction and shadows

////////////////////////////////////////////////////////////////////////////


#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.1415

class Vector3 { // Defining a class for vectors
public:
    explicit Vector3(double x = 0, double y = 0, double z = 0) {
        coordinates[0] = x;
        coordinates[1] = y;
        coordinates[2] = z;
    };
    double operator[](int i) const { return coordinates[i]; };
    double& operator[](int i) { return coordinates[i]; };
    double NormSquared() { // squared norm
        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
    }
    Vector3 Normalize() { // to normalize a vector
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

Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a double
    return Vector3(a * b[0], a * b[1], a * b[2]);
}

Vector3 operator*(const Vector3& a, double b) { //idem
    return Vector3(a[0] * b, a[1] * b, a[2] * b);
}

Vector3 operator/(const Vector3& a, double b) { //division of a vector by a double
    return Vector3(a[0] / b, a[1] / b, a[2] / b);
}


double Dot(const Vector3& a, const Vector3& b) { //dot product between two vectors
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}


class Ray { //Class for rays (characterized by their center and their direction)
public:
    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
    }
    Vector3 C, u;
};

class Sphere { //class for spheres
public:
    Sphere(const Vector3& O, double R, const Vector3& albedo) : O(O), R(R), albedo(albedo) {
    }
    bool Intersection(const Ray& r, Vector3& P, Vector3& N, double& t) {
        // solves a*t² + b*t + c = 0
        double a = 1;
        double b = 2 * Dot(r.u, r.C - O);
        double c = (r.C - O).NormSquared() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; //no solution, thus no intersection

        double sqDelta = sqrt(delta);
        double t2 = (-b + sqDelta) / (2 * a);
        if (t2 < 0) return false; //no positive solution, thus no intersection

        double t1 = (-b - sqDelta) / (2 * a);
        if (t1 > 0) // t is the smallest positive value between t1 and t2
            t = t1;

        else
            t = t2;

        P = r.C + t * r.u; // Point of intersection
        N = (P - O).Normalize(); //Normal at the point of intersection

        return true;

    };
    Vector3 O;
    double R;
    Vector3 albedo;
};

class Scene { // Class for the scene, which can contain several spheres
public:
    Scene() {};
    bool Intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t) {
        t = 1E9;
        bool intersected = false;
        for (int i = 0; i < objects.size(); i++) { // for each sphere in the scene
            Vector3 Pobj, Nobj;
            double tobj;
            if (objects[i].Intersection(r, Pobj, Nobj, tobj) && tobj < t) {  // if there is an intersection closer than existing ones, then we take this one into account and not the others
                intersected = true;
                t = tobj;
                P = Pobj;
                N = Nobj;
                albedo = objects[i].albedo;
            }

        }
        return intersected;
    }
    std::vector<Sphere> objects;

};

int main() {
    int W = 512;
    int H = 512;

    Vector3 C(0, 0, 55);
    Scene scene;


    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0)); // bright red
    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(0.8, 0.2, 0.9)); // magenta
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



    double fov = 70 * M_PI / 180;
    double I = 10000000000; //light intensity
    Vector3 L(0, 20, 40); //light source

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
            u = u.Normalize();

            Ray r(C, u);
            double t;
            Vector3 P, N, albedo;
            Vector3 color(0, 0, 0);
            if (scene.Intersection(r, P, N, albedo, t)) {// in case of intersection with an object of the scene
                double normPL = sqrt((L - P).NormSquared());//norm of PL
                Vector3 Pshadow, Nshadow, albedoshadow;
                double tshadow;
                Ray Rayshadow(P + 0.001 * N, (L - P) / normPL); // Ray starting from the point of intersection and directed towards the light source
                if (scene.Intersection(Rayshadow, Pshadow, Nshadow, albedoshadow, tshadow) && tshadow < normPL) { //if there is an intersection before reaching the light source
                    color = Vector3(0., 0., 0.); //not lit = shadow
                }
                else {
                    double prov = std::max(0., Dot(N, (L - P) / normPL));
                    color = I / (4 * M_PI * normPL * normPL) * albedo / M_PI * prov;
                }

            }

            color[0] = std::pow(color[0], 0.45);//applying gamma
            color[1] = std::pow(color[1], 0.45);
            color[2] = std::pow(color[2], 0.45);

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("scenewshadows.png", W, H, 3, &image[0], 0);

    return 0;
}

