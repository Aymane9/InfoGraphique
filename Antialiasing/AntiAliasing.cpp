////////////////////////////////////////////////////////////////////////////

// Adding indirect lighting

////////////////////////////////////////////////////////////////////////////

//#define _CRT_SECURE_NO_WARNINGS 1
//#include <vector>
//#include <algorithm>
//#include <cmath>
//#include <iostream>
//
//#define STB_IMAGE_WRITE_IMPLEMENTATION
//#include "stb_image_write.h"
//
//#define STB_IMAGE_IMPLEMENTATION
//#include "stb_image.h"
//#include <random>
//
//#define M_PI 3.1415
//
//static std::default_random_engine engine(10);
//static std::uniform_real_distribution<double> uniform(0, 1);
//
//class Vector3 { //Defining a class for 3D vectors
//public:
//    explicit Vector3(double x = 0, double y = 0, double z = 0) {
//        coordinates[0] = x;
//        coordinates[1] = y;
//        coordinates[2] = z;
//    };
//    double operator[](int i) const { return coordinates[i]; };
//    double& operator[](int i) { return coordinates[i]; };
//    double NormSquared() { //squared norm
//        return coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] + coordinates[2] * coordinates[2];
//    }
//    Vector3 Normalize() { //normalize a vector
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
//Vector3 operator-(const Vector3& a) { //inverse of a vector
//    return Vector3(-a[0], -a[1], -a[2]);
//}
//
//Vector3 operator*(double a, const Vector3& b) { //multiplication of a vector by a scalar
//    return Vector3(a * b[0], a * b[1], a * b[2]);
//}
//
//Vector3 operator*(const Vector3& a, double b) { //same as above
//    return Vector3(a[0] * b, a[1] * b, a[2] * b);
//}
//
//Vector3 operator/(const Vector3& a, double b) { //division of a vector by a scalar
//    return Vector3(a[0] / b, a[1] / b, a[2] / b);
//}
//
//Vector3 operator*(const Vector3& a, const Vector3& b) {
//    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
//}
//
//Vector3 CrossProduct(const Vector3& a, const Vector3& b) { //cross product of two vectors
//    return Vector3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
//}
//
//double DotProduct(const Vector3& a, const Vector3& b) { //dot product of two vectors
//    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
//}
//
//Vector3 TermByTermProduct(const Vector3& a, const Vector3& b) { //term-by-term product of two vectors
//    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
//}
//
//Vector3 RandomInUnitSphere(const Vector3& N) { //Returns a vector zN+xT1+y*T2, where (N,T1,T2) is a coordinate system, and x,y,z are random variables follorandomng a cosine probability distribution (the radius is more likely to be close to N)
//    double u1 = uniform(engine); // random number between 0 and 1
//    double u2 = uniform(engine);
//    double x = cos(2 * M_PI * u1) * sqrt(1 - u2);
//    double y = sin(2 * M_PI * u1) * sqrt(1 - u2);
//    double z = sqrt(u2);
//    Vector3 T1;
//    if (N[0] < N[1] && N[0] < N[2]) {
//        T1 = Vector3(0, N[2], -N[1]);
//    }
//    else {
//        if (N[1] < N[2] && N[1] < N[0]) {
//            T1 = Vector3(N[2], 0, -N[0]);
//        }
//        else {
//            T1 = Vector3(N[1], -N[0], 0);
//        }
//    }
//        T1 = T1.Normalize();
//        Vector3 T2 = CrossProduct(N, T1);
//        return x * T1 + y * T2 + z * N;
//}
//
//class Ray { //Class for rays (characterized by their origin and direction)
//public:
//    Ray(const Vector3& origin, const Vector3& direction) : origin(origin), direction(direction) {
//    }
//    Vector3 origin, direction;
//};
//
//class Sphere { //Class for spheres (added an attribute to indicate if it's a mirror or not, and another for transparency)
//public:
//    Sphere(const Vector3& center, double radius, const Vector3& albedo, bool isMirror = false, bool isTransparent = false) : center(center), radius(radius), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent) {
//    }
//    bool Intersect(const Ray& ray, Vector3& P, Vector3& N, double& t) {
//        // solve a*t² + b*t + c = 0
//        double a = 1;
//        double b = 2 * DotProduct(ray.direction, ray.origin - center);
//        double c = (ray.origin - center).NormSquared() - radius * radius;
//        double delta = b * b - 4 * a * c;
//        if (delta < 0) return false; //no solution, no intersection
//
//        double sqrtDelta = sqrt(delta);
//        double t2 = (-b + sqrtDelta) / (2 * a);
//        if (t2 < 0) return false; //no positive solution, no intersection
//
//        double t1 = (-b - sqrtDelta) / (2 * a);
//        if (t1 > 0) // t is the smallest positive value between t1 and t2
//            t = t1;
//        else
//            t = t2;
//
//        P = ray.origin + t * ray.direction; //Intersection point
//        N = (P - center).Normalize(); //Normal at the intersection point
//
//        return true;
//    }
//    Vector3 center;
//    double radius;
//    Vector3 albedo;
//    bool isMirror, isTransparent;
//};
//
//
//
//class Scene { // classe de la scene, qui peut comporter plusieurs sphères
//public:
//    Scene() {};
//    bool intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& miroir, bool& transp) {
//        t = 1E9;
//        bool intersecte = false;
//        for (int i = 0; i < objects.size(); i++) {// pour chacune des sphères de la scène
//            Vector3 Pobjet, Nobjet;
//            double tobjet;
//            if (objects[i].Intersect(r, Pobjet, Nobjet, tobjet) && tobjet < t) { //s'il y a une intersection plus proche que celles existantes, alors on prend en compte celle-ci et pas les autres
//                intersecte = true;
//                t = tobjet;
//                P = Pobjet;
//                N = Nobjet;
//                albedo = objects[i].albedo;
//                miroir = objects[i].isMirror;
//                transp = objects[i].isTransparent;
//
//            }
//
//        }
//        return intersecte;
//    }
//
//    Vector3 getColor(const Ray& r, int rebond) { //pour obtenir la couleur
//        if (rebond > 5) return Vector3(0., 0., 0.); // si on dépasse 5 rebonds, on renvoit la couleur noire
//        else {
//            double t;
//            bool miroir, transp;
//            Vector3 P, N, albedo;
//            Vector3 color(0, 0, 0);
//            if (intersection(r, P, N, albedo, t, miroir, transp)) {// en cas d'intersection avec un objet de la scène
//                if (miroir) { // si c'est un miroir
//                    Vector3 reflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N; // direction de la réflexion
//                    Ray reflectedRay(P + 0.00001 * N, reflectedDirection); //rayon réfléchi, partant du point d'intersection et allant dans la direction réfléchie
//                    return getColor(reflectedRay, rebond + 1);
//                }
//                else {
//                    if (transp) { // si  il est transparent
//
//                        double n1 = 1., n2 = 1.4; //indices optiques
//                        Vector3 N2 = N;
//                        if (DotProduct(r.direction, N) > 0) { //si on sort de la sphère, on inverse les n et la normale
//                            std::swap(n1, n2);
//                            N2 = -N;
//                        }
//                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - DotProduct(r.direction, N2) * DotProduct(r.direction, N2));
//                        if (angle < 0) { // si l'angle est plus petit que 0, il y a reflexion
//                            Vector3 reflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N;// direction de la réflexion
//                            Ray reflectedRay(P + 0.00001 * N, reflectedDirection);//rayon réfléchi, partant du point d'intersection et allant dans la direction réfléchie
//                            return getColor(reflectedRay, rebond + 1);
//                        }
//                        Vector3 Tt = n1 / n2 * (r.direction - DotProduct(r.direction, N2) * N2); //Composante tangentielle de la direction de réfraction
//                        Vector3 Tn = -sqrt(angle) * N2; //Composante normale de la direction de réfraction
//                        Vector3 refractedDirection = Tt + Tn; //direction de réfraction
//                        Ray refractedRay(P - 0.0001 * N2, refractedDirection); //rayon réfracté, partant du point d'intersection et allant dans la direction réfractée
//                        return getColor(refractedRay, rebond + 1);
//                    }
//                    else {
//                        //direct lighting
//                        double normePL = sqrt((lightPosition - P).NormSquared());//norme de PL
//                        Vector3 Pshadow, Nshadow, albedoshadow;
//                        double tshadow;
//                        bool miroirshadow, transshadow;
//                        Ray Rayonshadow(P + 0.001 * N, (lightPosition - P) / normePL); // Rayon partant du point d'intersection et dirigé vers la source de lumière
//                        if (intersection(Rayonshadow, Pshadow, Nshadow, albedoshadow, tshadow, miroirshadow, transshadow) && tshadow < normePL) {//s'il y a intersection avant d'arriver à la source de lumière
//                            color = Vector3(0., 0., 0.); //pas éclairé = shadow
//                        }
//                        else {
//                            double prov = std::max(0., DotProduct(N, (lightPosition - P) / normePL));
//                            color = lightIntensity / (4 * M_PI * normePL * normePL) * albedo / M_PI * prov;
//                        }
//                        //indirect lighting
//                        Vector3 random = RandomInUnitSphere(N);
//                        Ray Rayonrandom(P + 0.00001 * N, random); //rayon aléatoire
//                        color = color + TermByTermProduct(albedo, getColor(Rayonrandom, rebond + 1));
//
//
//                    }
//                }
//            }
//            return color;
//        }
//    }
//    std::vector<Sphere> objects;
//    Vector3 lightPosition;
//    Vector3 lightIntensity;
//
//};
//
//
//
//int main() {
//    int W = 512;
//    int H = 512;
//
//    Vector3 cameraPosition(0, 0, 55);
//    Scene scene;
//
//    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0), false, true);
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
//
//    double fov = 70 * M_PI / 180;
//    scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
//    scene.lightPosition = Vector3(-10, 20, 40);
//    int nmbRays = 3; 
//
//    std::vector<unsigned char> image(W * H * 3, 0);
//    for (int i = 0; i < H; i++) {
//        for (int j = 0; j < W; j++) {
//
//
//            Vector3 u(j - W / 2, i - H / 2, -W / (2. * tan(fov / 2)));
//            u = u.Normalize();
//
//            Ray r(cameraPosition, u);
//
//            Vector3 color(0, 0, 0);
//            for (int k = 0; k < nmbRays; k++)
//                color = color + scene.getColor(r, 0);
//            color = color / nmbRays;
//
//
//            color[0] = std::pow(color[0], 0.45);//application du gamma
//            color[1] = std::pow(color[1], 0.45);
//            color[2] = std::pow(color[2], 0.45);
//
//
//            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
//            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
//            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
//        }
//    }
//    stbi_write_png("indirectlight.png", W, H, 3, &image[0], 0);
//
//    return 0;
//}
//
//
//
//

////////////////////////////////////////////////////////////////////////////

// Adding AntiAliasing

////////////////////////////////////////////////////////////////////////////

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

Vector3 RandomInUnitSphere(const Vector3& N) { //Returns a vector zN+xT1+y*T2, where (N,T1,T2) is a coordinate system, and x,y,z are random variables follorandomng a cosine probability distribution (the radius is more likely to be close to N)
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

class Scene { // Scene class, which can contain multiple spheres
public:
    Scene() {};
    bool intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& mirror, bool& transparent) {
        t = 1E9;
        bool intersect = false;
        for (int i = 0; i < objects.size(); i++) {// for each sphere in the scene
            Vector3 Pobj, Nobj;
            double tobj;
            if (objects[i].Intersect(r, Pobj, Nobj, tobj) && tobj < t) { // if there is a closer intersection than the existing ones, then we take this one into account and not the others
                intersect = true;
                t = tobj;
                P = Pobj;
                N = Nobj;
                albedo = objects[i].albedo;
                mirror = objects[i].isMirror;
                transparent = objects[i].isTransparent;

            }

        }
        return intersect;
    }

    Vector3 getColor(const Ray& r, int bounce) { // to get the color
        if (bounce > 5) return Vector3(0., 0., 0.); // if we exceed 5 bounces, we return the black color
        else {
            double t;
            bool mirror, transparent;
            Vector3 P, N, albedo;
            Vector3 color(0, 0, 0);
            if (intersection(r, P, N, albedo, t, mirror, transparent)) {// in case of intersection randomth an object in the scene
                if (mirror) { // if it's a mirror
                    Vector3 ReflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N; // reflection direction
                    Ray ReflectedRay(P + 0.00001 * N, ReflectedDirection); //reflected ray, starting from the intersection point and going in the reflected direction
                    return getColor(ReflectedRay, bounce + 1);
                }
                else {
                    if (transparent) { // if it is transparent

                        double n1 = 1., n2 = 1.4; //optical indices
                        Vector3 N2 = N;
                        if (DotProduct(r.direction, N) > 0) { //if we exit the sphere, we reverse the n and the normal
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - DotProduct(r.direction, N2) * DotProduct(r.direction, N2));
                        if (angle < 0) { // if the angle is smaller than 0, there is reflection
                            Vector3 ReflectedDirection = r.direction - 2 * DotProduct(r.direction, N) * N;// reflection direction
                            Ray ReflectedRay(P + 0.00001 * N, ReflectedDirection);//reflected ray, starting from the intersection point and going in the reflected direction
                            return getColor(ReflectedRay, bounce + 1);
                        }
                        Vector3 Tt = n1 / n2 * (r.direction - DotProduct(r.direction, N2) * N2); //Tangential component of the refracted direction
                        Vector3 Tn = -sqrt(angle) * N2; //Normal component of the refracted direction
                        Vector3 RefractedDirection = Tt + Tn; //refracted direction
                        Ray RefractedRay(P - 0.0001 * N2, RefractedDirection); //refracted ray, starting from the intersection point and going in the refracted direction
                        return getColor(RefractedRay, bounce + 1);
                    }
                    else {
                        //direct lighting
                        double normPL = sqrt((lightPosition - P).NormSquared());//PL norm
                        Vector3 Pshadow, Nshadow, albedoshadow;
                        double tshadow;
                        bool mirrorshadow, transparentshadow;
                        Ray ShadowRay(P + 0.001 * N, (lightPosition - P) / normPL); // Ray starting from the intersection point and directed towards the light source
                        if (intersection(ShadowRay, Pshadow, Nshadow, albedoshadow, tshadow, mirrorshadow, transparentshadow) && tshadow < normPL) {// if there is an intersection before reaching the light source
                            color = Vector3(0., 0., 0.); //not illuminated = shadow
                        }
                        else {
                            double prov = std::max(0., DotProduct(N, (lightPosition - P) / normPL));
                            color = lightIntensity / (4 * M_PI * normPL * normPL) * albedo / M_PI * prov;
                        }
                        //indirect lighting
                        Vector3 random = RandomInUnitSphere(N);
                        Ray RandomRay(P + 0.00001 * N, random); //random ray
                        color = color + TermByTermProduct(albedo, getColor(RandomRay, bounce + 1));


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


    double fov = 70 * M_PI / 180;
    scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
    scene.lightPosition = Vector3(-10, 20, 40);
    int nmbRays = 100;

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            double u1 = uniform(engine);
            double u2 = uniform(engine);
            double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
            double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
            Vector3 u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2)));
            u = u.Normalize();

            Ray r(cameraPosition, u);

            Vector3 color(0, 0, 0);
            for (int k = 0; k < nmbRays; k++)
                color = color + scene.getColor(r, 0);
            color = color / nmbRays;


            color[0] = std::pow(color[0], 0.45);//application du gamma
            color[1] = std::pow(color[1], 0.45);
            color[2] = std::pow(color[2], 0.45);


            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("Antialiasing.png", W, H, 3, &image[0], 0);

    return 0;
}




