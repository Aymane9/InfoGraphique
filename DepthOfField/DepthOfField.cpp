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

class Vector3 {
public:
    explicit Vector3(double x = 0, double y = 0, double z = 0) {
        coordonnees[0] = x;
        coordonnees[1] = y;
        coordonnees[2] = z;
    };
    double operator[](int i) const { return coordonnees[i]; };
    double& operator[](int i) { return coordonnees[i]; };
    double NormSquared() {
        return coordonnees[0] * coordonnees[0] + coordonnees[1] * coordonnees[1] + coordonnees[2] * coordonnees[2];
    }
    Vector3 Normalize() {
        double norme = sqrt(NormSquared());
        return Vector3(coordonnees[0] / norme, coordonnees[1] / norme, coordonnees[2] / norme);
    }
private:
    double coordonnees[3];
};

Vector3 operator+(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}

Vector3 operator-(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}

Vector3 operator-(const Vector3& a) {
    return Vector3(-a[0], -a[1], -a[2]);
}

Vector3 operator*(double a, const Vector3& b) {
    return Vector3(a * b[0], a * b[1], a * b[2]);
}


Vector3 operator*(const Vector3& a, double b) {
    return Vector3(a[0] * b, a[1] * b, a[2] * b);
}

Vector3 operator*(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector3 operator/(const Vector3& a, double b) {
    return Vector3(a[0] / b, a[1] / b, a[2] / b);
}

Vector3 CrossProduct(const Vector3& a, const Vector3& b) {
    return Vector3(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);

}

double DotProduct(const Vector3& a, const Vector3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector3 TermByTermProduct(const Vector3& a, const Vector3& b) {
    return Vector3(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}



Vector3 RandomInUnitSphere(const Vector3& N) {
    double u1 = uniform(engine);
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



class Ray {
public:
    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
    }
    Vector3 C, u;
};

class Sphere {
public:
    Sphere(const Vector3& O, double R, const Vector3& albedo, bool mirror = false, bool Transparent = false) : O(O), R(R), albedo(albedo), mirror(mirror), Transparent(Transparent) {
    }
    bool intersection(const Ray& r, Vector3& P, Vector3& N, double& t) {
        double a = 1;
        double b = 2 * DotProduct(r.u, r.C - O);
        double c = (r.C - O).NormSquared() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false;

        double sqrtDelta = sqrt(delta);
        double t2 = (-b + sqrtDelta) / (2 * a);
        if (t2 < 0) return false;

        double t1 = (-b - sqrtDelta) / (2 * a);
        if (t1 > 0)
            t = t1;

        else
            t = t2;

        P = r.C + t * r.u;
        N = (P - O).Normalize();

        return true;

    };
    Vector3 O;
    double R;
    Vector3 albedo;
    bool mirror, Transparent;
};

class Scene {
public:
    Scene() {};
    bool intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& mirror, bool& transp, int& id) {
        t = 1E9;
        bool intersecte = false;
        for (int i = 0; i < objects.size(); i++) {
            Vector3 Pobjet, Nobjet;
            double tobjet;
            if (objects[i].intersection(r, Pobjet, Nobjet, tobjet) && tobjet < t) {
                intersecte = true;
                t = tobjet;
                P = Pobjet;
                N = Nobjet;
                albedo = objects[i].albedo;
                mirror = objects[i].mirror;
                transp = objects[i].Transparent;
                id = i;
            }

        }
        return intersecte;
    }

    Vector3 obtenircolor(const Ray& r, int bounce, bool lastdiffusion) {
        if (bounce > 5) return Vector3(0., 0., 0.);
        else {
            double t;
            bool mirror, transp;
            Vector3 P, N, albedo;
            Vector3 color(0, 0, 0);
            int id;
            if (intersection(r, P, N, albedo, t, mirror, transp, id)) {
                if (id == 0) {
                    if (bounce == 0 || !lastdiffusion)
                        return lightIntensity / (4 * M_PI * M_PI * objects[0].R * objects[0].R);
                    else
                        return Vector3(0, 0, 0);
                }

                if (mirror) {
                    Vector3 reflectedDirection = r.u - 2 * DotProduct(r.u, N) * N;
                    Ray reflectedRay(P + 0.00001 * N, reflectedDirection);
                    return obtenircolor(reflectedRay, bounce + 1, false);
                }
                else {
                    if (transp) {

                        double n1 = 1., n2 = 1.4;
                        Vector3 N2 = N;
                        if (DotProduct(r.u, N) > 0) {
                            std::swap(n1, n2);
                            N2 = -N;
                        }
                        double angle = 1 - n1 * n1 / (n2 * n2) * (1 - DotProduct(r.u, N2) * DotProduct(r.u, N2));
                        if (angle < 0) {
                            Vector3 reflectedDirection = r.u - 2 * DotProduct(r.u, N) * N;
                            Ray reflectedRay(P + 0.00001 * N, reflectedDirection);
                            return obtenircolor(reflectedRay, bounce + 1, false);
                        }
                        Vector3 Tt = n1 / n2 * (r.u - DotProduct(r.u, N2) * N2);
                        Vector3 Tn = -sqrt(angle) * N2;
                        Vector3 refractedDirection = Tt + Tn;
                        Ray refractedRay(P - 0.0001 * N2, refractedDirection);
                        return obtenircolor(refractedRay, bounce + 1, false);
                    }
                    else {

                        Vector3 w = RandomInUnitSphere((P - lightPosition).Normalize());
                        Vector3 xp = w * objects[0].R + objects[0].O;
                        Vector3 Pxp = xp - P;
                        double normePxp = sqrt(Pxp.NormSquared());
                        Pxp = Pxp.Normalize();
                        Vector3 Pshadow, Nshadow, albedoshadow;
                        double tshadow;
                        bool mirrorshadow, transshadow;
                        int idshadow;
                        Ray Rayshadow(P + 0.00001 * N, Pxp);
                        if (intersection(Rayshadow, Pshadow, Nshadow, albedoshadow, tshadow, mirrorshadow, transshadow, idshadow) && tshadow < normePxp - 0.0001) {
                            color = Vector3(0., 0., 0.);
                        }
                        else {
                            color = lightIntensity / (4 * M_PI * M_PI * objects[0].R * objects[0].R) * albedo / M_PI * std::max(0., DotProduct(N, Pxp)) * std::max(0., -DotProduct(w, Pxp)) / (normePxp * normePxp) / (std::max(0., -DotProduct((lightPosition - P).Normalize(), w)) / (M_PI * objects[0].R * objects[0].R));
                        }




                        Vector3 random = RandomInUnitSphere(N);
                        Ray Rayrandom(P + 0.00001 * N, random);
                        color = color + TermByTermProduct(albedo, obtenircolor(Rayrandom, bounce + 1, true));


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

    Vector3 C(0, 0, 55);
    Scene scene;
    scene.lightPosition = Vector3(-10, 20, 40);

    Sphere SLight(scene.lightPosition, 5, Vector3(1., 0.3, 0.2));
    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0));
    Sphere S2(Vector3(-10, 0, -20), 10, Vector3(0., 1., 0.));
    Sphere S3(Vector3(10, 0, 20), 10, Vector3(0.7, 0., 0.7));
    Sphere Sfront(Vector3(0, 0, 1000), 940, Vector3(1.0, 0.0, 0.0)); // red
    Sphere Sback(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.)); // blue
    Sphere Stop(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7)); // violet
    Sphere Sbottom(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2)); // yellow
    Sphere Sright(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.)); // green
    Sphere Sleft(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1)); // brown

    scene.objects.push_back(SLight);
    scene.objects.push_back(S1);
    scene.objects.push_back(S2);
    scene.objects.push_back(S3);
    scene.objects.push_back(Sleft);
    scene.objects.push_back(Sright);
    scene.objects.push_back(Stop);
    scene.objects.push_back(Sbottom);
    scene.objects.push_back(Sback);
    scene.objects.push_back(Sfront);


    double fov = 60 * M_PI / 180;
    scene.lightIntensity = Vector3(5E9, 5E9, 5E9);
    int nmbRay = 60;

    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {


            Vector3 color(0, 0, 0);
            for (int k = 0; k < nmbRay; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
                Vector3 u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2))); // antialiasing
                u = u.Normalize();

                double u3 = uniform(engine);
                double u4 = uniform(engine);
                double x3 = 1 * cos(2 * M_PI * u3) * sqrt(-2 * log(u4));
                double x4 = 1 * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

                

                Vector3 Focalpoint = C + 55 * u; 
                Vector3 RandOrigin = C + Vector3(x3, x4, 0); 
                Vector3 Raydir = (Focalpoint - RandOrigin).Normalize(); 

                Ray r(RandOrigin, Raydir);

                color = color + scene.obtenircolor(r, 0, false);

            }
            color = color / nmbRay;


            color[0] = std::pow(color[0], 0.45);// gamma correction
            color[1] = std::pow(color[1], 0.45);
            color[2] = std::pow(color[2], 0.45);


            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("softshadows.png", W, H, 3, &image[0], 0);

    return 0;
}

