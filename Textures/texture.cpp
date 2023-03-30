////////////////////////////////////////////////////////////////////////////

// Adding textures to the model

////////////////////////////////////////////////////////////////////////////

#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdio.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <random>
#include <list>

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
        double norm = sqrt(NormSquared());
        return Vector3(coordonnees[0] / norm, coordonnees[1] / norm, coordonnees[2] / norm);
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
    return z * N + x * T1 + y * T2;

}



class Ray { 
public:
    Ray(const Vector3& C, const Vector3& u) : C(C), u(u) {
    }
    Vector3 C, u;
};

class Object { 
public:
    Object() {};

    virtual bool intersection(const Ray& r, Vector3& P, Vector3& Normale, double& t, Vector3& color) = 0;
    Vector3 albedo;
    bool Mirror, Transparent;
};

class Sphere : public Object { 
public:
    Sphere(const Vector3& O, double R, const Vector3& albedo, bool Mirror = false, bool Transparent = false) : O(O), R(R) {
        this->albedo = albedo;
        this->Mirror = Mirror;
        this->Transparent = Transparent;
    };
    bool intersection(const Ray& r, Vector3& P, Vector3& N, double& t, Vector3& color) {
        // solve for a*t² + b*t + c = 0
        double a = 1;
        double b = 2 * DotProduct(r.u, r.C - O);
        double c = (r.C - O).NormSquared() - R * R;
        double delta = b * b - 4 * a * c;
        if (delta < 0) return false; 

        double racinedelta = sqrt(delta);
        double t2 = (-b + racinedelta) / (2 * a);
        if (t2 < 0) return false; 

        double t1 = (-b - racinedelta) / (2 * a);
        if (t1 > 0) 
            t = t1;

        else
            t = t2;

        P = r.C + t * r.u; 
        N = (P - O).Normalize(); 
        color = this->albedo;





        return true;

    };
    Vector3 O;
    double R;

};


class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
    };
    int vtxi, vtxj, vtxk; 
    int uvi, uvj, uvk;  
    int ni, nj, nk;  
    int group;       
};

class BoundingBox {
public:
    bool intersection(const Ray& r) {
        double x1 = (mini[0] - r.C[0]) / r.u[0];
        double x2 = (maxi[0] - r.C[0]) / r.u[0];
        double xmin = std::min(x1, x2);
        double xmax = std::max(x1, x2);

        double y1 = (mini[1] - r.C[1]) / r.u[1];
        double y2 = (maxi[1] - r.C[1]) / r.u[1];
        double ymin = std::min(y1, y2);
        double ymax = std::max(y1, y2);

        double z1 = (mini[2] - r.C[2]) / r.u[2];
        double z2 = (maxi[2] - r.C[2]) / r.u[2];
        double zmin = std::min(z1, z2);
        double zmax = std::max(z1, z2);

        double max = std::min(xmax, std::min(ymax, zmax));
        double min = std::max(xmin, std::max(ymin, zmin));
        if (max < 0) return false;
        return max > min;
    }
    Vector3 mini, maxi;
};

class Node {
public:
    Node* fg, * fd;
    BoundingBox b;
    int debut, fin;

};


class TriangleMesh : public Object {
public:
    ~TriangleMesh() {}
    TriangleMesh(const Vector3& albedo, bool Mirror = false, bool Transparent = false) {
        this->albedo = albedo;
        this->Mirror = Mirror;
        this->Transparent = Transparent;
        node = new Node;
    };

    BoundingBox boundingBox(int debut, int fin) {
        bb.mini = Vector3(1E9, 1E9, 1E9);
        bb.maxi = Vector3(-1E9, -1E9, -1E9);
        for (int i = debut; i < fin; i++) {
            for (int j = 0; j < 3; j++) {
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxi][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxi][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxj][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxj][j]);
                bb.mini[j] = std::min(bb.mini[j], vertices[indices[i].vtxk][j]);
                bb.maxi[j] = std::max(bb.maxi[j], vertices[indices[i].vtxk][j]);
            }
        }
        return bb;
    }

    void Pivot(Node* n, int debut, int fin) {

        n->debut = debut;
        n->fin = fin;
        n->b = boundingBox(n->debut, n->fin); 
        Vector3 lenghts = n->b.maxi - n->b.mini;
        int dim;
        if (lenghts[0] >= lenghts[1] && lenghts[0] >= lenghts[2]) { 
            dim = 0;
        }
        else {
            if (lenghts[1] >= lenghts[0] && lenghts[1] >= lenghts[2]) {
                dim = 1;
            }
            else {
                dim = 2;
            }
        }
        double center = (n->b.mini[dim] + n->b.maxi[dim]) / 2; 
        int pivot = n->debut; 
        for (int i = n->debut; i < n->fin; i++) { 
            double barycentre = (vertices[indices[i].vtxi][dim] + vertices[indices[i].vtxj][dim] + vertices[indices[i].vtxk][dim]) / 3; 
            if (barycentre < center) {  
                std::swap(indices[i], indices[pivot]);
                pivot++; 
            }
        }


        if (pivot == debut || pivot == fin || (fin - debut < 5)) return; 

        n->fg = new Node;
        n->fd = new Node;

        Pivot(n->fg, n->debut, pivot); 
        Pivot(n->fd, pivot, n->fin);

    }

    void readOBJ(const char* obj) {

        char matfile[255];
        char grp[255];

        FILE* f;
        f = fopen(obj, "r");
        int curGroup = -1;
        while (!feof(f)) {
            char line[255];
            if (!fgets(line, 255, f)) break;

            std::string linetrim(line);
            linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
            strcpy(line, linetrim.c_str());

            if (line[0] == 'u' && line[1] == 's') {
                sscanf(line, "usemtl %[^\n]\n", grp);
                curGroup++;
            }

            if (line[0] == 'v' && line[1] == ' ') {
                Vector3 vec;

                Vector3 col;
                if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
                    col[0] = std::min(1., std::max(0., col[0]));
                    col[1] = std::min(1., std::max(0., col[1]));
                    col[2] = std::min(1., std::max(0., col[2]));

                    vertices.push_back(vec);
                    vertexcolors.push_back(col);

                }
                else {
                    sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                    vertices.push_back(vec);
                }
            }
            if (line[0] == 'v' && line[1] == 'n') {
                Vector3 vec;
                sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
                normals.push_back(vec);
            }
            if (line[0] == 'v' && line[1] == 't') {
                Vector3 vec;
                sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
                uvs.push_back(vec);
            }
            if (line[0] == 'f') {
                TriangleIndices t;
                int i0, i1, i2, i3;
                int j0, j1, j2, j3;
                int k0, k1, k2, k3;
                int nn;
                t.group = curGroup;

                char* consumedline = line + 1;
                int offset;

                nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
                if (nn == 9) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                    if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                    if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                    if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                    if (nn == 6) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                        if (j0 < 0) t.uvi = uvs.size() + j0; else   t.uvi = j0 - 1;
                        if (j1 < 0) t.uvj = uvs.size() + j1; else   t.uvj = j1 - 1;
                        if (j2 < 0) t.uvk = uvs.size() + j2; else   t.uvk = j2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                        if (nn == 3) {
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            indices.push_back(t);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                            if (i0 < 0) t.vtxi = vertices.size() + i0; else t.vtxi = i0 - 1;
                            if (i1 < 0) t.vtxj = vertices.size() + i1; else t.vtxj = i1 - 1;
                            if (i2 < 0) t.vtxk = vertices.size() + i2; else t.vtxk = i2 - 1;
                            if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                            if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                            if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                            indices.push_back(t);
                        }
                    }
                }

                consumedline = consumedline + offset;

                while (true) {
                    if (consumedline[0] == '\n') break;
                    if (consumedline[0] == '\0') break;
                    nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                    TriangleIndices t2;
                    t2.group = curGroup;
                    if (nn == 3) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                        if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                        if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                        if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                        indices.push_back(t2);
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        k2 = k3;
                    }
                    else {
                        nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (j0 < 0) t2.uvi = uvs.size() + j0; else  t2.uvi = j0 - 1;
                            if (j2 < 0) t2.uvj = uvs.size() + j2; else  t2.uvj = j2 - 1;
                            if (j3 < 0) t2.uvk = uvs.size() + j3; else  t2.uvk = j3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            j2 = j3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                            if (nn == 2) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                if (k0 < 0) t2.ni = normals.size() + k0; else   t2.ni = k0 - 1;
                                if (k2 < 0) t2.nj = normals.size() + k2; else   t2.nj = k2 - 1;
                                if (k3 < 0) t2.nk = normals.size() + k3; else   t2.nk = k3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                k2 = k3;
                                indices.push_back(t2);
                            }
                            else {
                                nn = sscanf(consumedline, "%u%n", &i3, &offset);
                                if (nn == 1) {
                                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                    consumedline = consumedline + offset;
                                    i2 = i3;
                                    indices.push_back(t2);
                                }
                                else {
                                    consumedline = consumedline + 1;
                                }
                            }
                        }
                    }
                }

            }

        }
        fclose(f);

    }
    bool intersection(const Ray& r, Vector3& P, Vector3& Normale, double& t, Vector3& color) {
        if (!node->b.intersection(r)) return false; 
        t = 1E9;
        bool intersect = false;
        std::list<Node*> l;
        l.push_back(node);
        while (!l.empty()) {
            Node* c = l.front(); 
            l.pop_front(); 
            if (c->fg) {
                if (c->fg->b.intersection(r)) { 
                }
                if (c->fd->b.intersection(r)) {
                    l.push_front(c->fd);
                }
            }
            else { 
                for (int i = c->debut; i < c->fin; i++) {
                    const Vector3& A = vertices[indices[i].vtxi]; 
                    const Vector3& B = vertices[indices[i].vtxj];
                    const Vector3& C = vertices[indices[i].vtxk];
                    Vector3 e1 = B - A; 
                    Vector3 e2 = C - A; 
                    Vector3 N = CrossProduct(e1, e2); 
                    Vector3 AC = r.C - A;
                    Vector3 ACvectU = CrossProduct(AC, r.u);
                    double UscalN = DotProduct(r.u, N);
                    double beta = -DotProduct(e2, ACvectU) / UscalN;
                    double gamma = DotProduct(e1, ACvectU) / UscalN;
                    double alpha = 1 - beta - gamma;
                    double tObject = -DotProduct(AC, N) / UscalN;
                    if (beta >= 0 && gamma >= 0 && alpha >= 0 && beta <= 1 && gamma <= 1 && alpha <= 1 && tObject > 0) { //Dans le cas d'une intersection
                        intersect = true;
                        if (tObject < t) { 
                            t = tObject;
                            Normale = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk]; //on prend la normale définie dans le fichier obj
                            Normale = Normale.Normalize();
                            P = r.C + t * r.u;
                            Vector3 UV = alpha * uvs[indices[i].uvi] + beta * uvs[indices[i].uvj] + gamma * uvs[indices[i].uvk]; // coordonnees sur l'image du point d'intersection
                            int W = Wtex[0];
                            int H = Htex[0];
                            UV = TermByTermProduct(UV, Vector3(W, H, 0)); 
                            int uvx = UV[0] + 0.5; 
                            int uvy = UV[1] + 0.5;
                            uvx = uvx % W; 
                            uvy = uvy % H;
                            if (uvx < 0) uvx += W; 
                            if (uvy < 0) uvy += H;
                            uvy = H - uvy - 1; 
                            color = Vector3(std::pow(textures[0][(uvy * W + uvx) * 3] / 255., 2.2),
                                std::pow(textures[0][(uvy * W + uvx) * 3 + 1] / 255., 2.2),
                                std::pow(textures[0][(uvy * W + uvx) * 3 + 2] / 255., 2.2)); 















                        }

                    }
                }
            }

        }
        return intersect;
    }

    void loadTexture(const char* filename) {
        int W, H, C;
        unsigned char* texture = stbi_load(filename, &W, &H, &C, 3);
        Wtex.push_back(W); 
        Htex.push_back(H); 
        textures.push_back(texture);
    }

    std::vector<TriangleIndices> indices;
    std::vector<Vector3> vertices;
    std::vector<Vector3> normals;
    std::vector<Vector3> uvs;
    std::vector<Vector3> vertexcolors;
    std::vector<unsigned char*> textures;
    std::vector<int> Wtex, Htex;
    BoundingBox bb;
    Node* node;

};

class Scene { 
public:
    Scene() {};
    bool intersection(const Ray& r, Vector3& P, Vector3& N, Vector3& albedo, double& t, bool& Mirror, bool& transp, int& id) {
        t = 1E9;
        bool intersect = false;
        for (int i = 0; i < Objects.size(); i++) {
            Vector3 PObject, NObject, CoulObject;
            double tObject;
            if (Objects[i]->intersection(r, PObject, NObject, tObject, CoulObject) && tObject < t) { //s'il y a une intersection plus proche que celles existantes, alors on prend en compte celle-ci et pas les autres
                intersect = true;
                t = tObject;
                P = PObject;
                N = NObject;
                albedo = CoulObject;
                Mirror = Objects[i]->Mirror;
                transp = Objects[i]->Transparent;
                id = i;

            }

        }
        return intersect;
    }

    Vector3 GetColor(const Ray& r, int bounce, bool lastdiffusion) { 
        if (bounce > 5) return Vector3(0., 0., 0.); 
        else {
            double t;
            bool Mirror, transp;
            Vector3 P, N, albedo, color;
            int id;
            if (intersection(r, P, N, albedo, t, Mirror, transp, id)) {
                if (id == 0) { 
                    if (bounce == 0 || !lastdiffusion) {
                        double Rl = dynamic_cast<Sphere*>(Objects[0])->R;
                        return Vector3(I, I, I) / (4 * M_PI * M_PI * Rl * Rl);
                    }
                    else
                        return Vector3(0, 0, 0);
                }

                if (Mirror) { 
                    Vector3 reflectedDirection = r.u - 2 * DotProduct(r.u, N) * N; 
                    Ray reflectedRay(P + 0.00001 * N, reflectedDirection); 
                    return GetColor(reflectedRay, bounce + 1, false);
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
                            return GetColor(reflectedRay, bounce + 1, false);
                        }
                        Vector3 Tt = n1 / n2 * (r.u - DotProduct(r.u, N2) * N2); 
                        Vector3 Tn = -sqrt(angle) * N2;
                        Vector3 refractedDirection = Tt + Tn; 
                        Ray refractedRay(P - 0.0001 * N2, refractedDirection); 
                        return GetColor(refractedRay, bounce + 1, false);
                    }
                    else {
                        double Rl = dynamic_cast<Sphere*>(Objects[0])->R;
                        Vector3 Ol = dynamic_cast<Sphere*>(Objects[0])->O;
                        Vector3 w = RandomInUnitSphere((P - L).Normalize()); 
                        Vector3 xp = w * Rl + Ol;
                        Vector3 Pxp = xp - P;
                        double normPxp = sqrt(Pxp.NormSquared());
                        Pxp = Pxp / normPxp;
                        Vector3 Pshadow, Nshadow, albedoshadow;
                        double tshadow;
                        bool Mirrorshadow, transshadow;
                        int idshadow;
                        Ray Rayshadow(P + 0.00001 * N, Pxp);
                        if (intersection(Rayshadow, Pshadow, Nshadow, albedoshadow, tshadow, Mirrorshadow, transshadow, idshadow) && tshadow < normPxp - 0.0001) { //s'il y a intersection avant d'arriver à la source de lumière
                            color = Vector3(0., 0., 0.);
                        }
                        else {
                            color = I / (4 * M_PI * M_PI * Rl * Rl) * albedo / M_PI * std::max(0., DotProduct(N, Pxp)) * std::max(0., -DotProduct(w, Pxp)) / (normPxp * normPxp) / (std::max(0., -DotProduct((L - P).Normalize(), w)) / (M_PI * Rl * Rl));
                        }

                        
                        Vector3 random = RandomInUnitSphere(N);
                        Ray Rayrandom(P + 0.00001 * N, random);
                        color = color + TermByTermProduct(albedo, GetColor(Rayrandom, bounce + 1, true));


                    }
                }
            }
            return color;
        }
    }
    std::vector<Object*> Objects;
    Vector3 L;
    double I;

};


int main() {
    int W = 512;
    int H = 512;

    Vector3 C(0, 0, 30);
    Scene scene;
    scene.I = 5E9;
    scene.L = Vector3(-10, 20, 40);
    Sphere Slumiere(scene.L, 5, Vector3(1., 0.3, 0.2));
    Sphere S1(Vector3(0, 0, 0), 10, Vector3(1.0, 0.0, 0.0));
    Sphere Sgauche(Vector3(0, 0, 1000), 940, Vector3(1.0, 0.0, 0.0));
    Sphere Sdroite(Vector3(0, 0, -1000), 940, Vector3(0., 0.4, 1.));
    Sphere Shaut(Vector3(0, 1000, 0), 940, Vector3(0.7, 0., 0.7));
    Sphere Sbas(Vector3(0, -1000, 0), 990, Vector3(1., 1., 0.2));
    Sphere Splafond(Vector3(1000, 0, 0), 940, Vector3(0., 1., 0.));
    Sphere Ssol(Vector3(-1000, 0, 0), 940, Vector3(0.8, 0.4, 0.1));
    TriangleMesh m(Vector3(1., 1., 1.));

    m.readOBJ("gramophone1_01.obj");
    m.loadTexture("map_Gramophone1_BaseColor.png");
    m.loadTexture("map_Gramophone1_Metallic.png");
    m.loadTexture("map_Gramophone1_Normal.png");
    m.loadTexture("map_Gramophone1_Roughness.png");




    for (int i = 0; i < m.vertices.size(); i++) {
        m.vertices[i][0] *= 15;
        m.vertices[i][1] *= 15;
        m.vertices[i][2] *= 15;
        m.vertices[i][1] -= 10;
    }


    m.Pivot(m.node, 0, m.indices.size());

    scene.Objects.push_back(&Slumiere);
    scene.Objects.push_back(&Sgauche);
    scene.Objects.push_back(&Sdroite);
    scene.Objects.push_back(&Shaut);
    scene.Objects.push_back(&Sbas);
    scene.Objects.push_back(&Splafond);
    scene.Objects.push_back(&Ssol);
    scene.Objects.push_back(&m);




    double fov = 70 * M_PI / 180;

    int nmbRay = 100;

    std::vector<unsigned char> image(W * H * 3, 0);
#pragma omp parallel for schedule(dynamic,1)
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            Vector3 color(0, 0, 0);
            for (int k = 0; k < nmbRay; k++) {
                double u1 = uniform(engine);
                double u2 = uniform(engine);
                double x1 = 0.25 * cos(2 * M_PI * u1) * sqrt(-2 * log(u2));
                double x2 = 0.25 * sin(2 * M_PI * u1) * sqrt(-2 * log(u2));
                Vector3 u(j - W / 2 + x2 + 0.5, i - H / 2 + x1 + 0.5, -W / (2. * tan(fov / 2))); 
                u = u.Normalize();

                double u3 = uniform(engine);
                double u4 = uniform(engine);
                double x3 = 0.01 * cos(2 * M_PI * u3) * sqrt(-2 * log(u4));
                double x4 = 0.01 * sin(2 * M_PI * u3) * sqrt(-2 * log(u4));

                Vector3 Focalpoint = C + 55 * u;
                Vector3 RandOrigin = C + Vector3(x3, x4, 0);
                Vector3 Raydir = (Focalpoint - RandOrigin).Normalize();


                Ray r(RandOrigin, Raydir);

                color = color + scene.GetColor(r, 0, false);
            }
            color = color / nmbRay;


            color[0] = std::pow(color[0], 0.45);
            color[1] = std::pow(color[1], 0.45);
            color[2] = std::pow(color[2], 0.45);


            image[((H - i - 1) * W + j) * 3 + 0] = std::min(255., color[0]);
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(255., color[1]);
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(255., color[2]);
        }
    }
    stbi_write_png("Gramophone.png", W, H, 3, &image[0], 0);

    return 0;
}
