#include <math.h>
struct Color{
    float r,g,b; // should be in [0,1]
    Color (float r_=0, float g_=0, float b_=0): r(r_),g(g_),b(b_) {}
    Color SetClamp(){
        r = r>1 ? 1 : (r<0 ? 0 : r);
        g = g>1 ? 1 : (g<0 ? 0 : g);
        b = b>1 ? 1 : (b<0 ? 0 : b);
        return Color(r,g,b);
    }
    void print(){
        printf("%f %f %f\n",r,g,b );
    }
};
Color operator+ (const Color& a, const Color& b){
    return Color(a.r+b.r, a.g+b.g, a.b+b.b);
}

Color operator- (const Color& a, const Color& b){
    return Color(a.r-b.r, a.g-b.g, a.b-b.b);
}

Color operator* (const Color&a, float b){
    return Color(a.r*b, a.g*b, a.b*b);
}

Color Product(const Color &a, const Color &b){
    return Color(a.r*b.r, a.g*b.g, a.b*b.b);
}
bool operator== (const Color&a, const Color&b){
    return (a.r==b.r) && (a.g==b.g) && (a.b==b.b);
}


struct Material{
    Color amb, dif, spe, tran;
    float ns, ior;

    Material (Color amb_=Color(), Color dif_=Color(1,1,1), Color spe_=Color(),
              float ns_=5, Color tran_=Color(0,0,0), float ior_=1):
              amb(amb_), dif(dif_), spe(spe_) ,ns(ns_), tran(tran_), ior(ior_) {}
    void print(){
        amb.print();
        dif.print();
        spe.print();
        printf("%f\n", ns);
        tran.print();
        printf("%f\n", ior);
    }
};
Material operator*(const float a,const Material &b){
    return Material(b.amb*a, b.dif*a, b.spe*a, b.ns*a, b.tran*a, b.ior*a);
}
Material operator+(const Material &a,const Material &b){
    return Material(a.amb+b.amb, a.dif+b.dif, a.spe+b.spe, a.ns+b.ns, a.tran+b.tran, a.ior+b.ior);
}
bool operator==(const Material &a,const Material &b){
    return (a.amb==b.amb) && (a.dif==b.dif) && (a.spe==b.spe) && (a.ns==b.ns) && (a.tran==b.tran) && (a.ior==b.ior);
}

struct PVector
{
    float x, y, z;
    PVector (float x_=0, float y_=0, float z_=0): x(x_), y(y_), z(z_) {}


    PVector Normalize(){
        float length = sqrt(x*x + y*y + z*z);
        x = x/length;
        y = y/length;
        z = z/length;
        return PVector(x,y,z);
    }
    float GetLength(){
        return sqrt(x*x + y*y + z*z);
    }
    void print(){
        printf("%f %f %f\n",x,y,z );
    }

};

PVector CrossProduct(const PVector &a, const PVector &b){
    return PVector(a.y*b.z - a.z*b.y,
                   a.z*b.x - a.x*b.z,
                   a.x*b.y - a.y*b.x);
}

float InnerProduct(const PVector &a, const PVector &b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

PVector operator+ (const PVector& a, const PVector& b){
    return PVector(a.x+b.x, a.y+b.y, a.z+b.z);
}

PVector operator- (const PVector& a, const PVector& b){
    return PVector(a.x-b.x, a.y-b.y, a.z-b.z);
}

PVector operator* (const PVector&a, float b){
    return PVector(a.x*b, a.y*b, a.z*b);
}
bool operator==(const PVector& a, const PVector& b){
    return (a.x==b.x) && (a.y==b.y) && (a.z==b.z);
}


struct Sphere{
    float r;
    PVector p;
    Material m;
    Sphere (PVector p_=PVector(), float r_=0, Material m_=Material()
            ): p(p_), r(r_), m(m_) {}
};
struct Triangle{
    PVector v1, v2, v3;
    PVector n; // normal
    int n2, n3;
    Material m;
    Triangle( PVector v1_=PVector(), PVector v2_=PVector(), PVector v3_=PVector(), PVector n_=PVector(),
              int n2_=-1, int n3_=-1, Material m_=Material()):
             v1(v1_), v2(v2_), v3(v3_), n(n_), n2(n2_), n3(n3_), m(m_){}
    // PVector n1, n2, n3; // normal
    // Triangle( PVector v1_=PVector(), PVector v2_=PVector(), PVector v3_=PVector(), PVector n1_=PVector(),
    //           PVector n2_=PVector(), PVector n3_=PVector()):
    //          v1(v1_), v2(v2_), v3(v3_), n1(n1_), n2(n2_), n3(n3_) {}
};
bool operator==(const Triangle& t1, const Triangle&t2){
    return (t1.v1==t2.v1) && (t1.v2==t2.v2) && (t1.v3==t2.v3);
}
struct Hit{
    Sphere s;
    Triangle t;
    PVector hit_pt;
    Material m;
    PVector n; // normal for normal triangle
    Hit(Sphere s_=Sphere(), Triangle t_=Triangle(), PVector h_=PVector(),Material m_=Material(), PVector n_=PVector()) :
        s(s_), t(t_), hit_pt(h_), m(m_), n(n_){}
};
struct SpotLight{
    Color c;
    PVector p, d; // position and direction
    float ang1, ang2;
    SpotLight(Color c_=Color(), PVector p_=PVector(), PVector d_=PVector(), float a1_=0, float a2_=0):
              c(c_), p(p_), d(d_), ang1(a1_), ang2(a2_) {}
};
