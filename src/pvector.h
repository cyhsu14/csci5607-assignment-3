#include <math.h>

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

};

PVector CrossProduct(const PVector &a, const PVector &b){
    return PVector(a.y*b.z - a.z*b.y,
                   a.z*b.x - a.x*b.z,
                   a.x-b.y - a.y*b.x);
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
