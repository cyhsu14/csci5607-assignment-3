//File parsing example
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <algorithm>
#include <vector>
#include "image.h"
// #include "pvector.h"
#include "utils.h"

#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "stb_image_write.h"

using namespace std;
PVector cam_pos,cam_dir,cam_up;
Color bg = Color();
typedef list<Sphere> SphereList;
typedef list<Material> MatList;
typedef list< pair<Color, PVector> > LightList;
typedef vector<PVector> VerticesList;
typedef vector<PVector> NormalsList;
typedef list<Triangle> TriangleList;
typedef list<SpotLight> SLightList;
SphereList s;
LightList ptl;
LightList drl;
SLightList spl;
VerticesList vl;
NormalsList nl;
TriangleList trl;
int objNum=0;
Color ambl = Color();
Color ptl_c = Color();
PVector ptl_p = PVector();
Material air = Material(Color(),Color(),Color(),0, Color(),1);
int max_depth = 5;

Color FindIntersection(PVector &pr, int depth, PVector &start, Material pre);
float GetArea(PVector v1, PVector v2, PVector v3){
    return CrossProduct(v1-v2,v3-v2).GetLength()*0.5;
}
Color GetColor(Hit hit, int depth, PVector &start, Material pre){
    PVector u,v;
    LightList::iterator iterl;
    SphereList::iterator iters;
    SLightList::iterator itersp;
    TriangleList::iterator itert;
    float cos_theta,t;
    Color out,c;
    bool is_block;

    // ambient light
    out = Product(ambl, hit.m.amb);

    PVector surface_normal;

    // point light
    PVector pc;
    PVector norm;

    for(iterl = ptl.begin();iterl != ptl.end();iterl++){
        is_block = false;
        u = (hit.hit_pt - iterl->second).Normalize(); // uni vector from hit pt to light
        float d_light= (hit.hit_pt - iterl->second).GetLength(); // distance to light

        for (iters = s.begin();iters != s.end();iters++){
            if(iters->p == hit.s.p) continue; // if same sphere, ignore

            pc = iters->p - iterl->second; // distance from light to sphere center
            cos_theta = InnerProduct(u, pc)/pc.GetLength(); // cos(theta)
            float d = sqrt(1 - cos_theta*cos_theta)*pc.GetLength(); // distance from hit pt to another sphere

            if(d < iters->r){ // calculate the point intersecting with other sphere between hit pt and light
                t = InnerProduct(u,iters->p ) - InnerProduct(u, iterl->second)  - sqrt(iters->r*iters->r - d*d);
                if(t<d_light && t>=0.01){
                    is_block=true;
                    break;
                }
            }
        }
        for(itert=trl.begin();itert!=trl.end();itert++){
            if(*itert == hit.t)continue;

            u = u*(-1); // ray is opposite of light direction
            // see if really intersect with triangle
            if(itert->n2 != -1 && itert->n3!=-1){ // normal triangle
                surface_normal = CrossProduct(itert->v1 - itert->v2, itert->v3 - itert->v2).Normalize();
            }
            else{ // triangle
                surface_normal = itert->n;
            }
            t = InnerProduct(itert->v2 - hit.hit_pt, surface_normal)/InnerProduct(u,surface_normal);

            if(t<d_light && t>=0.01){
                float gamma = GetArea(hit.hit_pt+u*t,itert->v1,itert->v2)/GetArea(itert->v1,itert->v2,itert->v3);
                float alpha = GetArea(hit.hit_pt+u*t,itert->v2,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                float beta = GetArea(hit.hit_pt+u*t,itert->v1,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                if((alpha+beta+gamma)<= 1.01){
                    is_block = true;

                    break;
                }
            }
        }

        if(hit.s.r != -1) // if is sphere
            surface_normal = (hit.hit_pt - hit.s.p).Normalize();
        else { // is triangle
            if(hit.t.n2 != -1 && hit.t.n3!=-1){ // normal triangle
                surface_normal = hit.n;
            }
            else{ // triangle
                cos_theta = InnerProduct(hit.t.n, (iterl->second - hit.hit_pt).Normalize());
                surface_normal = cos_theta>=0 ? hit.t.n : hit.t.n*(-1);
            }
        }
        if(!is_block){
            // diffuse light for all point light
            PVector ray_n = (iterl->second - hit.hit_pt).Normalize();
            cos_theta =  InnerProduct(surface_normal, ray_n)>0 ? InnerProduct(surface_normal, ray_n) : 0;

            PVector current = hit.hit_pt - iterl->second; // vector from current pixel to light
            Color intensity = iterl->first*(1 / (current.GetLength()*current.GetLength()));
            c = Product(hit.m.dif,intensity);
            out = out + c*cos_theta;
            // specular light
            ray_n = ray_n * (-1);

            PVector reflected = ray_n - surface_normal*InnerProduct(surface_normal, ray_n)*2;
            PVector viewing = (start - hit.hit_pt).Normalize();
            float cos_alpha = InnerProduct(reflected.Normalize(), viewing);
            cos_alpha = cos_alpha>0 ? cos_alpha : 0;
            out = out + Product(intensity, hit.m.spe)*pow(cos_alpha, hit.m.ns);
            // surface_normal.print();
            // out=Color(1,0,0);
        }



    }

    // directional light
    for(iterl = drl.begin();iterl != drl.end();iterl++){
        is_block = false;

        for (iters = s.begin();iters != s.end();iters++){

            pc = iters->p - hit.hit_pt;
            cos_theta = InnerProduct(iterl->second*(-1), pc)/pc.GetLength(); // cos(theta)
            float d = sqrt(1 - cos_theta*cos_theta)*pc.GetLength(); // distance from hit pt to another sphere
            if(d < iters->r){ // calculate the point intersecting with other sphere between hit pt and light
                t = InnerProduct(iterl->second*(-1),iters->p ) - InnerProduct(iterl->second*(-1), hit.hit_pt)  - sqrt(iters->r*iters->r - d*d);
                if(t>=0.01){
                    is_block=true;
                    break;
                }
            }
        }
        for(itert=trl.begin();itert!=trl.end();itert++){
            // check normal
            // see if really intersect with triangle

            if(itert->n2 != -1 && itert->n3!=-1){// normal triangle
                surface_normal = CrossProduct(itert->v1 - itert->v2, itert->v3 - itert->v2).Normalize();
            }
            else{ // triangle
                surface_normal = hit.t.n;
            }
            t = InnerProduct(itert->v2 - hit.hit_pt, surface_normal) /InnerProduct(iterl->second*(-1),surface_normal);

            if(t>=0.01){
                float gamma = GetArea(hit.hit_pt+iterl->second*(-t),itert->v1,itert->v2)/GetArea(itert->v1,itert->v2,itert->v3);
                float alpha = GetArea(hit.hit_pt+iterl->second*(-t),itert->v2,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                float beta = GetArea(hit.hit_pt+iterl->second*(-t),itert->v1,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                if((alpha+beta+gamma)<= 1.01){
                    is_block = true;
                    break;
                }
            }

        }
        if(hit.s.r != -1) // if is sphere
            surface_normal = (hit.hit_pt - hit.s.p).Normalize();
        else { // is triangle
            if(hit.t.n2 != -1 && hit.t.n3!=-1){// normal triangle
                surface_normal = hit.n;
            }
            else{ // triangle
                cos_theta = InnerProduct(hit.t.n, iterl->second*(-1));
                surface_normal = cos_theta>=0 ? hit.t.n : hit.t.n*(-1);
            }
        }
        if(!is_block){
            // diffuse light for all point light
            c = Product(hit.m.dif,iterl->first);

            PVector ray_n = (iterl->second*(-1)).Normalize();
            cos_theta =  InnerProduct(surface_normal, ray_n)>0 ? InnerProduct(surface_normal, ray_n) : 0;
            out = out + c*cos_theta;

            // specular light
            ray_n = iterl->second;

            PVector reflected = ray_n - surface_normal*InnerProduct(surface_normal, ray_n)*2;
            PVector viewing = (start - hit.hit_pt).Normalize();
            float cos_alpha = InnerProduct(reflected, viewing)/reflected.GetLength();
            cos_alpha = cos_alpha>0 ? cos_alpha : 0;
            out = out + Product(iterl->first, hit.m.spe)*pow(cos_alpha, hit.m.ns);

        }

    }

    // spot light
    for(itersp=spl.begin();itersp!=spl.end();itersp++){
        is_block = false;
        //sphere
        float d_light= (hit.hit_pt - itersp->p).GetLength(); // distance to light
        Color intensity = itersp->c*(1 / (d_light*d_light));
        cos_theta = InnerProduct((hit.hit_pt - itersp->p).Normalize(), itersp->d);
        float spe_coe=1;
        if(cos_theta > cos(itersp->ang1*M_PI/180)){ // angle < ang1: point light
            c = Product(hit.m.dif,intensity);

        }
        else if(cos_theta < cos(itersp->ang1*M_PI/180) && cos_theta > cos(itersp->ang2*M_PI/180)){
            // between ang1, ang2: fall off
            c = Product(hit.m.dif,intensity*((itersp->ang2-(acos(cos_theta)*180/M_PI))/(itersp->ang2 - itersp->ang1)));
            spe_coe = (itersp->ang2-(acos(cos_theta)*180/M_PI))/(itersp->ang2 - itersp->ang1);

        }
        else {
            is_block=true;
        }
        for (iters = s.begin();iters != s.end();iters++){
            pc = iters->p - hit.hit_pt;
            cos_theta = InnerProduct(itersp->d*(-1), pc)/pc.GetLength(); // cos(theta)

            float d = sqrt(1 - cos_theta*cos_theta)*pc.GetLength(); // distance from (hit_pt - light.p) to another sphere
            if(d < iters->r){ // calculate the point intersecting with other sphere between hit pt and light
                t = InnerProduct(itersp->d*(-1),iters->p ) - InnerProduct(itersp->d*(-1), hit.hit_pt)  - sqrt(iters->r*iters->r - d*d);
                if(t<d_light && t>=0.01){
                    is_block=true;
                    break;
                }
            }
        }

        //triangle
        for(itert=trl.begin();itert!=trl.end();itert++){
            // check normal
            // see if really intersect with triangle

            if(itert->n2 != -1 && itert->n3!=-1){// normal triangle
                surface_normal = CrossProduct(itert->v1 - itert->v2, itert->v3 - itert->v2).Normalize();
            }
            else{ // triangle
                surface_normal = hit.t.n;
            }

            t = InnerProduct(itert->v2 - hit.hit_pt, surface_normal) /InnerProduct(itersp->d*(-1),surface_normal);

            if(t<d_light && t>=0.01){
                float gamma = GetArea(hit.hit_pt+itersp->d*(-t),itert->v1,itert->v2)/GetArea(itert->v1,itert->v2,itert->v3);
                float alpha = GetArea(hit.hit_pt+itersp->d*(-t),itert->v2,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                float beta = GetArea(hit.hit_pt+itersp->d*(-t),itert->v1,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
                if((alpha+beta+gamma)<= 1.01){
                    is_block = true;
                    break;
                }
            }

        }
        if(hit.s.r != -1) // if is sphere
            surface_normal = (hit.hit_pt - hit.s.p).Normalize();
        else { // is triangle

            if(hit.t.n2 != -1 && hit.t.n3!=-1){// normal triangle
                surface_normal = hit.n;
            }
            else{ // triangle
                cos_theta = InnerProduct(hit.t.n, iterl->second*(-1));
                surface_normal = cos_theta>=0 ? hit.t.n : hit.t.n*(-1);
            }
        }
        if(!is_block){
            // diffuse light for all point light
            PVector ray_n = (itersp->p - hit.hit_pt).Normalize();
            float theta =  InnerProduct(surface_normal, ray_n)>0 ? InnerProduct(surface_normal, ray_n) : 0;
            out = out + c*theta;

            // specular light
            ray_n = ray_n * (-1);
            PVector reflected = ray_n - surface_normal*InnerProduct(surface_normal, ray_n)*2;
            PVector viewing = (start - hit.hit_pt).Normalize();
            float cos_alpha = InnerProduct(reflected, viewing)/reflected.GetLength();
            cos_alpha = cos_alpha>0 ? cos_alpha : 0;

            out = out + Product(intensity*spe_coe, hit.m.spe)*pow(cos_alpha, hit.m.ns);


        }
    }
    /**/
    if(hit.s.r == -1) // if is triangle
        if(hit.t.n2 != -1 && hit.t.n3!=-1){// normal triangle
            surface_normal = hit.n;
        }
        else{ // triangle
            cos_theta = InnerProduct(hit.t.n, (hit.hit_pt - start).Normalize());
            surface_normal = cos_theta>=0 ? hit.t.n : hit.t.n*(-1);
        }
    else { // is sphere
        surface_normal = (hit.hit_pt - hit.s.p).Normalize();

    }
    // Mirror
    PVector ray_n = (hit.hit_pt - start).Normalize();
    ray_n = (ray_n - surface_normal*2*InnerProduct(ray_n, surface_normal)).Normalize();
    // printf("%f\n", InnerProduct(surface_normal, ray_n - surface_normal*InnerProduct(ray_n,surface_normal)));
    // assert(InnerProduct(surface_normal, ray_n - surface_normal*InnerProduct(ray_n,surface_normal)) <= 0.1 );
    c = FindIntersection(ray_n, depth+1, hit.hit_pt, pre);

    c = Product(c, hit.m.spe);
    out = out + c;

    // Refraction
    ray_n = (start - hit.hit_pt).Normalize();
    float ior_out = (pre == hit.m) ? 1 : hit.m.ior;
    float cos_in = InnerProduct(ray_n, surface_normal);
    float sin_out = sqrt(1 - cos_in*cos_in)*pre.ior/ior_out;
    float cos_out = sqrt(1 - sin_out*sin_out);
    PVector T = (surface_normal*(cos_in*pre.ior/ior_out - cos_out) - ray_n*(pre.ior/ior_out)).Normalize();

    c = FindIntersection(T, depth+1, hit.hit_pt, hit.m);
    c = Product(c, hit.m.tran);

    out = out + c;

    return out;
}
Color FindIntersection(PVector &pr,int depth ,PVector &start, Material pre){
    if(depth > max_depth) return Color(0,0,0);

    SphereList::iterator iters;
    TriangleList::iterator itert;
    Color out = bg;
    float t;
    float cos_theta;
    float t_smallest = -1;
    Hit hit = Hit();
    for(iters = s.begin(); iters!=s.end(); iters++){
        // lightening
        PVector pc = iters->p - start; // distance from start to sphere center
        cos_theta = InnerProduct(pr, pc)/(pr.GetLength()*pc.GetLength()); // cos(theta)
        float d = sqrt(1 - cos_theta*cos_theta)*pc.GetLength(); // distance from ray to center

        if(d < iters->r){
            t = InnerProduct(pr,iters->p ) - InnerProduct(pr, start)  - sqrt(iters->r*iters->r - d*d);
            if((t_smallest == -1 || t_smallest>t) && t>=0.01){
                t_smallest = t;
                hit.s = *iters;
                hit.n = PVector(0,0,0);
                hit.m = iters->m;
            }
        }

    }
    // triangle
    for(itert=trl.begin();itert!=trl.end();itert++){
        // check normal
        PVector norm;
        if(itert->n2 == -1 && itert->n3==-1){ // not normal triangle
            norm = itert->n;
        }
        else {
            norm = CrossProduct(itert->v1 - itert->v2, itert->v3 - itert->v2).Normalize();
        }

        t = -InnerProduct(start-itert->v2, norm)/InnerProduct(pr,norm);

        // see if really intersect with triangle
        float gamma = GetArea(start+pr*t,itert->v1,itert->v2)/GetArea(itert->v1,itert->v2,itert->v3);
        float alpha = GetArea(start+pr*t,itert->v2,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);
        float beta = GetArea(start+pr*t,itert->v1,itert->v3)/GetArea(itert->v1,itert->v2,itert->v3);

        if((t_smallest == -1 || t_smallest>t) && t>=0.01 && (alpha+beta+gamma)<= 1.01){
            t_smallest = t;
            hit.s.r = -1;
            hit.t = *itert;

            if(itert->n2 != -1 && itert->n3!=-1){
                hit.n = (itert->n*alpha + (nl.at(itert->n2))*beta + (nl.at(itert->n3))*gamma).Normalize();
            }
            else hit.n = PVector(0,0,0);
            hit.m = itert->m ;
            // hit.m.print();
        }
    }

    // get hit point and information
    if(t_smallest != -1){
        hit.hit_pt = start + pr*t_smallest;
        out = GetColor(hit,depth,start, pre);
        // printf("%d\n", depth);
    }

    out.SetClamp();

    return out;

}



int main(int argc, char*argv[]){
    string line;

    // string fileName = "input/ambient_sphere.scn";
    string fileName;
    fileName = argv[1];
    // string fileName = "input/spheres2.scn";

    // open the file containing the scene description
    ifstream input(fileName);

    // check for errors in opening the file
    if(input.fail()){
        cout << "Can't open file '" << fileName << "'" << endl;
        return 0;
    }

    // determine the file size (this is optional -- feel free to delete the 6 lines below)
    streampos begin,end;
    begin = input.tellg();
    input.seekg(0, ios::end);
    end = input.tellg();
    cout << "File '" << fileName << "' is: " << (end-begin) << " bytes long.\n\n";
    input.seekg(0, ios::beg);

    // set default number
    float px=0,py=0,pz=0,dx=0,dy=0,dz=1,ux=0,uy=1,uz=0,ha=45;
    int width=640,height=480;
    string outfileName = "raytraced.bmp";
    MatList mat;
    mat.push_back(Material());
    // amblient light and point light
    // float spl_r=0,spl_g=0,spl_b=0;

    float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, ior;
    int n=5;
    int max_vertices = -1;
    int max_normals = -1;
    float x,y,z,r;
    // Sphere s;
    Image *img = new Image(width, height);


    //Loop through reading each line
    string command;
    while(input >> command) { //Read first word in the line (i.e., the command type)

        if (command[0] == '#'){
            getline(input, line); //skip rest of line
            cout << "Skipping comment: " << command  << line <<  endl;
            continue;
        }


        if (command == "sphere"){ //If the command is a sphere command
           input >> x >> y >> z >> r;
           // matttt

           s.push_back(Sphere(PVector(x,y,z), r, mat.back()));

           printf("Sphere as position (%f,%f,%f) with radius %f\n",x,y,z,r);
        }
        else if (command == "background"){ //If the command is a background command
            float bg_r, bg_g, bg_b;
            input >> bg_r >> bg_g >> bg_b;
            bg = Color(bg_r,bg_g,bg_b);
            printf("Background color of (%f,%f,%f)\n",bg_r,bg_g,bg_b);
        }
        else if (command == "output_image"){ //If the command is an output_image command
            input >> outfileName;
            printf("Render to file named: %s\n", outfileName.c_str());
        }
        else if(command == "camera"){
            input >> px >> py >> pz >> dx >> dy >> dz >> ux >> uy >> uz >> ha;
            printf("Camera position: %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",px,py,pz,dx,dy,dz,ux,uy,uz,ha);
        }
        else if(command == "film_resolution"){
            input >> width >> height;
            img = new Image(width, height);
        }
        else if(command == "material"){
            input >> ar >> ag >> ab >> dr >> dg >> db >> sr >> sg >> sb >> ns >> tr >> tg >> tb >> ior;
            mat.push_back(Material(Color(ar, ag, ab), Color(dr, dg, db), Color(sr, sg, sb), ns, Color(tr, tg,tb), ior));
        }
        else if(command == "directional_light"){
            float ptr, ptg, ptb, ptx, pty, ptz;
            input >> ptr >> ptg >> ptb >> ptx >> pty >> ptz;
            drl.push_back(make_pair(Color(ptr,ptg,ptb), PVector(ptx,pty,ptz).Normalize()));
        }
        else if(command == "point_light"){
            float ptr, ptg, ptb, ptx, pty, ptz;
            input >> ptr >> ptg >> ptb >> ptx >> pty >> ptz;
            ptl.push_back(make_pair(Color(ptr,ptg,ptb), PVector(ptx,pty,ptz)));
        }
        else if(command == "ambient_light"){
            float aml_r, aml_g, aml_b;
            input >> aml_r >> aml_g >> aml_b;
            ambl = Color(aml_r,aml_g,aml_b);
        }
        else if(command == "spot_light"){
            float r, g, b, px, py, pz, dx, dy, dz, angle1, angle2;
            input >> r >> g >> b >> px >> py >> pz >> dx >> dy >> dz >> angle1 >> angle2;
            spl.push_back(SpotLight(Color(r,g,b),PVector(px,py,pz),PVector(dx,dy,dz).Normalize(),angle1,angle2));

        }
        else if(command == "max_depth"){
            input >> max_depth;
        }
        else if(command == "max_vertices"){
            input >> max_vertices;
        }
        else if(command == "max_normals"){
            input >> max_normals;
        }
        else if(command == "vertex"){
            if(max_vertices == -1){
                printf("Max vertices not defined. Please check your input.\n");
                break;
            }
            input >> x >> y >> z;

            vl.push_back(PVector(x,y,z));

        }
        else if(command == "normal"){
            if(max_normals == -1){
                printf("Max normals not defined. Please check your input.\n");
                break;
            }
            input >> x >> y >> z;
            nl.push_back(PVector(x,y,z).Normalize());
        }
        else if(command == "triangle"){
            int i1, i2, i3;
            input >> i1 >> i2 >> i3;
            PVector v1 = (vl.at(i1));
            PVector v2 = (vl.at(i2));
            PVector v3 = (vl.at(i3));
            PVector normal = CrossProduct(v1-v2,v3-v2).Normalize();
            trl.push_back(Triangle(v1,v2,v3,normal,-1,-1,mat.back()));
        }
        else if(command == "normal_triangle"){
            int i1, i2, i3;
            input >> i1 >> i2 >> i3;
            PVector v1 = (vl.at(i1));
            PVector v2 = (vl.at(i2));
            PVector v3 = (vl.at(i3));
            input >> i1 >> i2 >> i3;
            PVector n1 = (nl.at(i1));

            trl.push_back(Triangle(v1,v2,v3,n1,i2,i3,mat.back()));
        }
        else {
          getline(input, line); //skip rest of line
          cout << "WARNING. Do not know command: " << command << endl;
        }

    }

    printf("Normal num: %d\n", nl.size());
    printf("Vertices num: %d\n", vl.size());
    // setting
    cam_pos = PVector(px,py,pz);
    cam_dir = PVector(dx,dy,dz).Normalize(); // from plane to camera
    cam_up = PVector(ux,uy,uz).Normalize();
    // make sure up is orthogonal to direction
    cam_up = (cam_up - cam_dir*InnerProduct(cam_up, cam_dir)).Normalize();


    PVector left = CrossProduct(cam_dir,cam_up).Normalize();
    PVector left_top = cam_pos + cam_dir*(1/tan(ha*M_PI/180))*(0.5*height) + cam_up*(0.5*height) + left*(0.5*width);
    left_top = left_top - PVector(0.5,0.5,0);


    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            PVector v = left_top - cam_up*j - left*i;
            v = (v - cam_pos).Normalize();
            Color col = FindIntersection(v ,0, cam_pos, air);
            Pixel p = Pixel(col.r*255, col.g*255,col.b*255);

            img->SetPixel(i,j,p);

        }
    }
    img->Write(&outfileName[0u]);

    return 0;
}
