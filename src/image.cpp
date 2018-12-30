#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <random>
#include <iostream>
#include <algorithm>
#include <list>
using namespace std;
/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 255;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){

	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;

    data.raw = new uint8_t[num_pixels*4];

    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){

	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);

	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}


	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;

}

Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){

	int lastc = strlen(fname);
	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}


void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    /* WORK HERE */
    Pixel p;
    int x,y;
    int x0,y0,x1,y1;
    double r,g,b;
    if( sampling_method == IMAGE_SAMPLING_POINT){
        p = GetPixel((int)floor(u),(int)floor(v));
    }
    else if( sampling_method == IMAGE_SAMPLING_BILINEAR){
        x0 = floor(u);
        x1 = x0+1;
        if(x1==Width())x1-=1;
        y0 = floor(v);
        y1 = y0+1;
        if(y1==Height())y1-=1;
        // area: 0 1    (x0,y0) (x1,y0)
        //       2 3    (x0,y1) (x1,y1)
        double area0, area1,area2,area3;
        area0 = (u-x0)*(v-y0);
        area1 = (x1-u)*(v-y0);
        area2 = (u-x0)*(y1-v);
        area3 = (x1-u)*(y1-v);

        if(u!=x0 && x1!=u && y0!=v && v!=y1){
            r = GetPixel(x1,y1).r*area0+GetPixel(x0,y1).r*area1+GetPixel(x1,y0).r*area2+GetPixel(x0,y0).r*area3;
            g = GetPixel(x1,y1).g*area0+GetPixel(x0,y1).g*area1+GetPixel(x1,y0).g*area2+GetPixel(x0,y0).g*area3;
            b = GetPixel(x1,y1).b*area0+GetPixel(x0,y1).b*area1+GetPixel(x1,y0).b*area2+GetPixel(x0,y0).b*area3;
        }
        else{
            if(y0==v){
                r = GetPixel(x1,y0).r*(u-x0)+GetPixel(x0,y0).r*(x1-u);
                g = GetPixel(x1,y0).g*(u-x0)+GetPixel(x0,y0).g*(x1-u);
                b = GetPixel(x1,y0).b*(u-x0)+GetPixel(x0,y0).b*(x1-u);
            }
            else if(y1==v){
                r = GetPixel(x1,y1).r*(u-x0)+GetPixel(x0,y1).r*(x1-u);
                g = GetPixel(x1,y1).g*(u-x0)+GetPixel(x0,y1).g*(x1-u);
                b = GetPixel(x1,y1).b*(u-x0)+GetPixel(x0,y1).b*(x1-u);
            }
            else if(x0==u){
                r = GetPixel(x0,y1).r*(v-y0)+GetPixel(x0,y0).r*(y1-v);
                g = GetPixel(x0,y1).g*(v-y0)+GetPixel(x0,y0).g*(y1-v);
                b = GetPixel(x0,y1).b*(v-y0)+GetPixel(x0,y0).b*(y1-v);
            }
            else if(x1==u){
                r = GetPixel(x1,y1).r*(v-y0)+GetPixel(x1,y0).r*(y1-v);
                g = GetPixel(x1,y1).g*(v-y0)+GetPixel(x1,y0).g*(y1-v);
                b = GetPixel(x1,y1).b*(v-y0)+GetPixel(x1,y0).b*(y1-v);
            }
        }
        p.r = (int)(r+0.5);
        p.g = (int)(g+0.5);
        p.b = (int)(b+0.5);
        p.SetClamp(p.r,p.g,p.b);
        if((x0==u&&y0==v) || (x1==u&&y0==0) || (x1==u&&y1==v) || (x1==u&&y1==v))p = GetPixel((int)u,(int)v);

    }
    else if( sampling_method == IMAGE_SAMPLING_GAUSSIAN){
        double sigma = 1;
        double sum = 0.0;
        int n=4;
        double** kernel = new double*[n];
        for(x = 0;x < n; x++) kernel[x] = new double[n];
        // generate kernel
        for (x = 0; x < n; x++) {
            for (y = 0; y < n; y++) {
                kernel[x][y] = exp(-(x*x+y*y)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
                sum += kernel[x][y];
            }
        }
        // normalising
        for (x = 0; x < n; x++) {
            for (y = 0; y < n; y++) {
                kernel[x][y] /= sum;
            }
        }

        // multiply with gaussian kernel

        if(floor(u) == 0){ x0 = 0; x1 = x0 + 3;}
        else if(floor(u) == Width()-1) {x1 = floor(u); x0 = x1 - 3;}
        else if( floor(u) == Width()-2) {x0 = floor(u)-2; x1 = x0 + 3;}
        else {x0 = floor(u)-1; x1 = x0 +3;}
        // printf("test\n" );
        if(floor(v) == 0){ y0 = 0; y1 = y0 + 3;}
        else if(floor(v) == Height()-1) {y1 = floor(v); y0 = y1 - 3;}
        else if( floor(v) == Height()-2) {y0 = floor(v)-2; y1 = y0 + 3;}
        else {y0 = floor(v)-1; y1 = y0 +3;}
        // printf("test2\n" );
        // printf("%d %d %d %d\n",x0,x1,y0,y1 );
        r=0; g=0; b=0;
        for(x = x0; x <= x1; x++){
            for(y = y0; y <= y1; y++){
                r+= GetPixel(x,y).r*kernel[x-x0][y-y0];
                g+= GetPixel(x,y).g*kernel[x-x0][y-y0];
                b+= GetPixel(x,y).b*kernel[x-x0][y-y0];
            }
        }
        p.r = (int)(r+0.5);
        p.g = (int)(g+0.5);
        p.b = (int)(b+0.5);
        p.SetClamp(p.r,p.g,p.b);

    }
	return p;
}
