#ifndef _CANNY_
#define _CANNY_
#include "CImg.h"
#include <vector>

using namespace std;
using namespace cimg_library;

typedef struct
{
  unsigned char *data; /* input image */
  int width;           
  int height;
  int *idata;          /* output for edges */
  int *magnitude;      /* edge magnitude as detected by Gaussians */
  float *xConv;        /* temporary for convolution in x direction */
  float *yConv;        /* temporary for convolution in y direction */
  float *xGradient;    /* gradients in x direction, as detected by Gaussians */
  float *yGradient;    /* gradients in y direction, as detected by Gaussians */
} CANNY;

class canny {   
private:
    CImg<unsigned char> img; //Original Image
    CANNY can;
    float lowThreshold;
    float = highThreshold;
    float gaussianKernelRadius;
    int gaussianKernelWidth;
    bool contrastNormalised;


public:
    canny(String);
    void setLowThreshold(float);
    void setHighThreshold(float);
    void setGaussianKernelRadius(float);
    void setContrastNormalised(float);
    void allocationErrorExit(unsigned char*, CANNY*);
    void run();
    CANNY *allocatebuffers(unsigned char*, int, int);
    void killbuffers(CANNY*);
    void normalizeContrast(can->data, width, height);
    int computeGradients(CANNY*, float, int);
    void performHysteresis(CANNY*, int, int);
    void follow(CANNY*, int, int, int, int);
    void normalizeContrast(unsigned char*, int, int);
    float hypotenuse(float, float);
    float gaussian(float, float);
};

#endif