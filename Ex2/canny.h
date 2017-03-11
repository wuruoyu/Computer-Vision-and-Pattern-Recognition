#ifndef _CANNY_
#define _CANNY_
#define cimg_use_jpeg
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
    CImg<unsigned char> img; // Original Image
    CImg<unsigned char> edge; // image after edge detected
    CANNY can;
    int width;
    int height;
    // manual parameter
    float lowThreshold;
    float highThreshold;
    float gaussianKernelRadius;
    int gaussianKernelWidth;
    bool contrastNormalised;
    // fixed parameter
    const float gaussianCutOff = 0.005f;


public:
    canny(char const*);
    void setLowThreshold(float);
    void setHighThreshold(float);
    void setGaussianKernelRadius(float);
    void setContrastNormalised(float);
    void allocationErrorExit(unsigned char*, CANNY*);
    void run();
    CANNY *allocatebuffers(unsigned char*, int, int);
    void killbuffers(CANNY*);
    void normalizeContrast(unsigned char*, int, int);
    int computeGradients(CANNY*, float, int);
    void performHysteresis(CANNY*, int, int);
    void follow(CANNY*, int, int, int, int);
    float hypotenuse(float, float);
    float gaussian(float, float);
    void showEdgeDetected(unsigned char*);
    void rgb2bmp();
};

#endif