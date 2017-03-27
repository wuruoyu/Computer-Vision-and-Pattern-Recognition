#ifndef _CANNY_
#define _CANNY_
#define cimg_use_jpeg                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
#include "CImg.h"
#include <vector>

using namespace std;
using namespace cimg_library;

class canny {   
private:
    CImg<unsigned char> colorImg;
    CImg<unsigned char> img; // Original Image
    CImg<unsigned char> edge; // image after edge detected
    int width;
    int height;

    unsigned char *data; /* input image */
    int *idata;          /* output for edges */
    int *magnitude;      /* edge magnitude as detected by Gaussians */
    float *xConv;        /* temporary for convolution in x direction */
    float *yConv;        /* temporary for convolution in y direction */
    float *xGradient;    /* gradients in x direction, as detected by Gaussians */
    float *yGradient;    /* gradients in y direction, as detected by Gaussians */

    // manual parameter
    float lowThreshold;
    float highThreshold;
    int gaussianKernelWidth;
    float gaussianKernelRadius;
    bool contrastNormalised;

    // fixed parameter
    const float gaussianCutOff = 0.005f;
    const float magnitudeScale = 100.0f;
    const float magnitudeLimit = 1000.f;


public:
    canny(CImg<unsigned char>);
    void setLowThreshold(float);
    void setHighThreshold(float);
    void setGaussianKernelWidth(float);
    void setGaussianKernelRadius(float);
    void setContrastNormalised(bool);
    void run();
    void allocatebuffers();
    void normalizeContrast();
    bool computeGradients();
    void performHysteresis();
    void follow(int, int, int, int);
    float hypotenuse(float, float);
    float gaussian(float, float);
    void showEdgeDetected();
    void rgb2gray();
    float ffabs(int);
    void killbuffers();
    CImg<unsigned char> getImg();
    CImg<unsigned char> getEdge();
};

#endif