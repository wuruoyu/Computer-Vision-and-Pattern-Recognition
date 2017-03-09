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
    Mat grayscaled; // Grayscale
    Mat gFiltered; // Gradient
    Mat sFiltered; //Sobel Filtered
    Mat angles; //Angle Map
    Mat non; // Non-maxima supp.
    Mat thres; //Double threshold and final

public:
    canny(String); //Constructor
	Mat toGrayScale();
	vector<vector<double>> createFilter(int, int, double); //Creates a gaussian filter
	Mat useFilter(Mat, vector<vector<double>>); //Use some filter
    Mat sobel(); //Sobel filtering
    Mat nonMaxSupp(); //Non-maxima supp.
    Mat threshold(Mat, int, int); //Double threshold and finalize picture

    // code0's interface
    cannyparam();
};

#endif