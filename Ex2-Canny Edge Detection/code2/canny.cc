/****** Canny step-edge detection ******************************************
*
* compile: 
*   biwicompile canny.cc
*
* execute: 
*   canny -i ~buc/pic/zurlim.rgi --sigma 1.5 --threshold 7 --nms nms.png
*
*/

#include <cmath>
#include <string>
#include <iostream>
#include "gauss_filter.h"    // command line parameter reading
#include "CImg.h"         // C++ image library
using namespace cimg_library;
using namespace std;

/** CannyDiscrete.
* 
* Version with discrete implementation of the non-maximum-suppression.
* 
* During the nms, the gradient and its orientation are treated discrete,
* e.g., only four directions are considered: E, NE, N, SE.
*/
void CannyDiscrete(CImg<float> in, float sigma, float threshold,
				   CImg<float> &outSmooth, CImg<float> &outGradient, 
				   CImg<float> &outOrientation, CImg<float> &outThreshold, 
				   CImg<float> &outNMS) {
					   const int nx = in._width;
					   const int ny = in._height;

					   /************ initialize memory ************/
					   outGradient = in; outGradient.fill(0.0f);
					   CImg<int> dirmax(outGradient);
					   CImg<float> derivative[4];
					   for(int i = 0; i < 4; i++) { derivative[i] = outGradient; }
					   outOrientation = outGradient;
					   outThreshold = outGradient;
					   outNMS = outGradient;

					   /************** smoothing the input image ******************/
					   CImg<float> filter;
					   gauss_filter(filter, sigma, 0);
					   outSmooth = in.get_convolve(filter).convolve(filter.get_transpose());


					   /************ loop over all pixels in the interior **********************/
					   float fct = 1.0 / (2.0*sqrt(2.0f));
					   for (int y = 1; y < ny-1; y++) {
						   for(int x = 1; x < nx-1; x++) {
							   //***** compute directional derivatives (E,NE,N,SE) ****//
							   float grad_E = (outSmooth(x+1,y  ) - outSmooth(x-1,y  ))*0.5; // E
							   float grad_NE = (outSmooth(x+1,y-1) - outSmooth(x-1,y+1))*fct; // NE
							   float grad_N = (outSmooth(x,  y-1) - outSmooth(x,  y+1))*0.5; // N
							   float grad_SE = (outSmooth(x+1,y+1) - outSmooth(x-1,y-1))*fct; // SE

							   //***** compute gradient magnitude *********//
							   float grad_mag = grad_E*grad_E + grad_N*grad_N;
							   outGradient(x,y) = grad_mag;

							   //***** compute gradient orientation (continuous version)*******//
							   float angle = 0.0f;
							   if (grad_mag > 0.0f) { angle =  atan2(grad_N, grad_E); }
							   if (angle < 0.0) angle += cimg::PI;
							   outOrientation(x,y) = angle*255.0/cimg::PI + 0.5; // -> outOrientation

							   //***** compute absolute derivations *******//
							   derivative[0](x,y) = grad_E = fabs(grad_E);
							   derivative[1](x,y) = grad_NE = fabs(grad_NE);
							   derivative[2](x,y) = grad_N = fabs(grad_N);
							   derivative[3](x,y) = grad_SE = fabs(grad_SE);

							   //***** compute direction of max derivative //
							   if ((grad_E>grad_NE) && (grad_E>grad_N) && (grad_E>grad_SE)) {
								   dirmax(x,y) = 0; // E
							   } else if ((grad_NE>grad_N) && (grad_NE>grad_SE)){
								   dirmax(x,y) = 1; // NE
							   } else if (grad_N>grad_SE) {
								   dirmax(x,y) = 2; // N
							   } else {
								   dirmax(x,y) = 3; // SE
							   }
							   // one may compute the contiuous dominant direction computation...
							   //outOrientation(x,y) = dirmax(x,y)*255.f/4;  
						   } } // for x,y

					   // directing vectors (E, NE, N, SE)
					   int dir_vector[4][2] = {{1,0}, {1,-1}, {0,-1}, {1,1}};	
					   // direction of max derivative of
					   // current pixel and its two neighbouring pixel (in direction of dir)
					   int dir, dir1, dir2; 

					   //***** thresholding and (canny) non-max-supression *//
					   for (int y = 2; y < ny-2; y++) {
						   for (int x = 2; x < nx-2; x++) {
							   dir = dirmax(x,y);
							   if (derivative[dir](x,y) < threshold) {
								   outThreshold(x,y) = 0.0f;
								   outNMS(x,y) = 0.0f;
							   } else {
								   outThreshold(x,y) = 255.0f;
								   int dx = dir_vector[dir][0];
								   int dy = dir_vector[dir][1];
								   dir1 = dirmax(x + dx,y + dy);
								   dir2 = dirmax(x - dx,y - dy);
								   outNMS(x,y) = 255.f*
									   ((derivative[dir](x,y) > derivative[dir1](x + dx, y + dy)) && 
									   (derivative[dir](x,y) >= derivative[dir2](x-dx,y-dy)));
							   } // -> outThreshold, outNMS
						   } } // for x, y...
} // CannyDiscrete


int main(int argc, char **argv) {
	//***** read command line **********//


	// image after non-max-suppression
	string infile = "Input.bmp";      // required input filename
	string outfileGradient = "Gradient_Pro.bmp";  // saving the (normalised) gradient to file?
	string outfileNMS = "Edge_Output.bmp";       // saving the binary canny edges to file?

	// canny parameters
	float sigma = 1.5f;
	float threshold = 7.0f;




	//***** read image *****************//
	CImg<float> in(infile.c_str());
	const int widthIn = in._width;
	const int heightIn = in._height;
	if ( widthIn == 0 || heightIn == 0 ) {
		cerr << "Error when loading input image." << endl;
		return -1;
	}
	in = in.get_channel(0); // ensure greyscale img!

	//***** declare output images ******//
	CImg<float> outS, outG, outO, outT, outNMS;

	//***** apply Canny filter *********//
	CannyDiscrete(in, sigma, threshold, outS, outG, outO, outT, outNMS);

	//***** display output images ******//
	char  header[100];
	sprintf(header, "gaussian smoothed image: sigma = %f", sigma);
	outS.display(header);
	float maxgrad = 0.0f;
	cimg_forXY(outG,x,y) { maxgrad = std::max(maxgrad, outG(x,y)); }
	std::cout << "normalising [0.." << maxgrad << "] to [0..255]" << std::endl;
	sprintf(header, "gradient magnitude [0..%f]", maxgrad);
	outG.display(header);
	outO.display("orientation map");
	sprintf(header, "thresholded with %f", threshold);
	outT.display(header);
	outNMS.display("non-maximum suppression");


	//***** write output images ********//
	if (outfileGradient.length()>0) { 
		std::cout << "saving gradient to " << outfileGradient << std::endl;
		outG *= (255.f/maxgrad);
		outG.save(outfileGradient.c_str()); 
	}
	if (outfileNMS.length()>0) { 
		std::cout << "saving gradient to " << outfileNMS << std::endl;
		outNMS.save(outfileNMS.c_str()); 
	}

	return 0;
}
