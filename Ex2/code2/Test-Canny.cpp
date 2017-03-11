// Test-Canny.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"


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
				   CImg<float> &outNMS) 
{
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


#include "non_maximum_suppression.h"
void hough (const CImg<float>& img, 
			CImg<float>& HoughSpace, CImg<float>& result,
			float in_thresh, float out_thresh) 
{
	const int WIDTH = img._width;
	const int HEIGHT = img._height;
	result.assign(WIDTH,HEIGHT);
	result.fill(0.0f);

	// Hough transform *******************************************************
	const float WIDTH2 = 0.5*WIDTH;
	const float HEIGHT2 = 0.5*HEIGHT;
	const float DIAGONAL = sqrt(WIDTH*WIDTH+HEIGHT*HEIGHT);
	const int OFFSET_N = (int)DIAGONAL;             // how many bins?
	const int THETA_N = 100;                        // how many bins?
	const float THETA_STEP = cimg::PI/(float)THETA_N;
	HoughSpace.assign (THETA_N, OFFSET_N);
	HoughSpace.fill (0.0f);

	cimg_forXY(img, x, y) {
		if (img(x,y)<in_thresh) continue;
		// cast the vote
		for (int t=0; t<THETA_N; t++) {
			float theta = (float)(t*THETA_STEP);
			float offset=(x-WIDTH2)*sin(theta)+(y-HEIGHT2)*cos(theta);
			int offset_int = (int)(OFFSET_N*(offset/DIAGONAL+0.5));
			if (offset_int<0) {
				printf ("%d %d %d %d\n", x, y, t, offset_int);
			}
			HoughSpace(t,offset_int)++;                  // hard voting
			//HoughSpace(t,offset_int) += img(x,y);      // soft voting
		}
	}
	printf ("Voting done\n");
	HoughSpace.display ("HoughSpace");

	// identify the lines *****************************************************
	float maxvote = HoughSpace(0,0);
	for (int i=0; i<THETA_N*OFFSET_N; i++) maxvote = max(maxvote, HoughSpace[i]);
	TVectorOfPairs nonmax;
	non_maximum_suppression(HoughSpace, nonmax, out_thresh*maxvote, 4);
	printf ("Suppression done: %d lines found\n", nonmax.size());

	// draw the lines *********************************************************
	for (unsigned i=0; i<nonmax.size(); i++) {
		float theta  = THETA_STEP*nonmax[i].first;
		float offset = DIAGONAL*(float(nonmax[i].second)/OFFSET_N-0.5f);
		printf ("line: theta=%f offset=%f strength=%f\n", 
			theta, offset, HoughSpace(nonmax[i].first, nonmax[i].second));

		/* draw line : two cases.. */
		if ((theta < cimg::PI/4) || (theta>(3*cimg::PI/4.0f))) {
			// solving    offset=(x-WIDTH)*sin(theta)+(y-HEIGHT2)*cos(theta)
			// for   y(x) = (offset-(x-WIDTH)*sin(theta)) / cos(theta) + HEIGHT2
			for (int x=0; x<WIDTH; x++) {
				int y = (int)((offset-(x-WIDTH2)*sin(theta)) / cos(theta) + HEIGHT2);
				if ((y<0) || (y>HEIGHT)) continue; /* line outside of image */
				result(x,y) = 255.0f;
			} 
		} else {
			for (int y=0; y<HEIGHT; y++) {
				int x = (int)((offset-(y-HEIGHT2)*cos(theta)) / sin(theta) + WIDTH2);
				if ((x<0) || (x>WIDTH)) continue; /* line outside of image */
				result(x,y) = 255.0f;
			}
		}
	} /* draw all lines */

	// rescale hough space for final drawing
	HoughSpace *= 255/maxvote;
} /* hough */



int TestCanny()
{
	// image after non-max-suppression
	string infile = "Input.bmp";      // required input filename
	string outfileGradient = "Gradient_Pro.bmp";  // saving the (normalised) gradient to file?
	string outfileNMS = "Edge_Output.bmp";       // saving the binary canny edges to file?

	// canny parameters
	float sigma = 1.5f;
	float threshold = 6.0f;

	//***** read image *****************//
	CImg<float> inColor(infile.c_str());
	CImg<float> in = inColor.get_RGBtoGray(); // ensure greyscale img!
	const int widthIn = in._width;
	const int heightIn = in._height;
	if ( widthIn == 0 || heightIn == 0 ) {
		cerr << "Error when loading input image." << endl;
		return -1;
	}

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




void TestHough()
{

	string infile="Edge_Output.bmp";
	string outfile="Hough_Output.bmp";
	// pre-scaling/rotation operations
	float rotate = 0.0f;
	float zoom   = 1.0f;
	// post processing stuf
	float in_thresh = 200.0f; // absolute threshold for gradient magnitude
	float out_thresh = 0.5f;  // relative to global hough space maximum


	// load image and ensure greyscale img!
	CImg<float> input(infile.c_str());
	input = input.get_channel(0);


	input.display ("Input Image");

	// do the transform 
	CImg<float> houghspace, output;
	hough (input, houghspace, output, in_thresh, out_thresh);

	output.display();

	if (outfile.size()>0) output.save (outfile.c_str());
}



int _tmain(int argc, _TCHAR* argv[])
{
	TestCanny();
	TestHough();
	return 0;
}



