
#include "non_maximum_suppression.h"
#include "overlay.h"

#include <vector>
#include <CImg.h>
using namespace cimg_library;
using namespace std;

void hough (const CImg<float>& img, 
            CImg<float>& HoughSpace, CImg<float>& result,
            float in_thresh, float out_thresh) {
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

int main (int argc, char**argv) {
  
  string infile;
  string outfile;
  // pre-scaling/rotation operations
  float rotate = 0.0f;
  float zoom   = 1.0f;
  // post processing stuf
  float in_thresh = 200.0f; // absolute threshold for gradient magnitude
  float out_thresh = 0.5f;  // relative to global hough space maximum

  //as >> parameter ("i",      infile, "input file name", true);
  //// some pre-processing of the image for experimentation
  //as >> parameter ("rot",    rotate, "pre-rotate the image", false);
  //as >> parameter ("zoom",   zoom, "pre-zoom the image", false);
  //// harris parameters
  //as >> parameter ("in_thresh", in_thresh, 
  //                 "absolut thresh on gradient magnitude", false);
  //as >> parameter ("out_thresh", out_thresh, 
  //                 "relative thresh in hough space", false);
  //// post-processing parameters (non-maximum suppression)
  //as >> parameter ("o",      outfile, "output file name", false);
  //as >> help ();
  //as.defaultErrorHandling();

  // load image and ensure greyscale img!
  CImg<float> input(infile.c_str());
  input = input.get_channel(0);
  if ((rotate != 0.0f) || (zoom!=1)) {
    input.rotate (rotate, input._width/2, input._height/2, zoom);
  }
  input.display ("Input Image");

  // do the transform 
  CImg<float> houghspace, output;
  hough (input, houghspace, output, in_thresh, out_thresh);
  CImgList<float> imgl (input, overlay(input, output), 
                     output, houghspace);
  imgl.display ("input -- combined -- line only -- hough space", 'x');
  if (outfile.size()>0) imgl.get_append('x').save (outfile.c_str());
  return 0;
}
