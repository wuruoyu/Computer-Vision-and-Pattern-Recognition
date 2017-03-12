#include "canny.h"
using namespace std;

int main(int argc, char* argv[]) {
	canny img(argv[1]);
	img.setLowThreshold(atof(argv[2]));
	img.setHighThreshold(atof(argv[3]));
	img.setGaussianKernelWidth(atof(argv[4]));
	img.setGaussianKernelRadius(atof(argv[5]));
	img.setContrastNormalised(atoi(argv[6]));
	img.run();
	return 0;
}