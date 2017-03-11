
#include "CImg.h"
using namespace cimg_library;


/** compute Gaussian derivatives filter weights
* \param sigma = bandwidth of the Gaussian 
* \param deriv = computing the 'deriv'-th derivatives of a Gaussian
* the width of the filter is automatically determined from sigma.
* g  = \frac{1}{\sqrt{2\pi}\sigma}   \exp(-0.5 \frac{x^2}{\sigma^2} )
* g' = \frac{x}{\sqrt{2\pi}\sigma^3} \exp(-0.5 \frac{x^2}{\sigma^2} )
*    = -\frac{x}{\sigma^2} g
* g''= (\frac{x^2}{\sigma^2} - 1) \frac{1}{\sigma^2} g
*/
void gauss_filter (CImg<float>& filter, float sigma=1.0f, int deriv=0) {
	float width = 3*sigma;               // may be less width?
	float sigma2 = sigma*sigma;
	filter.assign(int(2*width)+1);

	int i=0;
	for (float x=-width; x<=width; x+=1.0f) {
		float g = exp(-0.5*x*x/sigma2) / sqrt(2*cimg::PI) / sigma;
		if (deriv==1) g *= -x/sigma2;
		if (deriv==2) g *= (x*x/sigma2 - 1.0f)/sigma2;
		filter[i] = g ;
		//printf ("i=%f -> %f\n", x, filter[i]);
		i++;
	}
}
