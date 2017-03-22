//
//  canny.cpp
//  Canny Edge Detector
//
//  Created by Hasan Akgün on 21/03/14.
//  Copyright (c) 2014 Hasan Akgün. All rights reserved.
//

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include "canny.h"
#include "CImg.h"

using namespace std;
using namespace cimg_library;


canny::canny(CImg<unsigned char>& img)
{
	this->img = img;
	
	vector<vector<double>> filter = createFilter(3, 3, 2);

    this->grayscaled = toGrayScale(); //Grayscale the image
    this->gFiltered = useFilter(grayscaled, filter); //Gaussian Filter
    this->sFiltered = sobel(); //Sobel Filter

    non = nonMaxSupp(); //Non-Maxima Suppression
    thres = threshold(non, 20, 40); //Double Threshold and Finalize
	
	//thres.display();
}

CImg<unsigned char> canny::toGrayScale()
{
    grayscaled = CImg<unsigned char>(img.width(), img.height(), 1, 1, 0); //To one channel
    cimg_forXY(img, x, y) {
        grayscaled(x, y, 0) = 0.2126 * img(x, y, 0) + 0.7152 * img(x, y, 1) + 0.0722 * img(x, y, 2);
    };
    return grayscaled;
}

vector<vector<double>> canny::createFilter(int row, int column, double sigmaIn)
{
	vector<vector<double>> filter;

	for (int i = 0; i < row; i++)
	{
        vector<double> col;
        for (int j = 0; j < column; j++)
        {
            col.push_back(-1);
        }
		filter.push_back(col);
	}

	float coordSum = 0;
	float constant = 2.0 * sigmaIn * sigmaIn;

	// Sum is for normalization
	float sum = 0.0;

	for (int x = - row/2; x <= row/2; x++)
	{
		for (int y = -column/2; y <= column/2; y++)
		{
			coordSum = (x*x + y*y);
			filter[x + row/2][y + column/2] = (exp(-(coordSum) / constant)) / (M_PI * constant);
			sum += filter[x + row/2][y + column/2];
		}
	}

	// Normalize the Filter
	for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
            filter[i][j] /= sum;

	return filter;

}

CImg<unsigned char> canny::useFilter(CImg<unsigned char> img_in, vector<vector<double>> filterIn)
{
    int size = (int)filterIn.size()/2;
	CImg<unsigned char> filteredImg = CImg<unsigned char>(img_in.width() - 2*size, img_in.height() - 2*size, 1, 1, 0);
	for (int i = size; i < img_in.width() - size; i++)
	{
		for (int j = size; j < img_in.height() - size; j++)
		{
			double sum = 0;
            
			for (int x = 0; x < filterIn.size(); x++)
				for (int y = 0; y < filterIn.size(); y++)
				{
                    sum += filterIn[x][y] * (double)(img_in(i + x - size, j + y - size, 0));
				}
            
            filteredImg(i-size, j-size, 0) = sum;
		}

	}
	return filteredImg;
}

CImg<unsigned char> canny::sobel()
{
    //Sobel X Filter
    double x1[] = {-1.0, 0, 1.0};
    double x2[] = {-2.0, 0, 2.0};
    double x3[] = {-1.0, 0, 1.0};

    vector<vector<double>> xFilter(3);
    xFilter[0].assign(x1, x1+3);
    xFilter[1].assign(x2, x2+3);
    xFilter[2].assign(x3, x3+3);
    
    //Sobel Y Filter
    double y1[] = {1.0, 2.0, 1.0};
    double y2[] = {0, 0, 0};
    double y3[] = {-1.0, -2.0, -1.0};
    
    vector<vector<double>> yFilter(3);
    yFilter[0].assign(y1, y1+3);
    yFilter[1].assign(y2, y2+3);
    yFilter[2].assign(y3, y3+3);
    
    //Limit Size
    int size = (int)xFilter.size()/2;
    
	CImg<unsigned char> filteredImg(gFiltered.width() - 2*size, gFiltered.height() - 2*size, 1, 1, 0);
    
    angles = CImg<unsigned char>(gFiltered.width() - 2*size, gFiltered.height() - 2*size, 1, 1, 0); //AngleMap

	for (int i = size; i < gFiltered.width() - size; i++)
	{
		for (int j = size; j < gFiltered.height() - size; j++)
		{
			double sumx = 0;
            double sumy = 0;
            
			for (int x = 0; x < xFilter.size(); x++)
				for (int y = 0; y < xFilter.size(); y++)
				{
                    sumx += xFilter[x][y] * (double)(gFiltered(i + x - size, j + y - size, 0)); //Sobel_X Filter Value
                    sumy += yFilter[x][y] * (double)(gFiltered(i + x - size, j + y - size, 0)); //Sobel_Y Filter Value
				}
            double sumxsq = sumx*sumx;
            double sumysq = sumy*sumy;
            
            double sq2 = sqrt(sumxsq + sumysq);
            
            if(sq2 > 255) //Unsigned Char Fix
                sq2 =255;
            filteredImg(i-size, j-size, 0) = sq2;
 
            if(sumx==0) //Arctan Fix
                angles(i-size, j-size, 0) = 90;
            else
                angles(i-size, j-size, 0) = atan(sumy/sumx);
		}
	}
    
    return filteredImg;
}


CImg<unsigned char> canny::nonMaxSupp()
{
    CImg<unsigned char> nonMaxSupped(sFiltered.width()-2, sFiltered.height()-2, 1, 1, 0);
    for (int i=1; i< sFiltered.width() - 1; i++) {
        for (int j=1; j<sFiltered.height() - 1; j++) {
            float Tangent = angles(i, j);

            nonMaxSupped(i-1, j-1) = sFiltered(i,j);
            //Horizontal Edge
            if (((-22.5 < Tangent) && (Tangent <= 22.5)) || ((157.5 < Tangent) && (Tangent <= -157.5)))
            {
                if ((sFiltered(i,j) < sFiltered(i,j+1)) || (sFiltered(i,j) < sFiltered(i,j-1)))
                    nonMaxSupped(i-1, j-1) = 0;
            }
            //Vertical Edge
            if (((-112.5 < Tangent) && (Tangent <= -67.5)) || ((67.5 < Tangent) && (Tangent <= 112.5)))
            {
                if ((sFiltered(i,j) < sFiltered(i+1,j)) || (sFiltered(i,j) < sFiltered(i-1,j)))
                    nonMaxSupped(i-1, j-1) = 0;
            }
            
            //-45 Degree Edge
            if (((-67.5 < Tangent) && (Tangent <= -22.5)) || ((112.5 < Tangent) && (Tangent <= 157.5)))
            {
                if ((sFiltered(i,j) < sFiltered(i-1,j+1)) || (sFiltered(i,j) < sFiltered(i+1,j-1)))
                    nonMaxSupped(i-1, j-1) = 0;
            }
            
            //45 Degree Edge
            if (((-157.5 < Tangent) && (Tangent <= -112.5)) || ((22.5 < Tangent) && (Tangent <= 67.5)))
            {
                if ((sFiltered(i,j) < sFiltered(i+1,j+1)) || (sFiltered(i,j) < sFiltered(i-1,j-1)))
                    nonMaxSupped(i-1, j-1) = 0;
            }
        }
    }
    return nonMaxSupped;
}

CImg<unsigned char> canny::threshold(CImg<unsigned char> imgin,int low, int high)
{
    if(low > 255)
        low = 255;
    if(high > 255)
        high = 255;
    
    CImg<unsigned char> EdgeMat(imgin.width(), imgin.height(), 1, 1, 0);
    
    for (int i=0; i<imgin.width(); i++) 
    {
        for (int j = 0; j<imgin.height(); j++) 
        {
            EdgeMat(i,j, 0) = imgin(i,j,0);
            if(EdgeMat(i,j,0) > high)
                EdgeMat(i,j,0) = 255;
            else if(EdgeMat(i,j,0) < low)
                EdgeMat(i,j,0) = 0;
            else
            {
                bool anyHigh = false;
                bool anyBetween = false;
                for (int x=i-1; x < i+2; x++) 
                {
                    for (int y = j-1; y<j+2; y++) 
                    {
                        if(x <= 0 || y <= 0 || EdgeMat.width() || y > EdgeMat.height()) //Out of bounds
                            continue;
                        else
                        {
                            if(EdgeMat(x,y,0) > high)
                            {
                                EdgeMat(i,j,0) = 255;
                                anyHigh = true;
                                break;
                            }
                            else if(EdgeMat(x,y,0) <= high && EdgeMat(x,y,0) >= low)
                                anyBetween = true;
                        }
                    }
                    if(anyHigh)
                        break;
                }
                if(!anyHigh && anyBetween)
                    for (int x=i-2; x < i+3; x++) 
                    {
                        for (int y = j-1; y<j+3; y++) 
                        {
                            if(x < 0 || y < 0 || x > EdgeMat.width() || y > EdgeMat.height()) //Out of bounds
                                continue;
                            else
                            {
                                if(EdgeMat(x,y,0) > high)
                                {
                                    EdgeMat(i,j,0) = 255;
                                    anyHigh = true;
                                    break;
                                }
                            }
                        }
                        if(anyHigh)
                            break;
                    }
                if(!anyHigh)
                    EdgeMat(i,j,0) = 0;
            }
        }
    }
    return EdgeMat;
}

CImg<unsigned char> canny::getEdge() {
    return this->thres;
}