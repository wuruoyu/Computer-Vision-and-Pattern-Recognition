#include <array>
#include <cmath>
#include <iostream>
#include "CImg.h"
using namespace cimg_library;


class imageSet {
public:
    static void rgb2grey(CImg<unsigned char>& img) {
        CImg<unsigned char> imgTemp(img.width(), img.height(), 1, 1, 0);
        cimg_forXY(img, x, y) {
            imgTemp(x, y, 0) = 0.299 * img(x, y, 0) + 0.587 * img(x, y, 1) + 0.144 * img(x, y, 2);
        };
        img = imgTemp;
    }

    static std::array<double, 3> getMean(CImg<double>& img) {
        std::array<double, 3> mean = {0, 0, 0};
        cimg_forXY(img, x, y) {
            for (int channel = 0; channel < 3; channel ++)
                mean[channel] += img(x, y, channel);
        };
        for (int channel = 0; channel < 3; channel ++)
            mean[channel] /= (img.width() * img.height());
        return mean;
    }

    static std::array<double, 3> getVariance(CImg<double>& img) {
        std::array<double, 3> variance = {0, 0, 0};
        std::array<double, 3> mean = imageSet::getMean(img);
        cimg_forXY(img, x, y) {
            for (int channel = 0; channel < 3; channel ++)
                variance[channel] += std::pow((img(x, y, channel) - mean[channel]), 2);
        };
        for (int channel = 0; channel < 3; channel ++) {
            variance[channel] /= (img.width() * img.height());
            variance[channel] = std::pow(variance[channel], 1/2);
        }
        return variance;
    }

    static CImg<unsigned char> colorTransfer(CImg<double>& source, CImg<double> target) {
        // convert rgb space to lab space
        source.RGBtoLab();
        target.RGBtoLab();

        // get stats from source image and 
        std::array<double, 3> sourceMean = imageSet::getMean(source);
        std::array<double, 3> sourceVariance = imageSet::getVariance(source);
        std::array<double, 3> targetMean = imageSet::getMean(target);
        std::array<double, 3> targetVariance = imageSet::getVariance(target);

        cimg_forXY(target, x, y) {
            for (int channel = 0; channel < 3; channel ++) {
                // substract mean value from source
                target(x, y, channel) -= targetMean[channel];
                // scale
                target(x, y, channel) *= (targetVariance[channel] / sourceVariance[channel]);
                // add the mean value from target
                target(x, y, channel) += sourceMean[channel];
            }
        };

        target.LabtoRGB();
        return target;
    }

    static std::array<int, 256> getHistogram(CImg<double>& img ,int channel = 0) {
        std::array<int, 256> histogram = {0};
        cimg_forXY(img, x, y) {
            histogram[img(x, y, channel)] ++;
        };
        return histogram;
    }

    static CImg<double>& equalize(CImg<double>& img) {
        for (int channel = 0; channel < img.spectrum(); channel ++) {
            std::array<int, 256> histogram = getHistogram(img, channel);
            unsigned long cumulated = 0;
            for (int index = 0; index < 256; index ++) {
                cumulated += histogram[index];
                histogram[index] = cumulated;
            }
            cimg_forXY(img, x, y) {
                img(x, y, channel) = (255 * histogram[img(x, y, channel)]/(double)cumulated);
            };
        }
        return img;
    }
};

// color interface
class colorSet {
public:
    static std::array<unsigned char, 3> getRed() {
        return {255, 0, 0};
    }

    static std::array<unsigned char, 3> getGreen() {
        return {0, 255, 0};
    }

    static std::array<unsigned char, 3> getBlue() {
        return {0, 0, 255};
    }

    static std::array<unsigned char, 3> getWhite() {
        return {255, 255, 255};
    }

    static std::array<unsigned char, 3> getBlack() {
        return {0, 0, 0};
    }

    static std::array<unsigned char, 3> getYellow() {
        return {255, 255, 0};
    }

    static bool isWhite(CImg<unsigned char>& image, int x, int y) {
        std::array<unsigned char, 3> white = colorSet::getWhite();
        if (image(x, y, 0) == white[0] && image(x, y, 1) == white[1] && image(x, y, 2) == white[2])
            return true;
        return false;
    }

    static bool isBlack(CImg<unsigned char>& image, int x, int y) {
        std::array<unsigned char, 3> black = colorSet::getBlack();
        if (image(x, y, 0) == black[0] && image(x, y, 1) == black[1] && image(x, y, 2) == black[2])
            return true;
        return false;
    }
};


int main(int argc, char* argv[]) {
    CImg<double> source("../image/lena.jpg");
    imageSet::equalize(source).display();

    /*
    CImg<double> target("../image/Montreal.jpeg");
    CImg<double> transfered = imageSet::colorTransfer(source, target);
    transfered.display();
    */
}
