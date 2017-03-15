#include <array>
#include <cmath>
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

    static std::array<double, 3> getMean(CImg<unsigned char>& img) {
        std::array<double, 3> mean = {0, 0, 0};
        cimg_forXY(img, x, y) {
            for (int channel = 0; channel < 3; channel ++)
                mean[channel] += img(x, y, channel);
        };
        for (int channel = 0; channel < 3; channel ++)
            mean[channel] /= (img.width() * img.height());
        return mean;
    }

    static std::array<double, 3> getVariance(CImg<unsigned char>& img) {
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

    static CImg<unsigned char> colorTransfer(CImg<unsigned char> source, CImg<unsigned char>& target) {
        // convert rgb space to lab space
        source.LabtoRGB();
        target.LabtoRGB();

        // get stats from source image and 
        std::array<double, 3> sourceMean = imageSet::getMean(source);
        std::array<double, 3> sourceVariance = imageSet::getVariance(source);
        std::array<double, 3> targetMean = imageSet::getMean(target);
        std::array<double, 3> targetVariance = imageSet::getVariance(target);

        cimg_forXY(source, x, y) {
            for (int channel = 0; channel < 3; channel ++) {
                // substract mean value from source
                source(x, y, channel) -= sourceMean[channel];
                // scale
                source(x, y, channel) *= (sourceVariance[channel] / targetVariance[channel]);
                // add the mean value from target
                source(x, y, channel) += targetMean[channel];
            }
        };

        source.RGBtoLab();
        return source;
    }
}

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
    CImg<unsigned char> source("../image/vancouver.jpeg");
    CImg<unsigned char> target("../image/vancouver.jpeg");
    CImg<unsigned char> transfered = imageSet::colorTransfer(source, target);
    transfered.display();
}
