#include <array>
#include <cmath>
#include "CImg.h"
using namespace cimg_library;

template <class T>
class imageUtil {
private:
    CImg<T> img;

public:
    imageUtil(const char* filepath) {
        this->img.load_jpg(filepath);
    }

    void rgb2grey() {
        CImg<unsigned char> imgTemp(this->img.width(), this->img.height(), 1, 1, 0);
        cimg_forXY(this->img, x, y) {
            imgTemp(x, y, 0) = 0.299 * this->img(x, y, 0) + 0.587 * this->img(x, y, 1) + 0.144 * this->img(x, y, 2);
        };
        this->img = imgTemp;
    }

    std::array<double, 3> getMean() {
        std::array<double, 3> mean = {0, 0, 0};
        cimg_forXY(this->img, x, y) {
            for (int channel = 0; channel < 3; channel ++)
                mean[channel] += this->img(x, y, channel);
        };
        for (int channel = 0; channel < 3; channel ++)
            mean[channel] /= (this->img.width() * this->img.height());
        return mean;
    }

    std::array<double, 3> getVariance() {
        std::array<double, 3> variance = {0, 0, 0};
        std::array<double, 3> mean = this->img.getMean();
        cimg_forXY(this->img, x, y) {
            for (int channel = 0; channel < 3; channel ++)
                variance[channel] += std::pow((this->img(x, y, channel) - mean[channel]), 2);
        };
        for (int channel = 0; channel < 3; channel ++) {
            variance[channel] /= (this->img.width() * this->img.height());
            variance[channel] = std::pow(variance[channel], 1/2);
        }
        return variance;
    }

    CImg<T> colorTransfer(CImg<T> source) {
        CImg<T> target = this->img;

        // convert rgb space to lab space
        source.RGBtoLab();
        target.RGBtoLab();

        // get stats from source image and 
        std::array<double, 3> sourceMean = getMean(source);
        std::array<double, 3> sourceVariance = imageUtil::getVariance(source);
        std::array<double, 3> targetMean = imageUtil::getMean(target);
        std::array<double, 3> targetVariance = imageUtil::getVariance(target);

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
        return target
    }

    std::array<int, 256> getHistogram(int channel = 0) {
        std::array<int, 256> histogram = {0};
        cimg_forXY(this->img, x, y) {
            histogram[img(x, y, channel)] ++;
        };
        return histogram;
    }

    CImg<T>& equalize() {
        for (int channel = 0; channel < this->img.spectrum(); channel ++) {
            std::array<int, 256> histogram = this->img.getHistogram(channel);
            unsigned long cumulated = 0;
            for (int index = 0; index < 256; index ++) {
                cumulated += histogram[index];
                histogram[index] = cumulated;
            }
            cimg_forXY(this->img, x, y) {
                this->img(x, y, channel) = 
                    (T)(255 * histogram[this->img(x, y, channel)/cumulated]);
            };
        }
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
    imageUtil<double> target("../image/Montreal.jpeg");
    CImg<double> transfered = target.colorTransfer(source);
    transfered.display();
}
