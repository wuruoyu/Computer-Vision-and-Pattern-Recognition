#include <iostream>
#include <array>
#include "CImg.h"
using namespace cimg_library;

// define color
class colorSet {
    public:
        static std::array<int, 3> getRed() {
            return {255, 0, 0};
        }

        static std::array<int, 3> getGreen() {
            return {0, 255, 0};
        }

        static std::array<int, 3> getBlue() {
            return {0, 0, 255};
        }

        static std::array<int, 3> getWhite() {
            return {255, 255, 255};
        }

        static std::array<int, 3> getBlack() {
            return {0, 0, 0};
        }

        static bool isWhite(CImg<unsigned char>& image, int x, int y) {
            std::array<int, 3> white = colorSet::getWhite();
            if (image(x, y, 0) == white[0] && image(x, y, 1) == white[1] && image(x, y, 2) == white[2])
                return true;
            return false;
        }

        static bool isBlack(CImg<unsigned char>& image, int x, int y) {
            std::array<int, 3> black = colorSet::getBlack();
            if (image(x, y, 0) == black[0] && image(x, y, 1) == black[1] && image(x, y, 2) == black[2])
                return true;
            return false;
        }

};

class answer {
    public:
        static void questionOne(CImg<unsigned char>& image) {
            image.display();
        }
        
        static void questionTwo(CImg<unsigned char>& image) {
           cimg_forXY(image, x, y) {
               if (colorSet::isWhite(image, x, y)) {
                   std::array<int, 3> red = colorSet::getRed();
                   image(x, y, 0) = red[0];
                   image(x, y, 1) = red[1];
                   image(x, y, 2) = red[2];
               }
               if (colorSet::isBlack(image, x, y)) {
                   std::array<int, 3> green = colorSet::getGreen();
                   image(x, y, 0) = green[0];
                   image(x, y, 1) = green[1];
                   image(x, y, 2) = green[2];
               }               
           };
           image.display();
        }
};

int main(int argc, char* argv[]) {
    CImg<unsigned char> image("1.bmp");
    
    //answer::questionOne(image);
    //answer::questionTwo(image);
    answer::questionThree(image);
}
