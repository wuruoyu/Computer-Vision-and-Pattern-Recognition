#include <iostream>
#include <array>
#include "CImg.h"
using namespace cimg_library;

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

class answer {
    public:
        static void questionOne(CImg<unsigned char>& image) {
            image.display();
        }
        
        static void questionTwo(CImg<unsigned char>& image) {
            cimg_forXY(image, x, y) {
               if (colorSet::isWhite(image, x, y)) {
                   std::array<unsigned char, 3> red = colorSet::getRed();
                   image(x, y, 0) = red[0];
                   image(x, y, 1) = red[1];
                   image(x, y, 2) = red[2];
               }
               if (colorSet::isBlack(image, x, y)) {
                   std::array<unsigned char, 3> green = colorSet::getGreen();
                   image(x, y, 0) = green[0];
                   image(x, y, 1) = green[1];
                   image(x, y, 2) = green[2];
               }               
           };
           image.display();
        }

        static void questionThree(CImg<unsigned char>& image) {
            std::array<unsigned char, 3> blue = colorSet::getBlue();
            unsigned char blueChar[3] = {blue[0], blue[1], blue[2]};
            image.draw_circle(50, 50, 30, blueChar);
            image.display();
        }

        static void questionFour(CImg<unsigned char>& image) {
            std::array<unsigned char, 3> yellow = colorSet::getYellow();
            unsigned char yellowChar[3] = {yellow[0], yellow[1], yellow[2]};
            image.draw_circle(50, 50, 3, yellowChar);
            image.display();
        }

        static void test(CImg<unsigned char>& image) {
            answer::questionOne(image);
            answer::questionTwo(image);
            answer::questionThree(image);
            answer::questionFour(image);
        }
};


int main(int argc, char* argv[]) {
    CImg<unsigned char> image("1.bmp");
    answer::test(image);
}
