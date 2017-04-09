#include <array>
#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include "CImg.h"
#include "canny.h"
#include "stdlib.h"
using namespace cimg_library;

struct Point
{
    int theta;
    int rho;
    int weight;

    Point(int _x, int _y, int _weigt)
        :theta(_x), rho(_y), weight(_weigt){}
};

struct Line
{
    double a;
    double b;
    double c;

    Line(double _a, double _b, double _c)
        :a(_a), b(_b), c(_c){}
};

bool operator<(const Point& a, const Point& b) {
        return a.weight > b.weight;
}

class houghDetector {

private:
    CImg<unsigned char> edgeImg;
    CImg<double> houghSpace;
    CImg<unsigned char> img;
    std::priority_queue<Point> max;
    int rhoGridWidth;
    int thetaFilterWidth;
    int rhoFilterWidth;
    double verticalThreshold;

public:
    houghDetector(canny& cannySource) {
        img = cannySource.getImg();
        edgeImg = cannySource.getEdge();
        rhoGridWidth = std::pow((std::pow(edgeImg.width(), 2) + std::pow(edgeImg.height(), 2)), 0.5);
        houghSpace = CImg<double>(360, rhoGridWidth, 1, 1, 0);
        thetaFilterWidth = 30;
        rhoFilterWidth = 450;
        verticalThreshold = 0.1;
    }

    CImg<unsigned char> getEdgeImg() {
        return this->edgeImg;
    }

    void detect() {
        gray2Hough();
        findMax();
        hough2gray();
    }

    void gray2Hough() {
        cimg_forXY(edgeImg, x, y) {
            if (edgeImg(x, y) == 1) {
                cimg_forX(houghSpace, theta) {
                    double rTheta = (double)theta * 3.14 / 180.0;
                    int rho = (int)(x * cos(rTheta) + y * sin(rTheta));
                    if (rho >= 0 && rho < houghSpace.height())
                        houghSpace(theta, rho) ++;
                };
            }
        };
    }

    bool inRegion(int theta, int rho) {
        double a = cos(theta * 3.14 / 180.0);
        double b = sin(theta * 3.14 / 180.0);
        double c = -rho;

        const int minY = 0;
        const int maxY = img.height() - 1;
        const int minX = 0;
        const int maxX = img.width() - 1;

        if (b == 0){
            if (-c/a > minX && -c/a < maxX)
                return true;
            return false;
        }
        else {
            double k = -a/b;
            c = -c/b;
            if (c > minY && c < maxY)
                return true;
            if (k * maxX + c > minY && k * maxX + c < minY)
                return true;
            if (-c/k > minX && -c/k < maxX)
                return true;
            if (maxY/k - c/k > minX && maxY/k - c/k < maxX)
                return true;
            return false;
        }
    }

    bool isLocalMax(int theta, int rho) {
        for (int thetaShift = -thetaFilterWidth; thetaShift <= thetaFilterWidth; thetaShift ++) {
            for (int rhoShift = -rhoFilterWidth; rhoShift <= rhoFilterWidth; rhoShift ++) {
                int t = theta + thetaShift;
                int r = rho + rhoShift;
                if (thetaShift == 0 && rhoShift == 0)
                    continue;
                if (r < 0 || r >= houghSpace.height())
                    continue;
                if (t < 0)
                    t += 360;
                if (t >= 360)
                    t %= 360;
                if (houghSpace(t, r) > houghSpace(theta, rho))
                    return false;
            }
        }
        return true;
    }

    void findMax() {
        while (max.size() < 4) {
            max.push(Point(-1, -1, 0));
        }

        cimg_forXY(houghSpace, theta, rho) {
            if (houghSpace(theta, rho) > max.top().weight) {
                if (isLocalMax(theta, rho)) {
                    if (inRegion(theta, rho)) {
                        max.pop();
                        max.push(Point(theta, rho, houghSpace(theta, rho)));
                    }
                }
            }
        };
    }

    void drawLine(CImg<double>& result, std::vector<Line>& line) {
        for (int i = 0; i < line.size(); i ++) {
            const int minY = 0;
            const int maxY = result.height() - 1;
            const int minX = 0;
            const int maxX = result.width() - 1;

            const int color[] = {0, 255, 255};

            std::vector<std::pair<int, int>> crossPoint;

            if (line[i].a == 0) {
                crossPoint.push_back(std::pair<int, int>(minX, -line[i].c / line[i].b));
                crossPoint.push_back(std::pair<int, int>(maxX, -line[i].c / line[i].b));
            }
            else if (line[i].b == 0) {
                crossPoint.push_back(std::pair<int, int>(-line[i].c / line[i].a, minY));
                crossPoint.push_back(std::pair<int, int>(-line[i].c / line[i].a, maxY));
            }
            else {
                if (- line[i].c / line[i].b > minY && - line[i].c / line[i].b < maxY)
                    crossPoint.push_back(std::pair<int, int>(0, - line[i].c / line[i].b));
                if (- line[i].a / line[i].b * maxX - line[i].c / line[i].b > minY && - line[i].a / line[i].b * maxX - line[i].c / line[i].b < maxY)
                    crossPoint.push_back(std::pair<int, int>(maxX, - line[i].a / line[i].b * maxX - line[i].c / line[i].b));
                if (- line[i].c / line[i].a > minX && - line[i].c / line[i].a < maxX)
                    crossPoint.push_back(std::pair<int, int>(- line[i].c / line[i].a, 0));
                if (- line[i].b / line[i].a * maxY - line[i].c / line[i].a > minX && - line[i].b / line[i].a * maxY - line[i].c / line[i].a < maxX)
                    crossPoint.push_back(std::pair<int, int>(- line[i].b / line[i].a * maxY - line[i].c / line[i].a, maxY));
            }
            result.draw_line(crossPoint[0].first, crossPoint[0].second, crossPoint[1].first, crossPoint[1].second, color);
        }
    }

    bool ifVertical(Line& line1, Line& line2) {
        if (abs(line1.a * line2.a + line1.b * line2.b) < verticalThreshold) {
            return true;
        }
        return false;
    }

    void drawIntersection(CImg<double>& result, std::vector<Line> line) {
        const int color[] = {0, 255, 0};

        const int minY = 0;
        const int maxY = result.height() - 1;
        const int minX = 0;
        const int maxX = result.width() - 1;

        for (int i = 0; i < line.size(); i ++) {
            for (int j = i + 1; j < line.size(); j ++) {
                if (ifVertical(line[i], line[j])) {
                    int x = (-line[j].b * line[i].c + line[i].b * line[j].c) / (-line[j].a * line[i].b + line[i].a * line[j].b);
                    int y = 0;
                    if (line[i].b != 0)
                        y = - line[i].a / line[i].b * x - line[i].c / line[i].b;
                    else
                        y = -line[j].a / line[j].b * x - line[j].c / line[j].b;
                    if (x > minX && x < maxX && y > minY && y < maxY) {
                        result.draw_circle(x, y, 50, color);
                    }
                }
            }
        }
    }

    void hough2gray() {
        std::vector<Line> line;
        CImg<double> result(this->img);

        while (!max.empty()) {
            line.push_back(Line(cos(max.top().theta * 3.14 / 180.0), sin(max.top().theta * 3.14 / 180.0), -max.top().rho));
            max.pop();
        }

        // print out the line's coefficients
        for (int i = 0; i < line.size(); i ++) {
            cout << line[i].a << "x + " << line[i].b << "y + " << line[i].c << " = 0" << endl;
        }

        drawLine(result, line);

        drawIntersection(result, line);

        result.display();
    }

};


int main(int argc, char* argv[]) {
    CImg<double> source(argv[1]);    

    canny cannySource(source);
    cannySource.setLowThreshold(atoi(argv[2]));
    cannySource.setHighThreshold(atoi(argv[3]));
    cannySource.setGaussianKernelWidth(atoi(argv[4]));
    cannySource.setGaussianKernelRadius(atoi(argv[5]));

    cannySource.run();
    cannySource.showEdgeDetected();

    houghDetector sourceEdge(cannySource);
    sourceEdge.detect();
}
