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
    int verticalThreshold;

public:
    houghDetector(canny& cannySource) {
        img = cannySource.getImg();
        this->edgeImg = cannySource.getEdge();
        rhoGridWidth = std::pow((std::pow(edgeImg.width(), 2) + std::pow(edgeImg.height(), 2)), 0.5);
        houghSpace = CImg<double>(360, rhoGridWidth, 1, 1, 0);
        thetaFilterWidth = 10;
        rhoFilterWidth = 150;
        verticalThreshold = 0.05;
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
            if (edgeImg(x, y) == 255) {
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
                if (thetaShift == 0 && rhoShift == 0)
                    continue;
                if (theta + thetaShift < 0 || theta + thetaShift >= houghSpace.width())
                    continue;
                if (rho + rhoShift < 0 || rho + rhoShift >= houghSpace.height())
                    continue;
                if (houghSpace(theta + thetaShift, rho + rhoShift) > houghSpace(theta, rho))
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

            const int color[] = {0, 0, 255};

            std::vector<std::pair<double, double>> crossPoint;

            if (line[i].a == 0) {
                crossPoint.push_back(std::pair<double, double>(minX, -line[i].c / line[i].b));
                crossPoint.push_back(std::pair<double, double>(maxX, -line[i].c / line[i].b));
            }
            else if (line[i].b == 0) {
                crossPoint.push_back(std::pair<double, double>(-line[i].c / line[i].a, minY));
                crossPoint.push_back(std::pair<double, double>(-line[i].c / line[i].a, maxY));
            }
            else {
                if (- line[i].c / line[i].b > minY && - line[i].c / line[i].b < maxY)
                    crossPoint.push_back(std::pair<double, double>(0, - line[i].c / line[i].b));
                if (- line[i].a / line[i].b * maxX - line[i].c / line[i].b > minY && - line[i].a / line[i].b * maxX - line[i].c / line[i].b < maxY)
                    crossPoint.push_back(std::pair<double, double>(maxX, - line[i].a / line[i].b * maxX - line[i].c / line[i].b));
                if (- line[i].c / line[i].a > minX && - line[i].c / line[i].a < maxX)
                    crossPoint.push_back(std::pair<double, double>(- line[i].c / line[i].a, 0));
                if (- line[i].b / line[i].a * maxY - line[i].c / line[i].a > minX && - line[i].b / line[i].a * maxY - line[i].c / line[i].a < maxX)
                    crossPoint.push_back(std::pair<double, double>(- line[i].b / line[i].a * maxY - line[i].c / line[i].a, maxY));
            }

            result.draw_line(crossPoint[0].first, crossPoint[0].second, crossPoint[1].first, crossPoint[1].second, color);
        }
    }

    bool ifVertical(Line& line1, Line& line2) {
        if (line1.a== 0) {
            if (line2.b == 0)
                return true;
            return false;
        }
        else if (line2.a == 0) {
            if (line1.b == 0)
                return true;
            return false;
        }
        else {
            if (abs(line1.a * line1.b + line2.a * line2.b) < verticalThreshold)
                return true;
            return false;
        }
    }

    void drawIntersection(CImg<double>& result, std::vector<Line> line) {
        const int color[] = {0, 255, 0};

        const int minY = 0;
        const int maxY = result.height() - 1;
        const int minX = 0;
        const int maxX = result.width() - 1;

        for (int i = 0; i < line.size(); i ++) {
            for (int j = i; j < line.size(); j ++) {
                if (ifVertical(line[i], line[j])) {
                    int y = (line[2].a * line[1].c - line[1].a * line[2].c) / (line[1].a * line[2].b - line[2].a * line[1].b);
                    int x = -(line[1].b / line[1].a) * y - line[1].c / line[1].a;
                    if (x > minX && x < maxX && y > minY && y < maxY)
                        result.draw_circle(x, y, 5, color);
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
    int low = atoi(argv[2]);
    int high = atoi(argv[3]);

    canny cannySource(source, low, high);
    houghDetector sourceEdge(cannySource);
    sourceEdge.detect();
}
