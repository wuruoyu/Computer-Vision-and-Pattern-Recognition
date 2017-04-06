#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>
#include "CImg.h"
#include "stdlib.h"
#include "houghDetector.h"
using namespace cimg_library;
using namespace std;


bool operator<(const Point& a, const Point& b) {
    return a.weight > b.weight;
}

class perspectiveTransform {
private:
    std::vector<Point> originalCorner;
    std::vector<Point> changedCorner;
    std::vector< std::pair<Point, Point> > cornerMap;
    std::vector<double> transformMatrix;
    CImg<double> transformed;
    CImg<double> source;

public:
    perspectiveTransform(const std::vector<Point>& _originalCorner, const std::vector<Point>& _changedCorner, CImg<double> _source) {
        this->originalCorner = _originalCorner;
        this->changedCorner = _changedCorner;
        this->source = _source;

        if (originalCorner.size() != changedCorner.size()) {
            throw "points sets size are not equal";
        }
        if (originalCorner.size() != 4) {
            throw "points set size are not 4";
        }

        std::sort(originalCorner.begin(), originalCorner.end(), [](Point& p1, Point& p2)
                                    {return pow(p1.theta, 2) + pow(p1.rho, 2) < pow(p2.theta, 2) + pow(p2.rho, 2);});
        std::sort(changedCorner.begin(), changedCorner.end(), [](Point& p1, Point& p2)
                                    {return pow(p1.theta, 2) + pow(p1.rho, 2) < pow(p2.theta, 2) + pow(p2.rho, 2);});

        for (int i = 0; i < originalCorner.size(); i ++) {
            cornerMap.push_back(make_pair(originalCorner[i], changedCorner[i]));
        }
    }

    // Ax = b
    void getPerspectiveTransformMatrix() {
        CImg<double> A(8, 8, 1, 1);
        CImg<double> b(1, 8, 1, 1);

        // give values to A
        for (int row = 0; row < 4; row ++) {
            A(0, row) = cornerMap[row].second.theta;
            A(1, row) = 0;
            A(2, row) = -cornerMap[row].second.theta * cornerMap[row].first.theta;
            A(3, row) = cornerMap[row].second.rho;
            A(4, row) = 0;
            A(5, row) = -cornerMap[row].second.rho * cornerMap[row].first.theta;
            A(6, row) = 1;
            A(7, row) = 0;
        }
        for (int row = 4; row < 8; row ++) {
            A(0, row) = 0;
            A(1, row) = cornerMap[row - 4].second.theta;
            A(2, row) = -cornerMap[row - 4].second.theta * cornerMap[row - 4].first.rho;
            A(3, row) = 0;
            A(4, row) = cornerMap[row - 4].second.rho;
            A(5, row) = -cornerMap[row - 4].second.rho * cornerMap[row - 4].first.rho;
            A(6, row) = 0;
            A(7, row) = 1;
        }

        // give values to b
        for (int row = 0; row < 4; row ++) {
            b(0, row) = cornerMap[row].first.theta;
        }
        for (int row = 4; row < 8; row ++) {
            b(0, row) = cornerMap[row - 4].first.rho;
        }

        CImg<double> transformMatrixColumn = A.get_invert() * b;

        for (int i = 0; i < 8; i ++) {
            transformMatrix.push_back(transformMatrixColumn(0, i));
        }
        transformMatrix.push_back(1);
    }

    void bilinearInterpolation() {
        transformed = CImg<double>(source.width(), source.height(), 1, 3);
        cimg_forXY(transformed, u, v) {
            double x = (transformMatrix[0] * u + transformMatrix[3] * v + transformMatrix[6]) / (transformMatrix[2] * u + transformMatrix[5] * v + 1);
            double y = (transformMatrix[1] * u + transformMatrix[4] * v + transformMatrix[7]) / (transformMatrix[2] * u + transformMatrix[5] * v + 1);
            int x1 = floor(x);
            int x2 = x1 + 1;
            int y1 = floor(y);
            int y2 = y1 + 1;
            for (int channel = 0; channel < source.spectrum(); channel ++) {
                // if some transformed points fall outside of source image
                if (x1 < 0 || x2 >= source.width())
                    continue;
                if (y1 < 0 || y2 >= source.height())
                    continue;

                transformed(u, v, channel) = (source(x1, y1, channel) * (x2 - x) * (y2 - y) 
                                                + source(x2, y1, channel) * (x - x1) * (y2 - y)
                                                + source(x1, y2, channel) * (x2 - x) * (y - y1)
                                                + source(x2, y2, channel) * (x - x1) * (y - y1)) / ((x2 - x1) * (y2 - y1));
            }
        };   
    }

    void display() {
        this->transformed.display();
    }

    CImg<double> getTransformed() {
        return this->transformed;
    }
};


int main(int argc, char* argv[]) {
    CImg<double> source(argv[1]);    

    // get edge and intersection of the source image
    houghDetector sourceEdge(source);
    sourceEdge.detect();

    // define the mapped corner point
    std::vector<Point> changedPoint;
    changedPoint.push_back(Point(0, 0, -1));
    changedPoint.push_back(Point(sourceEdge.getImg().width(), 0, -1));
    changedPoint.push_back(Point(0, sourceEdge.getImg().height(), -1));
    changedPoint.push_back(Point(sourceEdge.getImg().width(), sourceEdge.getImg().height(), -1));

    try {
        // perform perspective transform
        perspectiveTransform rect(sourceEdge.getIntersection(), changedPoint, source);
        rect.getPerspectiveTransformMatrix();
        rect.bilinearInterpolation();
        (source, rect.getTransformed()).display();
    } catch(const char* msg) {
        cerr << argv[1] << " " << msg << endl;
    }
}
