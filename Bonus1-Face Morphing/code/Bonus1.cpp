#include <cmath>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>
#include "CImg.h"
#include "stdlib.h"
using namespace cimg_library;
using namespace std;

struct Point
{
    int theta;
    int rho;
    int weight;

    Point(int _x, int _y, int _weigt)
        :theta(_x), rho(_y), weight(_weigt){}
};

bool operator<(const Point& a, const Point& b) {
    return a.weight > b.weight;
}

struct Triangle
{
    vector<Point> points;

    Triangle(Point a, Point b, Point c) {
        points.push_back(a);
        points.push_back(b);
        points.push_back(c);
    }
};

// return the PI value
constexpr double pi() {
    return atan(1) * 4;
}

class Delaunay {
public:
    static bool isPointsOnOneLine(const Point& a, const Point& b, const Point& c) {
        double k1 = (a.rho - b.rho) / double(a.theta - b.theta);
        double k2 = (b.rho - c.rho) / double(b.theta - c.theta);
        if (abs(k1 - k2) < 0.000001) 
            return true;
        return false;
    }

    // judge if m and p2 are on the same side of line p0p1
    static bool isOnSameSide(const Point& p0, const Point& p1, const Point& p2, const Point& m) {
        // line p0p1 is vertical
        if (abs((p0.rho - p1.rho) / double(p0.theta - p1.theta)) < 0.000001) {
            if ((p2.theta - p0.theta) * (m.theta - p0.theta) > 0)
                return true;
            else
                return false;
        }
        else {
            // calculate the efficient of line ab
            double k = (p0.rho - p1.rho) / double(p0.theta - p1.theta);
            double b = p0.rho - k * p0.theta;
            if ((k * p2.theta + b - p2.rho) * (k * m.theta + b - m.rho) < 0)
                return false;
            return true;
        }
    }

    // calculate the magnitude of a vector
    static double getMagnitude(int vec[2]) {
        return sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
    }

    // calculate the angle abc
    static double angle(const Point& a, const Point& b, const Point& c) {
        int ba[2] = {b.theta - a.theta, b.rho - a.rho};
        int bc[2] = {b.theta - c.theta, b.rho - c.rho};

        return acos((inner_product(begin(ba), end(ba), begin(bc), 0.0)) / (getMagnitude(ba) * getMagnitude(bc)));
    }

    // examine if point m in the circumscribed circle formed by (a, b, c)
    // the explanation of algorithm is in the testing report
    static bool isInCircumscribedCircle(const Point& m, const Point& a, const Point& b, const Point& c) {
        // if m is on line ab
        if (isPointsOnOneLine(a, b, m))
            return true;
        
        // calculate the angle
        double acb = angle(a, c, b);
        double amb = angle(a, m, b);

        if (isOnSameSide(a, b, c, m)) {
            if (acb <= amb)
                return true;
            return false;
        }
        else {
            if (acb + amb >= pi())
                return true;
            return false;
        }
    }

    static vector<Triangle> triangulate(const vector<Point>& points) {
        if (points.size() < 3) {
            cerr << "too few points" << endl;
            exit(-1);
        }

        vector<Triangle> triangles;

        // for all the combination of points
        for (int i = 0; i < points.size() - 2; i ++) {
            for (int j = i + 1; j < points.size() - 1; j ++) {
                for (int k = j + 1; k < points.size(); k ++) {
                    // if a, b, c are on the same line
                    if (isPointsOnOneLine(points[i], points[j], points[k]))
                        continue;

                    // no other points are in circumscribed circle
                    bool delaunayFlag = true;
                    for (int h = 0; h < points.size(); h ++) {
                        if (h == i || h == j || h == k)
                            continue;
                        if (isInCircumscribedCircle(points[h], points[i], points[j], points[k])) {
                            delaunayFlag = false;  
                            break;
                        }
                    }
                    if (delaunayFlag) {
                        triangles.push_back(Triangle(points[i], points[j], points[k]));
                    }
                }
            }
        }
        return triangles;
    }
};


class imageMorphing {
private:
    CImg<double> source;
    CImg<double> dest;
    // we use the point's weight as index
    vector<Point> sourceControlPoint;
    vector<Point> destControlPoint;

    vector<Triangle> sourceTriangle;
    vector<Triangle> destTriangle;

    int sourcePointColor[3] = {0, 0, 255};
    int destPointColor[3] = {255, 0, 0};

    int interval;

public:
    imageMorphing(CImg<double> _source, CImg<double> _dest, int _interval = 20) {
        this->source = _source;
        this->dest = _dest;
        this->interval = _interval;

        // add the corner points into the vector
        sourceControlPoint.push_back(Point(0, 0, sourceControlPoint.size()));
        sourceControlPoint.push_back(Point(0, source.height(), sourceControlPoint.size()));
        sourceControlPoint.push_back(Point(source.width(), 0, sourceControlPoint.size()));
        sourceControlPoint.push_back(Point(source.width(), source.height(), sourceControlPoint.size()));

        destControlPoint.push_back(Point(0, 0, destControlPoint.size()));
        destControlPoint.push_back(Point(0, dest.height(), destControlPoint.size()));
        destControlPoint.push_back(Point(dest.width(), 0, destControlPoint.size()));
        destControlPoint.push_back(Point(dest.width(), dest.height(), destControlPoint.size()));
    }

    void getControlPoint() {
        CImg<double> sourceCopy = source;
        CImg<double> destCopy = dest;

        CImgDisplay sDisp(sourceCopy, "Source");
        CImgDisplay dDisp(destCopy, "Dest");

        while (!sDisp.is_closed()) {
            sDisp.wait();
            if ((sDisp.button() & 1) && sDisp.mouse_x() != -1 && sDisp.mouse_y() != -1) {
                sourceControlPoint.push_back(Point(sDisp.mouse_x(), sDisp.mouse_y(), sourceControlPoint.size()));
                sourceCopy.draw_circle(sourceControlPoint.back().theta, sourceControlPoint.back().rho, 5, sourcePointColor);
                sourceCopy.display(sDisp);
            }
        }

        while (!dDisp.is_closed()) {
            dDisp.wait();
            if ((dDisp.button() & 1) && dDisp.mouse_x() != -1 && dDisp.mouse_y() != -1) {
                destControlPoint.push_back(Point(dDisp.mouse_x(), dDisp.mouse_y(), destControlPoint.size()));
                destCopy.draw_circle(destControlPoint.back().theta, destControlPoint.back().rho, 5, destPointColor);
                destCopy.display(dDisp);
            }
        }

        if (sourceControlPoint.size() != destControlPoint.size()) {
            cerr << "the number of control points are not equal" << endl;
            exit(-1);
        }
    }

    void drawTriangle() {
        CImg<double> sourceCopy = source;
        CImg<double> destCopy = dest;

        for (int i = 0; i < sourceTriangle.size(); i ++) {
            sourceCopy.draw_line(sourceTriangle[i].points[0].theta, sourceTriangle[i].points[0].rho, 
                            sourceTriangle[i].points[1].theta, sourceTriangle[i].points[1].rho, sourcePointColor);
            sourceCopy.draw_line(sourceTriangle[i].points[0].theta, sourceTriangle[i].points[0].rho,
                            sourceTriangle[i].points[2].theta, sourceTriangle[i].points[2].rho, sourcePointColor);
            sourceCopy.draw_line(sourceTriangle[i].points[1].theta, sourceTriangle[i].points[1].rho,
                            sourceTriangle[i].points[2].theta, sourceTriangle[i].points[2].rho, sourcePointColor);
        }

        for (int i = 0; i < destTriangle.size(); i ++) {
            destCopy.draw_line(destTriangle[i].points[0].theta, destTriangle[i].points[0].rho, 
                            destTriangle[i].points[1].theta, destTriangle[i].points[1].rho, destPointColor);
            destCopy.draw_line(destTriangle[i].points[0].theta, destTriangle[i].points[0].rho,
                            destTriangle[i].points[2].theta, destTriangle[i].points[2].rho, destPointColor);
            destCopy.draw_line(destTriangle[i].points[1].theta, destTriangle[i].points[1].rho,
                            destTriangle[i].points[2].theta, destTriangle[i].points[2].rho, destPointColor);
        }

        (sourceCopy, destCopy).display();
    }

    void syncTriangle() {
        for (int i = 0; i < sourceTriangle.size(); i ++) {
            destTriangle.push_back(Triangle(destControlPoint[sourceTriangle[i].points[0].weight],
                                            destControlPoint[sourceTriangle[i].points[1].weight],
                                            destControlPoint[sourceTriangle[i].points[2].weight]));
        }
    }

    vector<Triangle> getMiddleTriangle(double ratio) {
        vector<Triangle> middleTriangle;
        for (int i = 0; i < sourceTriangle.size(); i ++) {
            int ax = sourceTriangle[i].points[0].theta * ratio + destTriangle[i].points[0].theta * (1 - ratio);
            int ay = sourceTriangle[i].points[0].rho * ratio + destTriangle[i].points[0].rho * (1 - ratio);
            int bx = sourceTriangle[i].points[1].theta * ratio + destTriangle[i].points[1].theta * (1 - ratio);
            int by = sourceTriangle[i].points[1].rho * ratio + destTriangle[i].points[1].rho * (1 - ratio);
            int cx = sourceTriangle[i].points[2].theta * ratio + destTriangle[i].points[2].theta * (1 - ratio);
            int cy = sourceTriangle[i].points[2].rho * ratio + destTriangle[i].points[2].rho * (1 - ratio);

            middleTriangle.push_back(Triangle(Point(ax, ay, -1),
                                              Point(bx, by, -1),
                                              Point(cx, cy, -1)));
        }
        return middleTriangle;
    }

    vector<double> getTransformMatrix(Triangle begin, Triangle end) {
        // Ax = b
        CImg<double> A(6, 6, 1, 1);
        CImg<double> b(1, 6, 1, 1);

        // give values to A
        for (int row = 0; row < 3; row ++) {
            A(0, row) = begin.points[row].theta;
            A(1, row) = begin.points[row].rho;
            A(2, row) = 1;
            A(3, row) = 0;
            A(4, row) = 0;
            A(5, row) = 0;
        }
        for (int row = 3; row < 6; row ++) {
            A(0, row) = 0;
            A(1, row) = 0;
            A(2, row) = 0;
            A(3, row) = begin.points[row - 3].theta;
            A(4, row) = begin.points[row - 3].rho;
            A(5, row) = 1;
        }

        // give value to b
        for (int row = 0; row < 3; row ++) {
            b(0, row) = end.points[row].theta;
        }
        for (int row = 3; row < 6; row ++) {
            b(0, row) = end.points[row - 3].rho;
        }

        CImg<double> tranformMatrixCImg = A.get_invert() * b;

        vector<double> tranformMatrix;
        for (int i = 0; i < 6; i ++) {
            tranformMatrix.push_back(tranformMatrixCImg(0, i));
        }

        return tranformMatrix;
    } 

    double getArea(Point a, Point b, Point c) {
        double ab = sqrt(pow((a.theta - b.theta), 2) + pow((a.rho - b.rho), 2));
        double bc = sqrt(pow((b.theta - c.theta), 2) + pow((b.rho - c.rho), 2));
        double ac = sqrt(pow((a.theta - c.theta), 2) + pow((a.rho - c.rho), 2));
        double s = (ab + bc + ac) / 2;
        return sqrt(s * (s - ab) * (s - bc) * (s - ac));
    }

    bool isPointInTriangle(Point x, Triangle tri) {
        double sum = getArea(x, tri.points[0], tri.points[1]) + getArea(x, tri.points[0], tri.points[2]) + getArea(x, tri.points[1], tri.points[2]);
        double triArea = getArea(tri.points[0], tri.points[1], tri.points[2]);
        if (abs(sum - triArea) <= 0.00001)
            return true;
        return false;
    }

    void morph() {
        // # of interval frames
        for (int i = interval; i >= 0; i --) {
            double ratio = i / double(interval);

            vector<Triangle> middleTriangle = getMiddleTriangle(ratio);
            vector< vector<double> > middle2source;
            vector< vector<double> > middle2dest; 

            // get transform matrix for all triangles
            for (int j = 0; j < middleTriangle.size(); j ++) {
                middle2source.push_back(getTransformMatrix(middleTriangle[j], sourceTriangle[j]));
                middle2dest.push_back(getTransformMatrix(middleTriangle[j], destTriangle[j]));
            }

            CImg<double> target(source.width(), source.height(), 1, 3, 0);
            
            cimg_forXY(target, x, y) {
                bool sourceFindFlag = false;
                bool destFindFlag = false;

                for (int j = 0; j < middleTriangle.size(); j ++) {
                    if (isPointInTriangle(Point(x, y, -1), middleTriangle[j])) {
                        int uSource = middle2source[j][0] * x + middle2source[j][1] * y + middle2source[j][2];
                        int vSource = middle2source[j][3] * x + middle2source[j][4] * y + middle2source[j][5];

                        int uDest = middle2dest[j][0] * x + middle2dest[j][1] * y + middle2dest[j][2];
                        int vDest = middle2dest[j][3] * x + middle2dest[j][4] * y + middle2dest[j][5]; 
                        
                        if (uSource >= 0 && uSource < source.width() && vSource >= 0 && vSource < source.height()) {
                            cimg_forC(source, channel) {
                                target(x, y, 0, channel) += source.linear_atXY(uSource, vSource, 0, channel) * ratio;
                            };
                            sourceFindFlag = true;
                        }

                        if (uDest >= 0 && uDest < dest.width() && vDest >= 0 && vDest < dest.height()) {
                            cimg_forC(dest, channel) {
                                target(x, y, 0, channel) += dest.linear_atXY(uDest, vDest, 0, channel) * (1 - ratio);
                            };
                            destFindFlag = true;
                        }
                        break;
                    }
                }
                // might lead to bug because of size issue
                if (!sourceFindFlag) {
                    cimg_forC(source, channel) {
                        target(x, y, 0, channel) += source.linear_atXY(x, y, 0, channel) * ratio;
                    };
                }
                if (!destFindFlag) {
                    cimg_forC(source, channel) {
                        target(x, y, 0, channel) += dest.linear_atXY(x, y, 0, channel) * (1 - ratio);
                    };
                }
            };
            target.save(("../transformProcess/" + to_string(interval - i) + ".jpg").c_str());
        }
    }

    void start() {
        getControlPoint();
        sourceTriangle = Delaunay::triangulate(sourceControlPoint);
        syncTriangle();
        drawTriangle();
        morph();
    }
};


int main(int argc, char* argv[]) {
    CImg<double> source(argv[1]);    
    CImg<double> dest(argv[2]);
    int interval = atoi(argv[3]);

    // resize source
    source = source.resize(dest.width(), dest.height(), 1, 3, 4);

    imageMorphing faceMorphing(source, dest, interval);
    faceMorphing.start();
}
