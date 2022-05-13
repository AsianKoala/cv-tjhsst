#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <list>
#include <string>
#include <utility>
#include <chrono>
#include <ctime>
using namespace std;

class Point {
    public:
    	double x;
		double y;
		Point(double nx, double ny) {
			x = nx;
			y = ny;
		}

		Point() {
			x = 0;
			y = 0;
		}

        double getY() const {
            return y;
        }
};

int main() {
    vector<Point> v = {Point(0,1), Point(0,2), Point(0,3), Point(0,0), Point(0,-1)};
    for(int i=0; i<v.size(); i++) cout<<v[i].x<<" , "<<v[i].y<<endl;
    sort(v.begin(), v.end(), [](const Point& lhs, const Point& rhs) {
        return lhs.getY() < rhs.getY();
    });
    for(int i=0; i<v.size(); i++) cout<<v[i].x<<" , "<<v[i].y<<endl;
}