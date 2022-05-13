// Neil Mehra
// Lab 02 part 1
// Jurj 5th Period
// 9/18/21
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

class Point {
	public: 
		double x, y;
		
		Point() {
			x = ((double)rand() / RAND_MAX);
			y = ((double)rand() / RAND_MAX);
		}
};

double getArea(Point a, Point b, Point c) {
	return abs(0.5 * (a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y)));
}

vector<Point> getPoints(double epsilon) {
	srand(time(NULL));
	Point a = Point();
	Point b = Point();
	Point c = Point();
	Point d = Point();
	while(abs(getArea(a,b,c)-(getArea(a,b,d)+getArea(a,c,d)+getArea(b,c,d)))<epsilon) {
		d = Point();
	}
	return vector<Point>{a,b,c,d};
}

int main() {
	freopen("points.txt","w",stdout);
	cout<<fixed<<showpoint<<setprecision(17);
    vector<Point> vec = getPoints(pow(10,-7));
	cout<<"("<<vec[0].x<<","<<vec[0].y<<")";
	for(int i=1; i<vec.size(); ++i) {
		cout<<" , "<<"("<<vec[i].x<<","<<vec[i].y<<")";
	}
}