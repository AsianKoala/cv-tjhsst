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
#include <stack>

using namespace std;

class Point {
	private: 
		double x;
		double y;

	public:
		Point(double nx, double ny) {
			x = nx;
			y = ny;
		}

		Point() {
			x = 0;
			y = 0;
		}

		void setX(double nx) {
			x = nx;
		}

		void setY(double ny) {
			y = ny;
		}

		double getX() const {
			return x; 
		}

		double getY() {
			return y;
		}

		Point getShiftedPixel(int shift) {
			return Point(x, -y+shift);
		}

		double distanceTo(Point other) {
			return hypot(x-other.x, y-other.y);
		}

		bool operator < (const Point& p) const {
			return (y > p.y) || (y==p.y && x < p.x);
		}

		string tostring() {
			return to_string(x) + ", " + to_string(y);
		}
};

class Line {
	private:
		double slope, intercept;
		bool vertical;

	public:
		Line(Point point1, Point point2) {
			double dx = point1.getX() - point2.getX();
			double dy = point1.getY() - point2.getY();
			if(dx==0) {
				vertical = true;
			} else {
				slope = dy / dx;
				intercept = point1.getY() - point1.getX() * slope;
				vertical = false;
			}
		}

		Line(Point point, double s) {
			slope = s;
			intercept = point.getY() - point.getX() * s;
			vertical = false;
		}

		double getSlope() {
			return slope;
		}

		double getPerpendicularSlope() {
			return -1.0/slope;
		}

		double getIntercept() {
			return intercept;
		}

		double eval(double x) {
			return x * slope + intercept;
		}

		bool isVertical() {
			return vertical;
		}

		Point intersect(Line l2) {
			if(vertical) {
				return Point(intercept, l2.eval(intercept));
			} else if(l2.vertical) {
				return Point(l2.intercept, eval(intercept));
			} else {
				double xIntersect = (l2.intercept - intercept) / (slope - l2.slope);
				return Point(xIntersect, eval(xIntersect));
			}
		}
};

class VecPPM {
	private:
		int dimx, dimy;
		ofstream file;
		vector<Point> points;

		bool in_bounds(Point p) {
			bool in_x = p.getX() >= 0 && p.getX() < dimx-1;
			bool in_y = p.getY() >= 0 && p.getY() < dimy-1;
			return in_x && in_y;
		}

		void setupvec() {
			// sort so biggest y and smallest x is at first position
			sort(points.begin(), points.end());

			// remove points that are out of bounds or remove points that are duplicates of next point (to avoid multiple 0 0 0 in same pixel spot in ppm)
			vector<Point>::iterator it = points.begin();
			while(it != points.end()) {
				if(!in_bounds((*it)) || (((*it).getX() == (*next(it)).getX()) && ((*it).getY() == (*next(it)).getY()))) {
					it = points.erase(it);
				} else {
					++it;
				}
			}
		}

		void add_pixels(int amt, int& row) {
			int i = 0;
			while(i<amt) {
				file<<"0 0 0 ";
				++i;
			}
			file<<endl;
			row++;
		}




	public:
		VecPPM(string filename, int x, int y) {
			file.open(filename);
			file<<"P3 "<<x<<" "<<y<<" 1"<<endl;
			dimx = x;
			dimy = y;
		}

		vector<Point>& getPoints() {
			return points;
		}

		void setPoints(vector<Point> v) {
			points = v;
		}

		void build() {
			int row = 0;
			int col = 0;


			setupvec();
			// max col should be dimx*6+1
			// max row should be dimy+1 (cause header)
			for(unsigned int i=0; i<points.size(); ++i, ++col) {
				int currPixelX = points[i].getShiftedPixel(dimy-1).getX();
				int currPixelY = points[i].getShiftedPixel(dimy-1).getY();
				

				// coordinatse are 0 indexed
				// check if current row needs to add any to be full
				if(col < dimx && currPixelY != row) {
					add_pixels(dimx-col, row);
					col = 0;
				}
				// now we point to new line

				// check if we need to skip a row
				while(row < currPixelY) {
					add_pixels(dimx, row);
				}
			
				// move to currpixel
				while(col < currPixelX) {
					file<<"0 0 0 ";
					col++;
				}
				file<<"1 1 1 ";
			}

			add_pixels(dimx-col, row);
			while(row < dimy) {
				add_pixels(dimx, row);
			}

			file.close();
		}
};

VecPPM ppm = VecPPM("grahamscan.ppm", 400, 400);
vector<Point> points;
vector<Point> drawpoints;
vector<Point> convexHull;
vector<Point> scaledHull;
void part1();
void FindHull(vector<Point>& sk, Point p, Point q);
void drawHull();

int N = 60;
stack<Point> graham;
vector<Point> scaledGraham;
void part2();
void drawGraham();

void swap(Point& a, Point& b) {
	Point temp = Point(b.getX(), b.getY());
	b = Point(a.getX(), a.getY());
	a = temp;
}

void drawshallowline(Point start, Point end) {
	int dx = end.getX() - start.getX();
	int dy = end.getY() - start.getY();
	int increment = 1;
	int y = start.getY();
	int e = dy - dx;

	if(dy < 0) {
		increment = -1;
		dy = -dy;
	}

	for(int x=start.getX(); x<end.getX(); x++) {
		drawpoints.push_back(Point(x,y));
		if(e >= 0) {
			y += increment;
			e -= dx;
		}
		e+=dy;
	}
}

void drawsteepline(Point start, Point end) {
	int dx = end.getX() - start.getX();
	int dy = end.getY() - start.getY();
	int increment = 1;
	int x = start.getX();
	int e = dx - dy;

	if(dx < 0) {
		increment = -1;
		dx = -dx;
	}

	for(int y=start.getY(); y<(int)end.getY(); y++) {
		drawpoints.push_back(Point(x, y));
		if(e >= 0) {
			x += increment;
			e -= dy;
		}
		e += dx;
	}
}

// handles drawing of any type of line
void drawline(Point start, Point end) {
	if(abs(end.getY() - start.getY()) < abs((end.getX() - start.getX()))) {
		if(start.getX() > end.getX()) {
			swap(start,end);
		}
		drawshallowline(start, end);
	} else {
		if(start.getY() > end.getY()) {
			swap(start,end);
		}
		drawsteepline(start, end);
	}
}

void drawcircle(Point c, int c_r) {
	int xmax = (int) (c_r * 0.70710678);
	int y = c_r;
	int y2 = y * y;
	int ty = (2 * y) - 1;
	int y2_new = y2;

	for(int x = 0; x <= xmax; x++) {
		if((y2 - y2_new) >= ty) {
			y2 -= ty;
			y--;
			ty -= 2;
		}
		
		drawpoints.push_back(Point((int)c.getX()+x, (int)c.getY()+y));
		drawpoints.push_back(Point((int)c.getX()-x, (int)c.getY()+y));
		drawpoints.push_back(Point((int)c.getX()+x, (int)c.getY()-y));
		drawpoints.push_back(Point((int)c.getX()-x, (int)c.getY()-y));
		drawpoints.push_back(Point((int)c.getX()+y, (int)c.getY()+x));
		drawpoints.push_back(Point((int)c.getX()-y, (int)c.getY()+x));
		drawpoints.push_back(Point((int)c.getX()+y, (int)c.getY()-x));
		drawpoints.push_back(Point((int)c.getX()-y, (int)c.getY()-x));

		y2_new -= (2 * x) - 3;
	}
}

void generatePoints() {
    // srand(time(NULL));
    // ofstream outfile("points.txt");
    // outfile<<fixed<<showpoint<<setprecision(26);
    // for(int i=0; i<N; i++) {
    //     Point a = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
    //     points.push_back(a);
    //     outfile<<a.getX()<<" "<<a.getY()<<endl;
    // }
	ifstream infile("points.txt");
	for(int i=0; i<N; i++) {
		string x, y; infile>>x>>y;
		points.push_back(Point(stod(x),stod(y)));
	}
}

// returns if C is right of oriented line from a->b
bool isRightOfOrientedLine(Point s, Point a, Point b) {
	double temp = (b.getX()-a.getX())*(s.getY()-a.getY())-(b.getY()-a.getY())*(s.getX()-a.getX());
	return temp<0;
}

double triangleArea(Point a, Point b, Point c) {
	return abs(0.5 * (a.getX()*(b.getY()-c.getY())+b.getX()*(c.getY()-a.getY())+c.getX()*(a.getY()-b.getY())));
}

bool pointInTriangle(Point s, Point a, Point b, Point c) {
    return abs((triangleArea(a, b, c)-(triangleArea(s,b,c) + triangleArea(s,a,b) + triangleArea(s,a,c))))<0.0000001;
}

double lineDist(Point p, Point p1, Point p2) {
	return abs((p.getY() - p1.getY()) * (p2.getX() - p1.getX()) - 
			(p2.getY() - p1.getY()) * (p.getX() - p1.getX()));
}

void part1() {
    generatePoints();

	sort(points.begin(), points.end(), [](const Point& lhs, const Point& rhs) {
		return lhs.getX() < rhs.getX();
	});

	Point a = points[0];
	Point b = points.back();

	convexHull.push_back(a);
	convexHull.push_back(b);
	
    vector<Point> s1;
    vector<Point> s2;

    for(int i=1; i<points.size()-1; i++) {
		if(isRightOfOrientedLine(points[i], a, b)) {
			s1.push_back(points[i]);
		} else {
			s2.push_back(points[i]);
		}
	}

    FindHull(s1, a, b);

    for(int i=0; i<points.size(); i++) {
		drawcircle(Point(points[i].getX() * 400, points[i].getY() * 400), 3);
	}

	drawHull();
	convexHull = {a,b};
	scaledHull.clear();

    FindHull(s2, b, a);
	drawHull();

	ppm.setPoints(drawpoints);
	ppm.build();
}


void FindHull(vector<Point>& sk, Point p, Point q) {
    double maxDistance = lineDist(sk[0], p, q);
	int f = 0;
	for(int i=0; i<sk.size(); i++) {
		double d = lineDist(sk[i], p, q);
		if(d > maxDistance) {
			maxDistance = d;
			f = i;
		}
	}

	Point c = sk[f];

	convexHull.push_back(c);

	vector<Point> s1, s2;
	for(int i=0; i<sk.size(); i++) {
		if(i != f) {
			if(!pointInTriangle(sk[i], p, c, q)) {
				if(isRightOfOrientedLine(sk[i], p, c)) {
					s1.push_back(sk[i]);
				} else if(isRightOfOrientedLine(sk[i], c, q)) {
					s2.push_back(sk[i]);
				}
			} 
		}
	}
	
	if(s1.size() != 0) FindHull(s1, p, c);
	if(s2.size() != 0) FindHull(s2, c, q);
}


void drawHull() {
	for(int i=0; i<convexHull.size(); i++) {
		scaledHull.push_back(Point(convexHull[i].getX() * 400, convexHull[i].getY() * 400));
	}
	
	sort(scaledHull.begin(), scaledHull.end(), [](const Point& lhs, const Point& rhs) {
		return lhs.getX() < rhs.getX();
	});

	for(int i=0; i<scaledHull.size()-1; i++) {
		drawline(scaledHull[i], scaledHull[i+1]);
	}
}

double det(Point a, Point b, Point c) {
	return (a.getX() - b.getX()) * (c.getY() - b.getY()) - (a.getY() - b.getY()) * (c.getX() - b.getX());
}
Point pivot;
int ccw(Point a, Point b, Point c) {
	double determinant = det(a,b,c);
	if(determinant < 0) return -1;
	else if(determinant > 0) return 1;
	else return 0;
}

double dist(Point a, Point b) {
	return (a.getX() - b.getX()) * (a.getX() - b.getX()) + (a.getY() - b.getY()) * (a.getY() - b.getY());
}

Point nextToTop() {
	Point temp = graham.top();
	graham.pop();
	Point ret = graham.top();
	graham.push(temp);
	return ret;
}

void part2() {
	generatePoints();

	int lowestY = 0;
	for(int i=1; i<N; i++) {
		if((points[i].getY() < points[lowestY].getY()) ||
		(points[i].getY() == points[lowestY].getY() && points[i].getX() < points[lowestY].getX())) {
			lowestY = i;
		}
	}

	Point temp = points[0];
	points[0] = points[lowestY];
	points[lowestY] = temp;

	pivot = points[0];
	sort(points.begin() + 1, points.end(), [](const Point& a, const Point& b) {
		int val = ccw(pivot, a, b);
		if(val == 0) return dist(pivot, a) < dist(pivot, b);
		return val == 1;
	});

	graham.push(points[0]);
	graham.push(points[1]);
	graham.push(points[2]);

	for(int i=3; i<N; i++) {
		while(ccw(nextToTop(), graham.top(), points[i]) <= 0) {
			graham.pop();
		}
		graham.push(points[i]);
	}

	drawGraham();
}

void drawGraham() {
	for(int i=0; i<points.size(); i++) {
		drawcircle(Point(points[i].getX() * 400, points[i].getY() * 400), 3);
	}

	while(!graham.empty()) {
		Point p = graham.top();
		graham.pop();
		scaledGraham.push_back(Point(p.getX() * 400, p.getY() * 400));
	}

	for(int i=0; i<scaledGraham.size()-1; i++) {
		drawline(scaledGraham[i], scaledGraham[i+1]);
	}

	drawline(scaledGraham[0], scaledGraham[scaledGraham.size()-1]);

	ppm.setPoints(drawpoints);
	ppm.build();
}



int main() {
    // part1();
	part2();
}