// Neil Mehra
// Lab 01
// Jurj 5th Period
// 9/7/21
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime> 
using namespace std;

// points have ints because they represent pixels but origin is bottom left not top left
class Point {
	public: 
		int x, y;
		
		Point(int nx, int ny) {
			x = nx;
			y = ny;
		}

		Point getpixel() {
			return Point(x, -y+799);
		}

		string tostring() {
			return to_string(x) + ", " + to_string(y);
		}

		bool operator < (const Point& p) const {
			return (y > p.y) || (y==p.y && x < p.x);
		}
};

// Circle is a data class used for in between operations, which is why it uses doubles
class Circle {
	public:
		double x, y, r;

		Circle(double dx, double dy, double dr) {
			x = dx;
			y = dy;
			r = dr;
		}
};

// custom ppm class that handles filtering -> drawing vector 
class VecPPM {
	private:
		int dimx, dimy;
		ofstream file;
		vector<Point> points;

		bool in_bounds(Point p) {
			bool in_x = p.x >= 0 && p.x < 799;
			bool in_y = p.y >= 0 && p.y < 799;
			return in_x && in_y;
		}

		void setupvec() {
			// sort so biggest y and smallest x is at first position
			sort(points.begin(), points.end());

			// remove points that are out of bounds or remove points that are duplicates of next point (to avoid multiple 0 0 0 in same pixel spot in ppm)
			vector<Point>::iterator it = points.begin();
			while(it != points.end()) {
				if(!in_bounds((*it)) || (*it).x == (*next(it)).x && (*it).y == (*next(it)).y) {
					it = points.erase(it);
				} else {
					++it;
				}
			}
		}

		void add_pixels(int amt, int& row) {
			int i = 0;
			while(i<amt) {
				file<<"1 1 1 ";
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

		vector<Point>& getpoints() {
			return points;
		}

		void build() {
			int row = 0;
			int col = 0;

			setupvec();
			// max col should be dimx*6+1
			// max row should be dimy+1 (cause header)
			for(int i=0; i<points.size(); ++i, ++col) {
				Point currPixel = points[i].getpixel();

				// coordinatse are 0 indexed
				// check if current row needs to add any to be full
				if(col < dimx && currPixel.y != row) {
					add_pixels(dimx-col, row);
					col = 0;
				}
				// now we point to new line

				// check if we need to skip a row
				while(row < currPixel.y) {
					add_pixels(dimx, row);
				}
			
				// move to currpixel
				while(col < currPixel.x) {
					file<<"1 1 1 ";
					col++;
				}
				file<<"0 0 0 ";
			}

			add_pixels(dimx-col, row);
			while(row < dimy) {
				add_pixels(dimx, row);
			}

			file.close();
		}
};

void swap(Point& a, Point& b) {
	Point temp = Point(b.x, b.y);
	b = Point(a.x, a.y);
	a = temp;
}

void drawshallowline(Point start, Point end, vector<Point>& vec) {
	int dx = end.x - start.x;
	int dy = end.y - start.y;
	int increment = 1;
	int y = start.y;
	int e = dy - dx;

	if(dy < 0) {
		increment = -1;
		dy = -dy;
	}

	for(int x=start.x; x<end.x; x++) {
		vec.push_back(Point(x,y));
		if(e >= 0) {
			y += increment;
			e -= dx;
		}
		e+=dy;
	}
}

void drawsteepline(Point start, Point end, vector<Point>& vec) {
	int dx = end.x - start.x;
	int dy = end.y - start.y;
	int increment = 1;
	int x = start.x;
	int e = dx - dy;

	if(dx < 0) {
		increment = -1;
		dx = -dx;
	}

	for(int y=start.y; y<end.y; y++) {
		vec.push_back(Point(x, y));
		if(e >= 0) {
			x += increment;
			e -= dy;
		}
		e += dx;
	}
}

void drawline(vector<Point>& points, Point start, Point end) {
	// check for steep/shallowness ignoring negative slope
	// then fix negative slope by swapping points if start driver > end driver 
	if(abs(end.y - start.y) < abs((end.x - start.x))) {
		if(start.x > end.x) {
			swap(start,end);
		}
		drawshallowline(start, end, points);
	} else {
		if(start.y > end.y) {
			swap(start,end);
		}
		drawsteepline(start, end, points);
	}
}

void drawcircle(vector<Point>& points, Point c, int c_r) {
	int xmax = (int) (c_r * 0.70710678);
	int y = c_r;
	int y2 = y * y;
	int ty = (2 * y) - 1;
	int y2_new = y2;

	for(int x = 0; x <= xmax+2; x++) {
		if((y2 - y2_new) >= ty) {
			y2 -= ty;
			y--;
			ty -= 2;
		}
		
		points.push_back(Point(c.x+x, c.y+y));
		points.push_back(Point(c.x-x, c.y+y));
		points.push_back(Point(c.x+x, c.y-y));
		points.push_back(Point(c.x-x, c.y-y));
		points.push_back(Point(c.x+y, c.y+x));
		points.push_back(Point(c.x-y, c.y+x));
		points.push_back(Point(c.x+y, c.y-x));
		points.push_back(Point(c.x-y, c.y-x));

		y2_new -= (2 * x) - 3;
	}
}

void calc_triangle(vector<Point>& points, Point a, Point b, Point c) {
	drawline(points, b, c);
	drawline(points, a, b);
	drawline(points, a, c);
}

void calc_incircle(vector<Point>& points, Point a, Point b, Point c) {
	double in_a = hypot(b.x - c.x, b.y - c.y);
	double in_b = hypot(c.x - a.x, c.y - a.y);
	double in_c = hypot(a.x - b.x, a.y - b.y);
	
	double xsum = in_a * a.x + in_b * b.x + in_c * c.x;
	double ysum = in_a * a.y + in_b * b.y + in_c * c.y;
	double dsum = in_a + in_b + in_c;

	double s = (in_a + in_b + in_c) / 2.0;

	// ints for result
	double inc_x = xsum / dsum;
	double inc_y = ysum / dsum;
	double inc_r = sqrt(((s - in_a) * (s - in_b) * (s - in_c)) / s);

	drawcircle(points, Point(inc_x, inc_y), inc_r);
}

Circle calc_cc(vector<Point>& points, Point a, Point b, Point c) {
	double cc_x, cc_y, cc_r;

	// generate coordinates for a circle that goes through a,b,c
	double A = b.x - a.x;
	double B = b.y - a.y;
	double C = c.x - a.x;
	double D = c.y - a.y;
	double E = A * (a.x + b.x) + B * (a.y + b.y);
	double F = C * (a.x + c.x) + D * (a.y + c.y);
	double G = 2 * (A * (c.y - b.y) - B * (c.x - b.x));
	double minx, miny, dx, dy;


	// avoid division by zero if points are on same line
	if(G==0) {
		G+=0.0003;
	}

	cc_x = (D * E - B * F) / G;
	cc_y = (A * F - C * E) / G;
	dx = cc_x - a.x;
	dy = cc_y - a.y;
	// distance betwen each point and cc is same so just use a
	cc_r = hypot(dx, dy);

	// truncate doubles to ints by passing into Point, and double into int parameter
	drawcircle(points, Point(cc_x, cc_y), cc_r);

	// return doubles for in betweens
	return Circle(cc_x, cc_y, cc_r);
}

void calc_ninecircle(vector<Point>& points, Circle ortho, Circle cc) {
	double x = (ortho.x + cc.x) / 2.0;
	double y = (ortho.y + cc.y) / 2.0;
	double r = cc.r / 2.0;
	drawcircle(points, Point(x, y), r);
}

void calc_euler(vector<Point>& points, Circle centroid, Circle cc) {
	// avvoid division by zero
	if(centroid.x == cc.x) {
		centroid.x += 0.003;
	}

	double dy = cc.y - centroid.y;
	double dx = cc.x - centroid.x;
	double m = dy / dx;
	double b = cc.y - m * cc.x;

	// now find farthest away points on euler line for most accurate pixel representation
	double miny = b;
	double maxy = 800 * m + b;
	drawline(points, Point(1, miny), Point(800, maxy));
}

Circle calc_centroid(Point a, Point b, Point c) {
	double centroid_x = (a.x + b.x + c.x) / 3.0;
	double centroid_y = (a.y + b.y + c.y) / 3.0;
	return Circle(centroid_x, centroid_y, 0);
}

Circle calc_ortho(Circle cc, Circle centroid) {
	double ortho_x = 3 * centroid.x - 2 * cc.x;
	double ortho_y = 3 * centroid.y - 2 * cc.y;
	return Circle(ortho_x, ortho_y, 0);
}

void calc_circles(vector<Point>& points, Point a, Point b, Point c) {

	// draw incircle
	calc_incircle(points, a, b, c);

	// now get circumcenter, tnts for final result
	Circle cc = calc_cc(points, a, b, c);

	// get centroid for euler line + ortho
	Circle centroid = calc_centroid(a, b, c);
	
	// calculate ortho center to get 9 circle center
	Circle ortho = calc_ortho(cc, centroid);

	calc_ninecircle(points, ortho, cc);

	calc_euler(points, centroid, cc);
}

void run_calcs(vector<Point>& points) {
	srand (time(NULL));

	// use a b and c for triangle verticies
	Point a = Point((rand()%100)*8, (rand()%100)*8);
	Point b = Point((rand()%100)*8, (rand()%100)*8);
	Point c = Point((rand()%100)*8, (rand()%100)*8);

	calc_triangle(points, a, b, c);
	calc_circles(points, a, b, c);
}

int main() {	
	VecPPM ppm = VecPPM("output.ppm", 800, 800);
	run_calcs(ppm.getpoints());	
	ppm.build();
}
