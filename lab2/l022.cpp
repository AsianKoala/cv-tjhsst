// Neil Mehra
// Lab 02 part 2
// Jurj 5
// 10/8/21
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <climits>
#include <algorithm>
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

		double getX() {
			return x; 
		}

		double getY() {
			return y;
		}

		Point getShiftedPixel() {
			return Point(x, -y+799);
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
			bool in_x = p.getX() >= 0 && p.getX() < 799;
			bool in_y = p.getY() >= 0 && p.getY() < 799;
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
				int currPixelX = points[i].getShiftedPixel().getX();
				int currPixelY = points[i].getShiftedPixel().getY();
				

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

// used for swapping intenrals of points when we need to swap for drawing steep/shallow lines with bresenham
void swap(Point& a, Point& b) {
	Point temp = Point(b.getX(), b.getY());
	b = Point(a.getX(), a.getY());
	a = temp;
}

void drawshallowline(Point start, Point end, vector<Point>& vec) {
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
		vec.push_back(Point(x,y));
		if(e >= 0) {
			y += increment;
			e -= dx;
		}
		e+=dy;
	}
}

void drawsteepline(Point start, Point end, vector<Point>& vec) {
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
		vec.push_back(Point(x, y));
		if(e >= 0) {
			x += increment;
			e -= dy;
		}
		e += dx;
	}
}

// handles drawing of any type of line
void drawline(vector<Point>& points, Point start, Point end) {
	if(abs(end.getY() - start.getY()) < abs((end.getX() - start.getX()))) {
		if(start.getX() > end.getX()) {
			swap(start,end);
		}
		drawshallowline(start, end, points);
	} else {
		if(start.getY() > end.getY()) {
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
		
		points.push_back(Point((int)c.getX()+x, (int)c.getY()+y));
		points.push_back(Point((int)c.getX()-x, (int)c.getY()+y));
		points.push_back(Point((int)c.getX()+x, (int)c.getY()-y));
		points.push_back(Point((int)c.getX()-x, (int)c.getY()-y));
		points.push_back(Point((int)c.getX()+y, (int)c.getY()+x));
		points.push_back(Point((int)c.getX()-y, (int)c.getY()+x));
		points.push_back(Point((int)c.getX()+y, (int)c.getY()-x));
		points.push_back(Point((int)c.getX()-y, (int)c.getY()-x));

		y2_new -= (2 * x) - 3;
	}
}

// gets the area of a triangle bounded by three points: a, b, c
double getArea(Point a, Point b, Point c) {
	return abs(0.5 * (a.getX()*(b.getY()-c.getY())+b.getX()*(c.getY()-a.getY())+c.getX()*(a.getY()-b.getY())));
}

// creates 3 random points: a, b, c and gives d a new random value until it fits within the constraints of the problem
// if the area of triangle ABC is approximately equasl to the subtriangles created by ABD, ACD, and BCD, then set d 
// equal to a point with random x and y values
vector<Point> getPart1Points(double epsilon) {
	srand(time(NULL));
	Point a = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
	Point b = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
	Point c = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
	Point d = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
	double x;
	double y;
	while(abs(getArea(a,b,c)-(getArea(a,b,d)+getArea(a,c,d)+getArea(b,c,d)))<epsilon) {
		x = ((double)rand())/RAND_MAX;
		y = ((double)rand())/RAND_MAX;
		d = Point(x,y);
	}
	return vector<Point>{a,b,c,d};
}

// calculates smallest square created by the extensions of each of these points
// returns a vector containing the 4 intersection points, which is then used to draw one of the intersections
vector<Point> calculate(ofstream& outputfile, vector<Point> points, bool write) {
	Point a = points[0];
	Point b = points[1];
	Point c = points[2];
	Point d = points[3];

	Line ac = Line(a, c);
	double a_dist = hypot(a.getX()-c.getX(), a.getY()-c.getY());
	double p_x = d.getX() + sqrt((a_dist*a_dist)/(1+ac.getPerpendicularSlope()*ac.getPerpendicularSlope()));
	double p_y = ac.getPerpendicularSlope() * (p_x - d.getX()) + d.getY();
	Point p = Point(p_x, p_y);
	Line bm = Line(b, p);
	Line pline = Line(p, bm.getSlope());
	Line dline = Line(d, bm.getSlope());
	Line aline = Line(a, bm.getPerpendicularSlope());
	Line cline = Line(c, bm.getPerpendicularSlope());
	Point top_left = aline.intersect(dline);
	Point top_right = aline.intersect(pline);
	Point bottom_left = cline.intersect(dline);
	Point bottom_right = cline.intersect(pline);
	
	double area = top_left.distanceTo(top_right) * top_left.distanceTo(top_right);
	Point intersections[] = { top_left, top_right, bottom_left, bottom_right };
	
	if(write) {
		outputfile<<"("<<intersections[0].getX()<<", "<<intersections[0].getY()<<") ,";

		for(int i=1; i<3; i++) {
			outputfile<<" ("<<intersections[i].getX()<<", "<<intersections[i].getY()<<") ,";
		}
		outputfile<<" ("<<intersections[3].getX()<<", "<<intersections[3].getY()<<")"<<" AREA: "<<area<<endl;
	}

	return vector<Point>{top_left, top_right, bottom_left, bottom_right};
}

void calc_and_draw_line(vector<Point>& points, Point a, Point b) {
	Line l = Line(a, b);
	double maxy = 800 * l.getSlope() + l.getIntercept();
	drawline(points, Point(1, l.getIntercept()), Point(800, maxy));
}

void scaleToEightHundred(vector<Point>& v) {
	for(unsigned int i=0; i<v.size(); i++) {
		double newX = v[i].getX() * 800;
		double newY = v[i].getY() * 800;
		v[i].setX(newX);
		v[i].setY(newY);
	}
}

// basically find the lowest area by reading the output, and then read those points in by outtputting the 
// smallest points to another file and reading those
// then draw everything
void draw_everything(vector<Point> initialPoints) {
	ifstream points_file;
	points_file.open("output.txt");
	ofstream helper_file;
	helper_file.open("helper.txt");
	string s;
	int counter = 0;
	double min_area = INT_MAX;
	int lowest_counter = 0;
	while(getline(points_file,s)) {
		counter++;
		if(s.find("AREA")!=string::npos) {
			int index = s.find("AREA");
			double area = stod(s.substr(index+6));
			if(area < min_area) {
				min_area = area;
				lowest_counter = counter;
			}
		}
	}
	points_file.close();
	points_file.open("output.txt");
	counter = 0;
	while(getline(points_file,s)) {
		counter++;
		if(counter==lowest_counter) {
			helper_file<<s.substr(0,s.find("AREA"))<<endl;
		}
	}
	helper_file.close();

	vector<Point> square_points;
	ifstream helper_reader;
	helper_reader.open("helper.txt");
	double read_x=0;
	double read_y=0;
	for(int i=0; i<11; i++) {
		string s;
		helper_reader>>s;
		if(s==",") {
			continue;
		}

		if(s.at(0)=='(') {
			read_x = stod(s.substr(1));
		} else {
			read_y = stod(s);
			square_points.push_back(Point(read_x, read_y));
		}
	}

	VecPPM ppm = VecPPM("output.ppm", 800, 800);
	vector<Point> drawpoints;

	scaleToEightHundred(initialPoints);
	scaleToEightHundred(square_points);

	calc_and_draw_line(drawpoints, square_points[0], square_points[1]);
	calc_and_draw_line(drawpoints, square_points[3], square_points[2]);
	calc_and_draw_line(drawpoints, square_points[2], square_points[0]);
	calc_and_draw_line(drawpoints, square_points[1], square_points[3]);

	drawcircle(drawpoints, initialPoints[0], 2);
	drawcircle(drawpoints, initialPoints[1], 2);
	drawcircle(drawpoints, initialPoints[2], 2);
	drawcircle(drawpoints, initialPoints[3], 2);

	ppm.setPoints(drawpoints);
	ppm.build();
}

// helper method used to swap array values in the permute method
void arrSwap(vector<Point>& v, int i, int j) {
	Point temp = v[i];
	v[i] = v[j];
	v[j] = temp;
}

// used to get all possible permutations of the array (24 in this case)
void permute(ofstream& outputfile, vector<Point> v, int k) {
	for(int i = k; i < 4; i++){
		arrSwap(v, i, k);
		permute(outputfile, v, k+1);
		arrSwap(v, k, i);
	}

	if(k==3) {
		calculate(outputfile, v, true);
	}
}

// runs the main part1() program
// finds 4 valid points that can be used for part2()
void part1() {
	ofstream pointstxt;
	pointstxt.open("points.txt");
	pointstxt<<setprecision(17);
    vector<Point> vec = getPart1Points(pow(10,-7));
	pointstxt<<"("<<vec[0].getX()<<","<<vec[0].getY()<<")";
	for(unsigned int i=1; i<vec.size(); ++i) {
		pointstxt<<" , "<<"("<<vec[i].getX()<<","<<vec[i].getY()<<")";
	}
}

void part2() {
	ifstream pointsfile;
	pointsfile.open("points.txt");
	ofstream outputfile;
	outputfile.open("output.txt");
	outputfile<<setprecision(17);
	vector<Point> points;

	for(int i=0; i<4; i++) {
		if(i!=0) {
			string comma;
			pointsfile>>comma;
		}
		double x, y;
		string s;
		pointsfile>>s;
		int startIndex = 1;
		int endIndex = s.find(",")-1;
		int secondStartIndex = endIndex+2;

		x = stod(s.substr(startIndex,endIndex));
		y = stod(s.substr(secondStartIndex));
		points.push_back(Point(x,y));
	}
	outputfile<<"("<<points[0].getX()<<", "<<points[0].getY()<<") ,";

	for(int i=0; i<3; i++) {
		outputfile<<" ("<<points[i].getX()<<", "<<points[i].getY()<<") ,";
	}
	outputfile<<" ("<<points[3].getX()<<", "<<points[3].getY()<<")"<<endl;

	permute(outputfile, points, 0);
	
	pointsfile.close();
	outputfile.close();

	draw_everything(points);
}

int main() {
    // part1();
    part2();
}