// Neil Mehra
// Period 5 Jurj

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

        double constGetX() const { 
            return x;
        }

        double constGetY() const {
            return y;
        }

		Point getShiftedPixel() {
			return Point(x, -y+799);
		}

		double getDistance(Point other) {
			return hypot(x-other.x, y-other.y);
		}

        bool operator < (const Point& p) const {
            return x < p.x;
        }

		string tostring() {
			return to_string(x) + ", " + to_string(y);
		}
};

class ClosestPointFinder {
    private:
        double minDist;
        vector<Point> l;
        ofstream file;
        Point closestFirst;
        Point closestSecond;

        void generatePoints(int amt) {
            srand(time(NULL));
            for(int i=0; i<amt; i++) {
                Point a = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
                l.push_back(a);
            }
        }

        void findMinDistBruteForce() {
            for(auto it=l.begin(); it!=l.end(); it++) {
                for(auto jt=next(it,1); jt!=l.end(); jt++) {
                    double dist = (*jt).getDistance(*it);
                    if(minDist > dist) {
                        minDist = dist;
                        closestFirst = *it;
                        closestSecond = *jt;
                    }
                }
            }
        }

        void findMinDistRecur() {
            sort(l.begin(), l.end());
            recur(l, 0, l.size()-1);
        }

        void findMinDistRecurFull() {
            sort(l.begin(), l.end());
            fullRecur(l, 0, l.size()-1);
        }

        pair<Point, Point> recur(vector<Point> &v, int start, int end) {
            if(abs(start-end)==1) {
                return make_pair(v[start], v[end]);
            }
            else if(abs(start-end)==2) {
                double first = v[start].getDistance(v[end]);
                double second = v[start].getDistance(v[start+1]);
                double third = v[start+1].getDistance(v[end]);
                double minThree = min(min(first, second), third);
                if(minThree == first) return make_pair(v[start], v[end]);
                if(minThree == second) return make_pair(v[start], v[start+1]);
                if(minThree == third) return make_pair(v[start+1], v[end]);
            }

            int half = (start+end) / 2;
            pair<Point, Point> one = recur(v, start, half);
            pair<Point, Point> two = recur(v, half, end);
            pair<Point, Point> closests;
            double d = -1;
            if(one.first.getDistance(one.second)) {
                d = one.first.getDistance(one.second);
                closests = one;
            } else {
                d = two.first.getDistance(two.second);
                closests = two;
            }
            double stripMiddle = (v[half].getX()+v[half+1].getX())/2;

            int leftStripIndex = half;
            int rightStripIndex = half;
            while(v[leftStripIndex].getX() < stripMiddle - d)  leftStripIndex--;
            while(v[rightStripIndex].getX() > stripMiddle + d) rightStripIndex++;

            for(int i=leftStripIndex; i<half; i++) {
                for(int j=half; j<rightStripIndex; j++) {
                    if(v[i].getDistance(v[j]) < d) {
                        d = v[i].getDistance(v[j]);
                        closests = make_pair(v[i], v[j]);
                    }
                }
            }

            if(d < minDist) {
                minDist = d;
                closestFirst = closests.first;
                closestSecond = closests.second;
            }
            return closests;
        }

        pair<Point, Point> bruteForceThree(vector<Point>& v, int start, int end) {
            double first = v[start].getDistance(v[end]);
            double second = v[start].getDistance(v[start+1]);
            double third = v[start+1].getDistance(v[end]);
            double minThree = min(min(first, second), third);
            if(minThree == first) return make_pair(v[start], v[end]);
            else if(minThree == second) return make_pair(v[start], v[start+1]);
            else return make_pair(v[start+1], v[end]);
        }

        double findSmallestDistance(pair<Point, Point> a, pair<Point, Point> b) {
            if(a.first.getDistance(a.second) < b.first.getDistance(b.second)) {
                return a.first.getDistance(a.second);
            } else {
                return b.first.getDistance(b.second);
            }
        }

        pair<Point, Point> findSmallestPair(pair<Point, Point> a, pair<Point, Point> b) {
            if(findSmallestDistance(a, b) == a.first.getDistance(a.second)) {
                return a;
            } else {
                return b;
            }
        }

        double getStripMiddle(vector<Point>& v, int half) {
            return (v[half].getX()+v[half+1].getX())/2;
        }

        double getLeftStripIndex(vector<Point>& v, double stripMiddle, double d, int half) {
            int leftStripIndex = half;
            while(leftStripIndex >= 0 && v[leftStripIndex].getX() < stripMiddle - d)  leftStripIndex--;
            return leftStripIndex;
        }

        double getRightStripIndex(vector<Point>& v, double stripMiddle, double d, int half) {
            int rightStripIndex = half;
            while(rightStripIndex < v.size() && v[rightStripIndex].getX() > stripMiddle + d) rightStripIndex++;
            return rightStripIndex;
        }


        vector<Point> getStrip(vector<Point>& v, int leftStripIndex, int rightStripIndex) {
            vector<Point> strip;
            for(int i=leftStripIndex; i<rightStripIndex; i++) {
                strip.push_back(v[i]);
            }
            sort(strip.begin(), strip.end(), [](const Point& lhs, const Point& rhs) {
                return lhs.constGetX() < rhs.constGetY();
            });
            return strip;
        }

        pair<Point, Point> fullRecur(vector<Point> &v, int start, int end) {
            if(abs(start-end)==1) {
                return make_pair(v[start], v[end]);
            }
            else if(abs(start-end)==2) {
                return bruteForceThree(v, start, end);
            }

            int half = (start+end) / 2;
            pair<Point, Point> one = fullRecur(v, start, half);
            pair<Point, Point> two = fullRecur(v, half, end);
            pair<Point, Point> closests = findSmallestPair(one, two);
            double d = findSmallestDistance(one, two);

            double stripMiddle = getStripMiddle(v, half);

            int leftStripIndex = getLeftStripIndex(v, stripMiddle, d, half);
            int rightStripIndex = getRightStripIndex(v, stripMiddle, d, half);

            vector<Point> strip = getStrip(v, leftStripIndex, rightStripIndex);
            for(int i=0; i<strip.size(); i++) {
                for(int j=i; j<i+15; j++) {
                    if(j<strip.size()&&i!=j&&strip[i].getDistance(strip[j]) < d) {
                        d = strip[i].getDistance(strip[j]);
                        closests = make_pair(strip[i], strip[j]);
                    }
                }
            }

            if(d < minDist) {
                minDist = d;
                closestFirst = closests.first;
                closestSecond = closests.second;
            }
            return closests;
        }
    public:
        ClosestPointFinder() {
            file.open("results.txt");
            file<<fixed<<showpoint<<setprecision(26);
            cout<<fixed<<showpoint<<setprecision(26);
            minDist = 1000000000;
            generatePoints(1000);
        }

        void part1() {
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistBruteForce();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Brute Force seconds: "<<seconds.count()<<endl;
            cout<<"Brute Force point 1: "<<closestFirst.tostring()<<endl;
            cout<<"Brute Force point 2: "<<closestSecond.tostring()<<endl;
            cout<<"Brute Force min dist: "<<minDist<<endl<<endl;

            file<<"Brute Force seconds: "<<seconds.count()<<endl;
            file<<"Brute Force point 1: "<<closestFirst.tostring()<<endl;
            file<<"Brute Force point 2: "<<closestSecond.tostring()<<endl;
            file<<"Brute Force min dist: "<<minDist<<endl<<endl;

        }

        void part2() {
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistRecur();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Recur seconds: "<<seconds.count()<<endl;
            cout<<"Recur point 1: "<<closestFirst.tostring()<<endl;
            cout<<"Recur point 2: "<<closestSecond.tostring()<<endl;
            cout<<"Recur min dist: "<<minDist<<endl<<endl;

            file<<"Recur seconds: "<<seconds.count()<<endl;
            file<<"Recur point 1: "<<closestFirst.tostring()<<endl;
            file<<"Recur point 2: "<<closestSecond.tostring()<<endl;
            file<<"Recur min dist: "<<minDist<<endl<<endl;
        }

        void part3() {
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistRecurFull();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Recur full seconds: "<<seconds.count()<<endl;
            cout<<"Recur full point 1: "<<closestFirst.tostring()<<endl;
            cout<<"Recur full point 2: "<<closestSecond.tostring()<<endl;
            cout<<"Recur full min dist: "<<minDist<<endl<<endl;

            file<<"Recur full seconds: "<<seconds.count()<<endl;
            file<<"Recur full point 1: "<<closestFirst.tostring()<<endl;
            file<<"Recur full point 2: "<<closestSecond.tostring()<<endl;
            file<<"Recur full min dist: "<<minDist<<endl<<endl;
        }
};


int main() {
    ClosestPointFinder finder = ClosestPointFinder();
    finder.part1();
    finder.part2();
    finder.part3();
}