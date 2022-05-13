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
#include <unordered_map>

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

		double getY() const {
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

        bool operator==(const Point& p) const {
            return x == p.x && y == p.y;
        }
};

class PointHashFunction {
    public:
        size_t operator()(const Point& p) const {
            // return (53 + hash<double>{}(p.getX())) * 53 + hash<double>{}(p.getY());
            return p.getX() + p.getY();
        }
};

class ClosestPointFinder {
    private:
        double minDist;
        vector<Point> l;
        ofstream file;
        ofstream otherfile;
        Point closestFirst;
        Point closestSecond;

        void reset() {
            minDist = 10000000;
            closestFirst = Point();
            closestSecond = Point();
        }

        void generatePoints(int amt) {
            srand(time(NULL));
            otherfile = ofstream("points.txt");
            otherfile<<fixed<<showpoint<<setprecision(26);
            for(int i=0; i<amt; i++) {
                Point a = Point(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX);
                l.push_back(a);
                otherfile<<a.getX()<<" "<<a.getY()<<endl;
            }
        }

        void readPoints(string name) {
            ifstream file(name);
            while(!file.eof()) {
                double x, y;
                file>>x>>y;
                l.push_back(Point((double)x,(double)y));
            }
            l.pop_back();     
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

        double getLeftStripIndex(vector<Point>& v, double stripMiddle, double d, int half, int start) {
            int leftStripIndex = half;
            while(leftStripIndex >= start && v[leftStripIndex].getX() > stripMiddle - d)  leftStripIndex--;
            return leftStripIndex;
        }

        double getRightStripIndex(vector<Point>& v, double stripMiddle, double d, int half, int end) {
            int rightStripIndex = half;
            while(rightStripIndex < end && v[rightStripIndex].getX() < stripMiddle + d) rightStripIndex++;
            return rightStripIndex;
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

            pair<Point, Point> closests = findSmallestPair(one, two);

            double d = findSmallestDistance(one, two);
            double stripMiddle = getStripMiddle(v, half);
            int leftStripIndex = getLeftStripIndex(v, stripMiddle, d, half, start);
            int rightStripIndex = getRightStripIndex(v, stripMiddle, d, half, end);


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

        pair<Point, Point> fullRecur(vector<Point>& v, int start, int end) {
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
            int leftStripIndex = getLeftStripIndex(v, stripMiddle, d, half, start);
            int rightStripIndex = getRightStripIndex(v, stripMiddle, d, half, end);

            vector<Point> strip(v.begin()+leftStripIndex, v.begin()+rightStripIndex);
            sort(strip.begin(), strip.end(), [](const Point& lhs, const Point& rhs) {
                return lhs.getY() < rhs.getY();
            });

            for(unsigned int i=0; i<strip.size(); i++) {
                for(unsigned int j=i; j<min((unsigned int)strip.size(),i+15); j++) {
                    if(i!=j&&strip[i].getDistance(strip[j]) < d) {
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


        void findMinDistPartFour() {
            // knuth shuffle
            for(int i=l.size()-1; i>0; i--) {
                int j = rand() % (i + 1);
                Point temp = l[i];
                l[i] = l[j];
                l[j] = temp;
            }
            pair<Point, Point> closest = partFour(l);
            closestFirst = closest.first;
            closestSecond =  closest.second;
            minDist = closestFirst.getDistance(closestSecond);
        }

        unordered_map<int, vector<Point>> MakeDictionary(vector<Point>& v, int num, unsigned long long N, double s) {
            unordered_map<int, vector<Point>> dict;
            for(int i=0; i<num; i++) {
                unsigned long long x = (unsigned long long)(v[i].getX() / s);
                unsigned long long y = (unsigned long long)(v[i].getY() / s);
                unsigned long long val = x * N + y;
                if(dict.find(val) != dict.end()) {
                    dict[val].push_back(v[i]);
                } else {
                    dict[val] = vector<Point>{v[i]};
                }
            }
            return dict;
        }

        pair<Point, Point> partFour(vector<Point>& v) {
            pair<Point, Point> closest = make_pair(v[0], v[1]);
            double delta = v[0].getDistance(v[1]);
            double s = delta / 2;
            unsigned long long N = 100000;
            unordered_map<int, vector<Point>> dict = MakeDictionary(v, 0, N, s);

            for(int i=0; i<v.size(); i++) {
                unsigned long long x = (unsigned long long)(v[i].getX() / s);
                unsigned long long y = (unsigned long long)(v[i].getY() / s);

                double lowestDistanceInClose = 1;
                Point closestPointInClose;
                for(unsigned long long s = x-2; s <= x+2; s++) {
                    for(unsigned long long t = y-2; t <= y+2; t++) {
                        unsigned long long nVal = s * N + t;
                        if(dict.find(nVal) != dict.end()) {
                            for(Point j : dict[nVal]) {
                                if(v[i].getDistance(j) < lowestDistanceInClose) {
                                    lowestDistanceInClose = v[i].getDistance(j);
                                    closestPointInClose = j;
                                }
                            }
                        }
                    }
                }

                if(lowestDistanceInClose < delta) {
                    delta = lowestDistanceInClose;
                    s = delta / 2;
                    closest = make_pair(v[i], closestPointInClose);
                    dict = MakeDictionary(v, i, N, s);
                } else {
                    unsigned long long val = x * N + y;
                    dict[val].push_back(v[i]);
                }
            }
            return closest;
        }

            /*
            first stage set delta = d(p1, p2)
            goal of stage: verify delta is distance between closest points, or find pair pi, pj, with d(pi, pj) < delta
            terminate stage when: reach pi, such that j < i, we have d(pi, pj) < delta
            then let delta = closest distance found so far for next stage
            verification of delta being smallest:
            divide unit square into squares with side length delta / 2
            there will be N^2 subsquares, with N=[1/(2 * delta)]
            for 0 <= s <= N-1
            and 0 <= t <= N-1
            subsquare Sst:
            x = s delta / 2 <= x < (s+1) delta / 2
            y = t delta / 2 <= y < (t+1) delta / 2
            any points that lie in same subsquare have distance less than delta
            any two points less than delta away must fall in same subsquare
            or be in very close subsquares
            subsquare Sst and Ss't' are close if |s-s'| <= 2 and |t-t'| <= 2
            a subsquare is close to itself
            there are at most 25 subsquares close to Sst, including Sst
            there will be less than 25 is Sst is at the edge of unit square
            for each point in P', keep track of subsquares containing it
            when next point p is considered, check which subsquares Sst p belongs to
            if p is going to cause minimum distancec to change, there must be some earlier point
            p' in P' at distance less than delta from it,
            hence that point p' is in one of the 25 squares around Sst containing p
            check each of 25 squares
            */


    public:
        ClosestPointFinder() {
            file.open("results.txt");
            file<<fixed<<showpoint<<setprecision(26);
            cout<<fixed<<showpoint<<setprecision(26);
            reset();
            generatePoints(128*2*2*2*2*2*2*2*2*2*2*2*2);
            readPoints("points.txt");
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
            reset();
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistRecur();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Recur seconds: "<<seconds.count()<<endl;
            cout<<"Recur point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            cout<<"Recur point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            cout<<"Recur min dist: "<<minDist<<endl<<endl;

            file<<"Recur seconds: "<<seconds.count()<<endl;
            file<<"Recur point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            file<<"Recur point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            file<<"Recur min dist: "<<minDist<<endl<<endl;
        }

        void part3() {
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistRecurFull();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Part 3 seconds: "<<seconds.count()<<endl;
            cout<<"Part 3 point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            cout<<"Part 3 point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            cout<<"Part 3 min dist: "<<minDist<<endl<<endl;

            file<<"Part 3 full seconds: "<<seconds.count()<<endl;
            file<<"Part 3 point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            file<<"Part 3 point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            file<<"Part 3 min dist: "<<minDist<<endl<<endl;
        }

        void part4() {
            reset();
            chrono::time_point<chrono::system_clock> start, end;
            chrono::duration<double> seconds;
            start = chrono::system_clock::now();
            findMinDistPartFour();
            end = chrono::system_clock::now();
            seconds = end - start;
            cout<<"Part 4 full seconds: "<<seconds.count()<<endl;
            cout<<"Part 4 point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            cout<<"Part 4 point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            cout<<"Part 4 min dist: "<<minDist<<endl<<endl;

            file<<"Part 4 full seconds: "<<seconds.count()<<endl;
            file<<"Part 4 point 1: "<<closestFirst.getX()<<" , "<<closestFirst.getY()<<endl;
            file<<"Part 4 point 2: "<<closestSecond.getX()<<" , "<<closestSecond.getY()<<endl;
            file<<"Part 4 min dist: "<<minDist<<endl<<endl;
        }

};


int main() {
    ClosestPointFinder finder = ClosestPointFinder();
    // finder.part1();
    // finder.part2();
    finder.part3();
    finder.part4();
}