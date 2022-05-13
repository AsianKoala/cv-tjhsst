// neil mehra
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iomanip>
#include <string>
using namespace std;

class Pixel {
    private:
        int r,g,b;

    public:
        Pixel() {
            r = 0;
            g = 0;
            b = 0;
        }

        Pixel(int rr, int gg, int bb) {
            r = rr;
            g = gg;
            b = bb;
        }

        int getR() {
            return r;
        }

        int getG() {
            return g;
        }

        int getB() {
            return b;
        }

        void setR(int rr) {
            r = rr;
        }

        void setG(int gg) {
            g = gg;
        }

        void setB(int bb) {
            b = bb;
        }

        double linearify(int number) {
            double scaled = ((double)number) / 255.0;
            if(scaled <= 0.04045) {
                return 255 * (scaled / 12.92);
            } else {
                return 255 * pow(((scaled + 0.055)/(1.055)),2.4);
            }
        }

        Pixel getGrayscalePixel() {
            double r_linear = linearify(r);
            double g_linear = linearify(g);
            double b_linear = linearify(b);
            double y_linear = 0.2126 * r_linear + 0.7152 * g_linear + 0.0722 * b_linear;
            return Pixel(y_linear, y_linear, y_linear);
        }

        string toString() {
            return to_string(r) + " " + to_string(g) + " " + to_string(b);
        }
};

string filename = "image.ppm";
string ppmType;
int width, height;
int rgbSize;
ifstream file;
Pixel BLACK = Pixel();
Pixel WHITE = Pixel(255,255,255);
vector<vector<Pixel>> initialImage;
vector<vector<Pixel>> grayscaleImage;
vector<vector<double>> gMag;
vector<vector<double>> gAngle;

double inputWeakThreshold = 72;
double inputHighThreshold = 90;

void outputImage(string name, vector<vector<Pixel>> image) {
    ofstream output(name);
    output<<ppmType<<" "<<width<<" "<<height<<" "<<rgbSize<<endl;
    for(int h=0; h<height; h++) {
        for(int w=0; w<width; w++) {
            output<<image[h][w].toString()<<" ";
        }
        output<<endl;
    }
}

void outputBinaryImage(string name, vector<vector<Pixel>> image) {
    ofstream output(name);
    output<<"P1"<<" "<<width<<" "<<height<<endl;
    for(int h=0; h<height; h++) {
        for(int w=0; w<width; w++) {
            if(image[h][w].getR() == 255) {
                output<<"0";
            } else {
                output<<"1";
            }
            output<<" ";
        }
        output<<endl;
    }
}

void init() {
    ifstream file(filename);
    file>>ppmType;
    string tempWidth, tempHeight;
    file>>tempWidth>>tempHeight;
    width = stoi(tempWidth);
    height = stoi(tempHeight);
    string tempRGBSize;
    file>>tempRGBSize;
    rgbSize = stoi(tempRGBSize);

    for(int h=0; h<height; h++) {
        initialImage.push_back(vector<Pixel>());
        grayscaleImage.push_back(vector<Pixel>());
        gMag.push_back(vector<double>());
        gAngle.push_back(vector<double>());
        for(int w=0; w<width; w++) {
            int r,g,b; file>>r>>g>>b;
            initialImage[h].push_back(Pixel(r,g,b));
            grayscaleImage[h].push_back(Pixel());
            gMag[h].push_back(0);
            gAngle[h].push_back(0);
        }
    }
}


void grayscale() {
    for(int h=0; h<height; h++) {
        for(int w=0; w<width; w++) {
            grayscaleImage[h][w] = initialImage[h][w].getGrayscalePixel();
        }
    }
    outputImage("imageg.ppm", grayscaleImage);
}

void sobel() {
    double GxKernal[3][3] = 
        {
            {1, 0, -1},
            {2, 0, -2},
            {-1, 0, 1}
        };

    double GyKernal[3][3] = 
        {
            {1, 2, 1},
            {0, 0, 0},
            {-1, -2, -1}
        };

    for(int y=1; y<height-1;  y++) {
        for(int x=1; x<width-1; x++) {
            double accumX = 0;
            double accumY = 0;
            for(int i=-1; i<=1; i++) {
                for(int j=-1; j<=1; j++) {
                    int newX = x + i;
                    int newY = y + j;
                    Pixel p = initialImage[newY][newX];

                    accumX += p.getR() * GxKernal[i+1][j+1];
                    accumY += p.getR() * GyKernal[i+1][j+1];
                }
            }

            gMag[y][x] = hypot(accumX, accumY);
            gAngle[y][x] = atan2(accumY, accumX);
        }
    }
}

vector<vector<double>> copyMag(vector<vector<double>>& inp) {
    vector<vector<double>> out;
    for(int i=0; i<inp.size(); i++) {
        out.push_back(vector<double>());
        for(int j=0; j<inp[i].size(); j++) {
            out[i].push_back(inp[i][j]);
        }
    }
    return out;
}

class Point {
    private:
        int x, y;
    public:
        Point(int nx, int ny) {
            x = nx;
            y = ny;
        }

        int getX() { return x; }
        int getY() { return y; }
};

bool in_bounds(int y, int x) {
    return y > 0 && y < gMag.size() && x > 0 && x < gMag[y].size();
}


bool epsilonEquals(double one, double two) {
    return abs(one - two) < 0.000000001;
}


double MAX_MAG = -1;

    // normalize gradient to [0,1]
void normalizeGMag() {
    for(int y=1; y<height-1; y++) {
        for(int x=1; x<width-1; x++) {
            if(gMag[y][x] > MAX_MAG) {
                MAX_MAG = gMag[y][x];
            }
        }
    }

    for(int i=0; i<gMag.size(); i++) {
        for(int j=0; j<gMag[i].size(); j++) {
            gMag[i][j] /= MAX_MAG;
        }
    }
}

void hysteresis(Point strongEdge, vector<vector<double>>& mag) {
    for(int y=-3+strongEdge.getY(); y<3+strongEdge.getY(); y++) {
        for(int x=-3+strongEdge.getX(); x<3+strongEdge.getX(); x++) {
            if(in_bounds(y, x)) {
                if(mag[y][x] > 0 && mag[y][x] < 1) {
                    mag[y][x] = 1;
                    hysteresis(Point(x,y), mag);
                }
            }
        }
    }
}

vector<vector<double>> doubleThresholdAndHysteresis(vector<vector<double>>& mag, double weak, double high) {
    double highThreshold = high / MAX_MAG;
    double weakThreshold = weak / MAX_MAG;

    vector<Point> strongEdges;
    vector<Point> weakEdges;

    // double thresholding
    for(int y=1; y<height-1; y++) {
        for(int x=1; x<width-1; x++) {
            if(mag[y][x] > highThreshold) {
                mag[y][x] = 1;
                strongEdges.push_back(Point(x,y));
            } else if(gMag[y][x] < weakThreshold) {
                mag[y][x] = 0;
            } else {
                weakEdges.push_back(Point(x,y));
            }
        }
    }

    // hysteresis
    for(int i=0; i<strongEdges.size(); i++) {
        hysteresis(strongEdges[i], mag);
    }

    return mag;
}

vector<vector<Pixel>> getFinalPPM(vector<vector<double>> mag) {
    vector<vector<Pixel>> ret;
    for(int i=0; i<mag.size(); i++) {
        ret.push_back(vector<Pixel>());
        for(int j=0; j<mag[i].size(); j++) {
            int value = mag[i][j] * 255;
            if(value > 255/2) value = 255; else value = 0;
            ret[i].push_back(Pixel(value, value, value));
        }
    }
    return ret;
}

bool inAngleThresh(double angle, double lower, double upper) {
    return angle >= lower && angle < upper;
}

void nonMaxSupression(vector<vector<double>>& mag) {
    for(int y=1; y<height-1; y++) {
        for(int x=1; x<width-1; x++) {
            double magnitude = mag[y][x];
            double pi = 3.14159265;
            double angle = gAngle[y][x] * (180.0 / pi);

            // only need to cover [0, pi] this way
            if(angle < 0) angle += 180.0;

            int ax, ay;

            // horizontal ine
            if(inAngleThresh(angle, 0, 22.5) || (157.5 <= angle && angle <= 180)) {
                ay = 0;
                ax = 1;
            } else if(inAngleThresh(angle, 22.5, 67.5)) {
                ay = 1;
                ax = 1;
            } else if(inAngleThresh(angle, 67.5, 112.5)) {
                ay = 1;
                ax = 0;
            } else if(inAngleThresh(angle, 112.5, 157.5)) {
                ay = 1;
                ax = -1;
            }

            double a = mag[y+ay][x+ax];
            double b = mag[y-ay][x-ax];

            if(magnitude < a || magnitude < b) {
                mag[y][x] = 0;
            }
        }
    }
}


void simpleSingleThresh(double thresh, vector<vector<double>>& mag) {
    thresh /= MAX_MAG;
    for(int y=1; y<height-1; y++) {
        for(int x=1; x<width-1; x++) {
            if(mag[y][x] < thresh) {
                mag[y][x] = 0;
            } else {
                mag[y][x] = 1;
            }
        }
    }
}

void part3() {
    init();
    // image g 
    grayscale();
    sobel();
    normalizeGMag();
    vector<vector<double>> image1Mag = copyMag(gMag);
    vector<vector<double>> image2Mag = copyMag(gMag);
    vector<vector<double>> imageFMag = copyMag(gMag);
    
    // image 1
    doubleThresholdAndHysteresis(image1Mag, inputWeakThreshold, inputHighThreshold);
    outputImage("image1.ppm", getFinalPPM(image1Mag));

    // image 2
    nonMaxSupression(image2Mag);
    simpleSingleThresh(65, image2Mag);
    outputImage("image2.ppm", getFinalPPM(image2Mag));

    // image f
    nonMaxSupression(imageFMag);
    doubleThresholdAndHysteresis(imageFMag, inputWeakThreshold, inputHighThreshold);
    outputImage("imagef.ppm", getFinalPPM(imageFMag));
}

void part2() {
    init();
    grayscale();
    sobel();
    normalizeGMag();
    vector<vector<double>> imageMag = copyMag(gMag);
    vector<vector<double>> image1Mag = copyMag(gMag);

    // imagem
    simpleSingleThresh(65, imageMag);
    outputImage("imagem.ppm", getFinalPPM(imageMag));

    // image1
    doubleThresholdAndHysteresis(image1Mag, inputWeakThreshold, inputHighThreshold);
    outputImage("image1.ppm", getFinalPPM(image1Mag));
}

void part1() {
    init();
    grayscale();
    sobel();
    simpleSingleThresh(65, gMag);
    outputImage("imagem.ppm", getFinalPPM(gMag));
}


void readArgs(int argc, char** argv) {
    for(int i=0; i<argc; i++) {
        string arg = string(argv[i]);
        if(arg == "-L") {
            inputWeakThreshold = stoi(argv[i+1]);
        } else if(arg == "-H") {
            inputHighThreshold = stoi(argv[i+1]);
        } else if(arg == "-F") {
            filename = argv[i+1];
        }
    }
}

int main(int argc, char** argv) {
    readArgs(argc, argv);
    // part1();
    // part2();
    part3();
}

