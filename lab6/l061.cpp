#include <bits/stdc++.h>

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

        string toString() {
            return to_string(r) + " " + to_string(g) + " " + to_string(b);
        }
};

string filename = "image.ppm";
vector<vector<Pixel>> image;
string ppmType;
int width, height, rgbSize;

void readppm() {
    // read header
    ifstream file(filename);
    file>>ppmType;
    string tempWidth, tempHeight;
    file>>tempWidth>>tempHeight;
    width = stoi(tempWidth);
    height = stoi(tempHeight);
    string tempRGBSize;
    file>>tempRGBSize;
    rgbSize = stoi(tempRGBSize);

    cout<<"header finished"<<endl;
    for(int h=0; h<height; h++) {
        cout<<"on height "<<h<<endl;
        image.push_back(vector<Pixel>());
        for(int w=0; w<width; w++) {
            int r,g,b; file>>r>>g>>b;
            image[h].push_back(Pixel(r,g,b));
        }
    }
}


void part1() {
    readppm();
}

int main() {
    part1();
}