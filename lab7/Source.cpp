#include <iostream>
#include <opencv2/opencv.hpp>
#include <fstream>

using namespace std;
using namespace cv;

Mat image;
Mat workingImage;
Mat grayImage;
Mat cannyImage;
Mat hsvImage;
vector<Vec3f> circles;
Scalar SCALAR_RED = Scalar(0, 0, 255);
Scalar SCALAR_BLUE = Scalar(255, 255, 0);
Scalar SCALAR_PURPLE = Scalar(128, 0, 128);
Scalar SCALAR_GREEN = Scalar(0, 255, 0);
Scalar SCALAR_YELLOW = Scalar(0, 255, 255);
int pennyCount = 0;
int nickelCount = 0;
int dimeCount = 0;
int quarterCount = 0;
int silverDollarCount = 0;
double dimeThreshold = 80;
double pennyThreshold = 100;
double nickelThreshold = 110;
double quarterThreshold = 125;
int amountOfCoins = 0;
double value = 0;
ofstream results;
string filename = "coins.jpg";
int cannyThreshold = 90;
int minDist = 70;
int lowRadius = 70;
int highRadius = 150;
int accumThreshold = 160;


void setupAndWriteImages() {	
	image = imread(filename);
	GaussianBlur(image, workingImage, Size(5, 5), 20);
	cvtColor(workingImage, hsvImage, COLOR_BGR2HSV);
	cvtColor(workingImage, workingImage, COLOR_BGR2GRAY);
	grayImage = workingImage.clone();
	Canny(grayImage, cannyImage, 30, 60);
	results.open("results.txt");
}


void getResult(Mat workingImage) {
	string pennyString = "pennies: " + to_string(pennyCount);
	string nickelString = "nickels: " + to_string(nickelCount);
	string dimeString = "dimes: " + to_string(dimeCount);
	string quarterString = "quarters: " + to_string(quarterCount);
	string silverDollarString = "silver dollars: " + to_string(silverDollarCount);
	string valueString = "total value: $" + to_string(value);

	results << pennyString << endl << nickelString << endl << dimeString << endl << quarterString << endl << silverDollarString << valueString << endl;
	results.close();

	putText(workingImage, pennyString, Point(100, 100), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);
	putText(workingImage, nickelString, Point(100, 150), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);
	putText(workingImage, dimeString, Point(100, 200), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);
	putText(workingImage, quarterString, Point(100, 250), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);
	putText(workingImage, silverDollarString, Point(100, 250), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);
	putText(workingImage, valueString, Point(100, 300), FONT_HERSHEY_SIMPLEX, 2, Scalar(0, 0, 0), 5);


	imwrite("imageg.jpg", grayImage);
	imwrite("imagef.jpg", cannyImage);
	imwrite("results.png", workingImage);

	namedWindow("window", 0);
	double ratio = ((double)image.rows) / ((double)image.cols);
	int x = 1000;
	int y = x * ratio;
	imshow("window", workingImage);
	resizeWindow("window", x, y);
	waitKey(0);
}

void decideCoin(Point center, double radius) {
	Scalar drawColor;
	if (radius < dimeThreshold) {
		dimeCount++;
		value += 0.10;
		drawColor = SCALAR_BLUE;
	}
	else if (radius < pennyThreshold) {
		pennyCount++;
		value += 0.01;
		drawColor = SCALAR_RED;
	}
	else if (radius < nickelThreshold) {
		nickelCount++;
		value += 0.05;
		drawColor = SCALAR_PURPLE;
	}
	else if(radius<quarterThreshold) {
		quarterCount++;
		value += 0.25;
		drawColor = SCALAR_GREEN;
	}
	else {
		silverDollarCount++;
		value += 1.0;
		drawColor = SCALAR_YELLOW;
	}

	circle(image, center, 3, drawColor, 3, LINE_AA);
	for (int j = radius; j < radius + 4; j++) {
		circle(image, center, j, drawColor, 3, LINE_AA);
	}
}
void work() {
	HoughCircles(workingImage, circles, HOUGH_GRADIENT, 2, minDist, cannyThreshold, accumThreshold, lowRadius, highRadius);
	amountOfCoins = circles.size();
	vector<Vec3i> pennyCircles;
	for (int i = 0; i < circles.size(); i++) {
		Vec3i c = circles[i];
		Point center = Point(c[0], c[1]);
		int radius = c[2];

		int brownHLow = 10;
		int brownHHigh = 16;
		int maxPennyThresh = 100;
		Vec3b hsv = hsvImage.at<Vec3b>(center);
		cout << hsv << endl;

		if (hsv[0] > brownHLow && hsv[0] < brownHHigh && radius < maxPennyThresh) {
			pennyCircles.push_back(c);
		}
	}

	double radsum = 0;
	for (int i = 0; i < pennyCircles.size(); i++) {
		radsum += pennyCircles[i][2];
	}
	double avgRad = radsum / pennyCircles.size();
	double thresholdRatio = avgRad / pennyThreshold;
	dimeThreshold *= thresholdRatio;
	pennyThreshold *= thresholdRatio;
	nickelThreshold *= thresholdRatio;
	quarterThreshold *= thresholdRatio;

	for (int i = 0; i < circles.size(); i++) {
		Vec3i c = circles[i];
		Point center = Point(c[0], c[1]);
		double rad = c[2];
		decideCoin(center, rad);
	}
}

void part1() {
	setupAndWriteImages();
	work();
	getResult(image);
}

int main() {
	part1();
}