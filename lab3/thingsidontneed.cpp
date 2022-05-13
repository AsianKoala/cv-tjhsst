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
			// compare using int cause when we actually draw we will cast point to int (because pixels on a string have integer locations)
            vector<Point>::iterator it = points.begin();
			while(it != points.end()) {
				if(!in_bounds((*it)) || (((int)(*it).getX() == (int)(*next(it)).getX()) && ((int)(*it).getY() == (int)(*next(it)).getY()))) {
					it = points.erase(it);
				} else {
					it++;
				}
			}
		}

		void add_pixels(int amt, int& row) {
			int i = 0;
			while(i<amt) {
				// file<<"255 255 255 ";
                file<<Color::white.tostring()<<" ";
				i++;
			}
			file<<endl;
			row++;
		}

        // used for swapping intenrals of points when we need to swap for drawing steep/shallow lines with bresenham
        void swap(Point& a, Point& b) {
            Point temp = Point(b.getX(), b.getY());
            b = Point(a.getX(), a.getY());
            a = temp;
        }

	public:
		VecPPM(string filename, int x, int y) {
			file.open(filename);
			file<<"P3 "<<x<<" "<<y<<" 255"<<endl;
			dimx = x;
			dimy = y;
		}

		void build() {
			int row = 0;
			int col = 0;


			setupvec();
			// max col should be dimx*6+1
			// max row should be dimy+1 (cause header)
			for(unsigned int i=0; i<points.size(); i++, col++) {
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
					// file<<"255 255 255 ";
                    file<<Color::white.tostring()<<" ";
					col++;
				}

                // file<<points[i].getColor().getR()<<" "<<points[i].getColor().getG()<<" "<<points[i].getColor().getB()<<" ";
                file<<points[i].getColor().tostring()<<" ";
			}

			add_pixels(dimx-col, row);
			while(row < dimy) {
				add_pixels(dimx, row);
			}

			file.close();
		}

        void drawshallowline(Point start, Point end, Color color) {
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
                points.push_back(Point(x, y, color));
                if(e >= 0) {
                    y += increment;
                    e -= dx;
                }
                e+=dy;
            }
        }

        void drawsteepline(Point start, Point end, Color color) {
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
                points.push_back(Point(x, y, color));
                if(e >= 0) {
                    x += increment;
                    e -= dy;
                }
                e += dx;
            }
        }

        // handles drawing of any type of line
        void drawline(Point start, Point end, Color color) {
            if(abs(end.getY() - start.getY()) < abs((end.getX() - start.getX()))) {
                if(start.getX() > end.getX()) {
                    swap(start,end);
                }
                drawshallowline(start, end, color);
            } else {
                if(start.getY() > end.getY()) {
                    swap(start,end);
                }
                drawsteepline(start, end, color);
            }
        }

        void drawcircle(Point c, int c_r, Color color) {
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
                
                points.push_back(Point((int)c.getX()+x, (int)c.getY()+y, color));
                points.push_back(Point((int)c.getX()-x, (int)c.getY()+y, color));
                points.push_back(Point((int)c.getX()+x, (int)c.getY()-y, color));
                points.push_back(Point((int)c.getX()-x, (int)c.getY()-y, color));
                points.push_back(Point((int)c.getX()+y, (int)c.getY()+x, color));
                points.push_back(Point((int)c.getX()-y, (int)c.getY()+x, color));
                points.push_back(Point((int)c.getX()+y, (int)c.getY()-x, color));
                points.push_back(Point((int)c.getX()-y, (int)c.getY()-x, color));
                y2_new -= (2 * x) - 3;
            }
        }
};