#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include<iostream>
using namespace std;
using namespace cv;
//定义几种颜色
#define BLACK	Scalar(0, 0, 255)
#define WHITE	Scalar(255, 0, 0)
#define RED		Scalar(0, 0, 255)
#define RED1	Scalar(0, 17, 226)
#define BLUE	Scalar(255, 0, 0)
#define GREEN	Scalar(0, 255, 0)
#define YELLOW	Scalar(0, 255, 255)
#define BLUE1	Scalar(255, 255, 0)
#define BLUE2	Scalar(225, 169, 36)
#define BLUE3	Scalar(175, 165, 1)
#define PINK	Scalar(255, 0, 255)
#define PINK1	Scalar(67, 13, 90)
#define PINK2	Scalar(100, 13, 244)
#define PINK3	Scalar(101, 67, 256)
#define COLOR1	Scalar(154, 157, 252)
#define draw_point(x,y,color,map)	circle(map, Point(x, y), 0, color, -1);  // 画半径为0的圆(画点）
#define draw_circle(x,y,r,color,map)	circle(map, Point(x, y), r, color, 1);  // 画半径为r的圆