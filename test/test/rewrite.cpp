#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

//不同操作系统选择不同的路径
// const string ROAD = "E:\\AR\\image\\测试图.BMP";
#define ROADBMP  "E:\\AR\\image\\测试图.BMP"
//#define ROADBMP  "/run/media/anysets/Files/Files/OpenCV++/test.BMP"

//读取图片
Mat img = imread(ROADBMP);


int main()
{
    imshow("img", img);
    waitKey();
    destroyAllWindows();
    
    return 0;
}