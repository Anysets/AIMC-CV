#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

//��ͬ����ϵͳѡ��ͬ��·��
// const string ROAD = "E:\\AR\\image\\����ͼ.BMP";
#define ROADBMP  "E:\\AR\\image\\����ͼ.BMP"
//#define ROADBMP  "/run/media/anysets/Files/Files/OpenCV++/test.BMP"

//��ȡͼƬ
Mat img = imread(ROADBMP);


int main()
{
    imshow("img", img);
    waitKey();
    destroyAllWindows();
    
    return 0;
}