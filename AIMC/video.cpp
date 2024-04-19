#include "image_process.h"

int main()
{
    VideoCapture capVideo("/home/anysets/Documents/Video/1.mp4");
    if (!capVideo.isOpened())
    {
        cout << "打开视频失败" << endl;
        return -1;
    }
    //创建每帧的文件
    Mat frame;
    while(true)
    {
        capVideo>>frame;//等价于cap.read(frame);
        if(frame.empty())//如果某帧为空则退出循环
            break;
 
        imshow("video", frame);
        waitKey(25);//每帧延时25毫秒
    }
    capVideo.release();//释放资


    return 0;
}