#define _CRT_SECURE_NO_WARNINGS
#include "image_process.h"
#include <fstream>
#include <vector>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/video.hpp>
#include <opencv2/imgproc/imgproc.hpp>
/*
这里是起点（0.0）************  col  *——>*************x值最大
************************************************************
************************************************************
************************************************************
************************************************************
******************假如这是一副图像*************************
row  *****************************************************
***********************************************************
***********************************************************
***********************************************************
***********************************************************
***********************************************************
y值最大*******************************************(94.60)
*/
#define IMAGE_H 80
#define IMAGE_W 188
#define image_h 40
#define image_w 94
#define start_line image_h - 1
uint8_t image[IMAGE_H][IMAGE_W];    //原始图像
uint8_t bin_image[image_h][image_w]; //最后使用的二值化图像
void dataTrans(const Mat& image, uint8_t imageArray[][IMAGE_W]) {
	Mat grayImage;
	cvtColor(image, grayImage, COLOR_BGR2GRAY); // 将彩色图像转换为灰度图像

	for (int i = 0; i < IMAGE_H; ++i) {
		for (int j = 0; j < IMAGE_W; ++j) {
			imageArray[i][j] = grayImage.at<uchar>(i, j);
		}
	}
}
#include <stdio.h>
void write_to_csv(uint8_t bin_image[][image_w]) {
    FILE* fp = fopen("E:\\AR\\image\\output.csv", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return;
    }

    // 将二值化后的数组写入CSV文件
    for (int i = 0; i < image_h; ++i) {
        for (int j = 0; j < image_w; ++j) {
            fprintf(fp, "%d,", bin_image[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
//大津法求阈值
uint8_t otsuThreshold(uint8_t image[][IMAGE_W]) {
    #define GrayScale 256
    int Pixel_Max = 0;
    int Pixel_Min = 255;
    uint16_t width = IMAGE_W;
    uint16_t height = IMAGE_H;
    int pixelCount[GrayScale];
    float pixelPro[GrayScale];
    int i, j, pixelSum = width * height / 4;
    uint8_t threshold = 0;

    for (i = 0; i < GrayScale; i++) {
        pixelCount[i] = 0;
        pixelPro[i] = 0;
    }

    uint32_t gray_sum = 0;
    // 统计灰度级中每个像素在整幅图像中的个数
    for (i = 0; i < height; i += 2) {
        for (j = 0; j < width; j += 2) {
            pixelCount[(int)image[i][j]]++;  // 将当前的点的像素值作为计数数组的下标
            gray_sum += (int)image[i][j];       // 灰度值总和
            if (image[i][j] > Pixel_Max)   Pixel_Max = image[i][j];
            if (image[i][j] < Pixel_Min)   Pixel_Min = image[i][j];
        }
    }

    // 计算每个像素值的点在整幅图像中的比例
    for (i = Pixel_Min; i < Pixel_Max; i++) {
        pixelPro[i] = (float)pixelCount[i] / pixelSum;
    }

    // 遍历灰度级[0,255]
    float w0, w1, u0tmp, u1tmp, u0, u1, u, deltaTmp, deltaMax = 0;
    w0 = w1 = u0tmp = u1tmp = u0 = u1 = u = deltaTmp = 0;
    for (j = Pixel_Min; j < Pixel_Max; j++) {
        w0 += pixelPro[j];  // 背景部分每个灰度值的像素点所占比例之和   即背景部分的比例
        u0tmp += j * pixelPro[j];  // 背景部分 每个灰度值的点的比例 *灰度值
        w1 = 1 - w0;
        u1tmp = gray_sum / pixelSum - u0tmp;
        u0 = u0tmp / w0;              // 背景平均灰度
        u1 = u1tmp / w1;              // 前景平均灰度
        u = u0tmp + u1tmp;            // 全局平均灰度
        deltaTmp = (float)(w0 * w1 * (u0 - u1) * (u0 - u1));
        if (deltaTmp > deltaMax) {
            deltaMax = deltaTmp;
            threshold = j;
        }
        if (deltaTmp < deltaMax) {
            break;
        }
    }
    //固定阈值
//if (threshold > 90 && threshold < 130)
//    last_threshold = threshold;
//else
//    threshold = last_threshold;

    return threshold;
}
// 大津法二值化
uint16_t white_point = 0;
void turn_to_bin(void) {
    uint8_t image_threshold = otsuThreshold(image);
    printf("阈值=%d\n", image_threshold);
    white_point = 0;    //先清零
    // 遍历原始图像数组，根据动态阈值进行二值化
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // 如果原始图像灰度值大于阈值，设置为白色像素
                white_point++;
            }
            else {
                image[i][j] = 0; // 如果原始图像灰度值小于等于阈值，设置为黑色像素
            }
        }
    }
    white_point /= 2;
    printf("white_point = %d\n", white_point);
}
//图像压缩一倍
void compressimage(void){
    int i, j, row, line;
    const float div_h = IMAGE_H / image_h, div_w = IMAGE_W / image_w;
    for (i = 0; i < image_h; i++){
        row = i * div_h + 0.5;
        for (j = 0; j < image_w; j++){
            line = j * div_w + 0.5;
            bin_image[i][j] = image[row][line];
        }
    }
}
//滤波
/***滤波***/
#define threshold_max	(255 * 5) // 此参数可根据自己的需求调节
#define threshold_min	(255 * 2) // 此参数可根据自己的需求调节
void image_filter(void) {
    uint16_t i, j;
    uint32_t num = 0;

    for (i = 1; i < image_h - 1; i++)
    {
        for (j = 1; j < (image_w - 1); j++)
        {
            // 统计八个方向的像素值
            num = bin_image[i - 1][j - 1] + bin_image[i - 1][j] + bin_image[i - 1][j + 1]
                + bin_image[i][j - 1] + bin_image[i][j + 1]
                + bin_image[i + 1][j - 1] + bin_image[i + 1][j] + bin_image[i + 1][j + 1];

            // 对像素进行处理
            if (num >= threshold_max && bin_image[i][j] == 0)
            {
                bin_image[i][j] = 255; // 白
            }
            if (num <= threshold_min && bin_image[i][j] == 255)
            {
                bin_image[i][j] = 0; // 黑
            }
        }
    }
}
//画个黑框
void image_draw_rectan(uint8_t image[image_h][image_w]) {
    for (uint8_t i = 0; i < image_h; i++) { //把第一列和最后一列画黑框
        image[i][0] = 0;
        image[i][1] = 0;
        image[i][image_w - 1] = 0;
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //把图像最上面画黑框
        image[0][i] = 0;
        //image[1][i] = 0;
    }
}
uint8_t start_left[2];    //row col
uint8_t start_right[2];   //row col
uint8_t get_start_point(uint8_t start_row) {
    //从中间往右边搜起点
    for (uint8_t col = image_w / 2; col < image_w; col++) {
        if (bin_image[start_row][col] == 255 && bin_image[start_row][col + 1] == 0
            && bin_image[start_row][col + 2] == 0) {
            start_right[0] = start_row;
            start_right[1] = col;
            break;
        }
    }
    //从中间往左边搜起点
    for (uint8_t col = image_w / 2; col > 0; col--) {
        if (bin_image[start_row][col] == 255 && bin_image[start_row][col - 1] == 0
            && bin_image[start_row][col - 2] == 0) {
            start_left[0] = start_row;
            start_left[1] = col;
            break;
        }
    }
    return 0;
}
//八邻域爬线
#define search_amount 100   //允许最多搜索到的点
uint8_t L_search_count, R_search_count; //搜索到的点
uint8_t left_points[search_amount][2],right_points[search_amount][2];    //row col
//{-1,-1},{0,-1},{+1,-1},
//{-1, 0},	     {+1, 0},
//{-1,+1},{0,+1},{+1,+1},
//这个是顺时针
int8_t L_dir_row[8] = {-1,-1,-1,0,1,1,1,0};
int8_t L_dir_col[8] = {-1,0,1,1,1,0,-1,-1};
//{-1,-1},{0,-1},{+1,-1},
//{-1, 0},	     {+1, 0},
//{-1,+1},{0,+1},{+1,+1},
//这个是逆时针
int8_t R_dir_row[8] = { -1,-1,-1,0,1,1,1,0 };
int8_t R_dir_col[8] = { 1,0,-1,-1,-1,0,1,1 };
uint8_t dir_left[search_amount], dir_right[search_amount];
void search(void) {
    L_search_count = R_search_count = 0;    //搜到的点数清零
    left_points[0][0] = start_left[0];  //y row
    left_points[0][1] = start_left[1];  //X col
    int8_t L_curr_row = start_left[0];   //初始化row坐标
    int8_t L_curr_col = start_left[1];  //初始化col坐标
    for(int j = 1; j < search_amount; j++){
        for (int i = 0; i < 8; i++) {
            if ((bin_image[L_curr_row + L_dir_row[i]][L_curr_col + L_dir_col[i]] == 0)
                && (bin_image[L_curr_row + L_dir_row[(i + 1) & 7]][L_curr_col + L_dir_col[(i + 1) & 7]] == 255)){  
                //由黑到白
                L_search_count++;
                L_curr_row += L_dir_row[(i + 1) & 7];
                L_curr_col += L_dir_col[(i + 1) & 7]; 
                left_points[L_search_count][0] = L_curr_row;   //存储白点
                left_points[L_search_count][1] = L_curr_col;   //存储白点 
                dir_left[L_search_count] = i;
                printf("第%d次，L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
                break;
            }
        }
    }
    right_points[0][0] = start_right[0];  //y row
    right_points[0][1] = start_right[1];  //X col
    int8_t R_curr_row = start_right[0];   //初始化row坐标
    int8_t R_curr_col = start_right[1];  //初始化col坐标
    for (int j = 1; j < search_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[R_curr_row + R_dir_row[i]][R_curr_col + R_dir_col[i]] == 0)
                && (bin_image[R_curr_row + R_dir_row[(i + 1) & 7]][R_curr_col + R_dir_col[(i + 1) & 7]] == 255)) {
                R_search_count++;
                R_curr_row += R_dir_row[(i + 1) & 7];
                R_curr_col += R_dir_col[(i + 1) & 7];
                right_points[R_search_count][0] = R_curr_row;   //存储白点
                right_points[R_search_count][1] = R_curr_col;   //存储白点
                dir_right[R_search_count] = i;
                printf("第%d次，R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
//提取边线
uint8_t left_line[image_h][2];  //左线数组 row col
uint8_t right_line[image_h][2];    //右线数组 row col
uint8_t center_line[image_h][2];   //中线数组 row col
void get_left(uint8_t total_L){
    uint8_t height = start_line;
    //初始化左线数组
    for (uint8_t i = image_h - 1; i > 0; i--) {
        left_line[i][0] = i;    //row
        left_line[i][1] = 0;    //col
    }
    //提取边线,一行只留第一个
    for (uint8_t i = 0; i < total_L; i++) { //遍历所有的点
        if (height == left_points[i][0]) {  
            left_line[height][1] = left_points[i][1];    //col
            height--;
            if (height < 0) break;
            continue;
        }
    }
}
void get_right(uint8_t total_L){
    uint8_t height = start_line;
    //初始化右线数组
    for (uint8_t i = image_h - 1; i > 0; i--) {
        right_line[i][0] = i;    //row
        right_line[i][1] = image_w - 1;    //col
    }
    //提取边线,一行只留第一个
    for (uint8_t i = 0; i < total_L; i++) { //遍历所有的点
        if (height == right_points[i][0]) {
            right_line[height][1] = right_points[i][1];    //col
            height--;
            if (height < 0) break;
            continue;
        }
    }
}
//最小二乘法
float Slope_Calculate(uint8_t begin, uint8_t end, uint8_t* border)
{
    float xsum = 0, ysum = 0, xysum = 0, x2sum = 0;
    int i = 0;
    float result = 0;
    static float resultlast;

    for (i = begin; i < end; i++)
    {
        xsum += i;
        ysum += border[i];
        xysum += i * (border[i]);
        x2sum += i * i;

    }
    if ((end - begin) * x2sum - xsum * xsum) //判断除数是否为零
    {
        result = ((end - begin) * xysum - xsum * ysum) / ((end - begin) * x2sum - xsum * xsum);
        resultlast = result;
    }
    else
    {
        result = resultlast;
    }
    return result;
}
//计算斜率和截距
void calculate_s_i(uint8_t start, uint8_t end, uint8_t* border, float* slope_rate, float* intercept)
{
    uint16_t i, num = 0;
    uint16_t xsum = 0, ysum = 0;
    float y_average, x_average;
    num = 0;
    xsum = 0;
    ysum = 0;
    y_average = 0;
    x_average = 0;
    for (i = start; i < end; i++)
    {
        xsum += i;
        ysum += border[i];
        num++;
    }
    //计算各个平均数
    if (num)
    {
        x_average = (float)(xsum / num);
        y_average = (float)(ysum / num);

    }
    /*计算斜率*/
    *slope_rate = Slope_Calculate(start, end, border);//斜率
    *intercept = y_average - (*slope_rate) * x_average;//截距
}
uint8_t element_start_line = 5; //从下面计数
uint8_t element_end_line = 40; //从下面计数
//void get_turning_point(void) {
//    //for (int i = element_start_line; i <= 60; i++) {
//    // /*   if(right_points[i][0] == )*/
//    //    printf("第%d次，第%d行\n",i, right_points[i][0]);
//    //}
//    for (int i = element_start_line; i < element_end_line; i++) {
//        if(dir_right[i] == 6 && dir_right[i + 1] == 6 && dir_right[i + 2] == 6)  //横线
//    }
//    
//}
void image_process(void)
{
    turn_to_bin();  //图像二值化
    compressimage(); //图像压缩
    image_filter();
    image_draw_rectan(bin_image);
    get_start_point(start_line);
    search();
    get_left(L_search_count);
    get_right(R_search_count);
    //get_turning_point();
}
int main(){
    //视频处理
    VideoCapture cap("/home/anysets/Documents/Video/1.mp4");
    // 检查视频是否成功打开
    if (!cap.isOpened()) {
        std::cerr << "Error opening video file" << std::endl;
        return -1;
    }
    Mat frame;
    int frame_number = 0;
    // 循环读取每一帧
    while (cap.read(frame)) {
         //保存每一帧为图片
        string filename = "frame_" + to_string(frame_number) + ".png";
        imwrite(filename, frame);
        waitKey(10);
        // 在这里对每一帧进行处理
        dataTrans(frame, image);
        image_process();
        // 创建一个彩色图像，与二值化数组的大小相同
       Mat Binimage_color(image_h, image_w, CV_8UC3, Scalar(0, 0, 0)); // 初始化为全黑色
       // 将二值化数组中的白色像素点绘制为指定颜色
       for (int i = 0; i < image_h; i++) {
           for (int j = 0; j < image_w; j++) {
               if (bin_image[i][j] == 255) { // 白色像素点
                   Binimage_color.at<Vec3b>(i, j) = Vec3b(255, 255, 255); // 设置为白色（BGR颜色）
               }
           }
       }
       for (int i = 0; i < image_h; i++) {
           draw_point(left_line[i][1], left_line[i][0],GREEN, Binimage_color);
           draw_point(right_line[i][1], right_line[i][0], GREEN, Binimage_color);
           //draw_point(left_points[i][1], left_points[i][0], PINK, Binimage_color);
           //draw_point(right_points[i][1], right_points[i][0], PINK, Binimage_color);
       }
       // 显示图像
       write_to_csv(bin_image);
       namedWindow("bin_image", WINDOW_NORMAL);
       imshow("bin_image", Binimage_color);
       resizeWindow("bin_imagel", 360, 240);
        namedWindow("output", WINDOW_NORMAL);
        imshow("output", frame);
        resizeWindow("output", 360, 240);
        // 等待一段时间，按下ESC键退出
        if (waitKey(3) == 27) {
            break;
        }
        frame_number++;
    }

 //   //图片处理
 //   const string ROAD = "E:\\AR\\image\\测试图.BMP";
 //   Mat roadImage = imread(ROAD);
 //   if (roadImage.empty())
 //   {
 //       cerr << "Failed to load image!" << endl;
 //       return -1;
 //   }
 //   dataTrans(roadImage, image);
 //   image_process();
 //   // 创建一个彩色图像，与二值化数组的大小相同
 //   Mat Binimage_color(image_h, image_w, CV_8UC3, Scalar(0, 0, 0)); // 初始化为全黑色
 //   // 将二值化数组中的白色像素点绘制为指定颜色
 //   for (int i = 0; i < image_h; i++) {
 //       for (int j = 0; j < image_w; j++) {
 //           if (bin_image[i][j] == 255) { // 白色像素点
 //               Binimage_color.at<Vec3b>(i, j) = Vec3b(255, 255, 255); // 设置为白色（BGR颜色）
 //           }
 //       }
 //   }
 //   //draw_point(start_left[1] + 1, start_left[0], PINK, Binimage_color);
 //   //draw_point(start_right[1] - 1, start_right[0], BLUE, Binimage_color);
 //   for (int i = 0; i < image_h; i++) {
 //       draw_point(left_line[i][1], left_line[i][0],GREEN, Binimage_color);
 //       draw_point(right_line[i][1], right_line[i][0], GREEN, Binimage_color);
 //       //draw_point(left_points[i][1], left_points[i][0], PINK, Binimage_color);
 //       //draw_point(right_points[i][1], right_points[i][0], PINK, Binimage_color);

 //   }
 //   // 显示图像
 //   write_to_csv(bin_image);
 //   namedWindow("二值化", WINDOW_NORMAL);
 //   imshow("二值化", Binimage_color);
 //   resizeWindow("二值化", 360, 240);
	//namedWindow("output", WINDOW_NORMAL);
	//imshow("output", roadImage);
	//resizeWindow("output", 360, 240);


	waitKey(0);
	destroyAllWindows();
	return 0;
}