#define _CRT_SECURE_NO_WARNINGS
#include "image_process.h"
#include <fstream>
#include <vector>
#include <opencv2/videoio/videoio.hpp>
#include <opencv2/video.hpp>
#include <opencv2/imgproc/imgproc.hpp>
/*
��������㣨0.0��************  col  *����>*************xֵ���
************************************************************
************************************************************
************************************************************
************************************************************
******************��������һ��ͼ��*************************
row  *****************************************************
***********************************************************
***********************************************************
***********************************************************
***********************************************************
***********************************************************
yֵ���*******************************************(94.60)
*/
#define IMAGE_H 80
#define IMAGE_W 188
#define image_h 40
#define image_w 94
#define start_line image_h - 1
uint8_t image[IMAGE_H][IMAGE_W];    //ԭʼͼ��
uint8_t bin_image[image_h][image_w]; //���ʹ�õĶ�ֵ��ͼ��
void dataTrans(const Mat& image, uint8_t imageArray[][IMAGE_W]) {
	Mat grayImage;
	cvtColor(image, grayImage, COLOR_BGR2GRAY); // ����ɫͼ��ת��Ϊ�Ҷ�ͼ��

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

    // ����ֵ���������д��CSV�ļ�
    for (int i = 0; i < image_h; ++i) {
        for (int j = 0; j < image_w; ++j) {
            fprintf(fp, "%d,", bin_image[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
//�������ֵ
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
    // ͳ�ƻҶȼ���ÿ������������ͼ���еĸ���
    for (i = 0; i < height; i += 2) {
        for (j = 0; j < width; j += 2) {
            pixelCount[(int)image[i][j]]++;  // ����ǰ�ĵ������ֵ��Ϊ����������±�
            gray_sum += (int)image[i][j];       // �Ҷ�ֵ�ܺ�
            if (image[i][j] > Pixel_Max)   Pixel_Max = image[i][j];
            if (image[i][j] < Pixel_Min)   Pixel_Min = image[i][j];
        }
    }

    // ����ÿ������ֵ�ĵ�������ͼ���еı���
    for (i = Pixel_Min; i < Pixel_Max; i++) {
        pixelPro[i] = (float)pixelCount[i] / pixelSum;
    }

    // �����Ҷȼ�[0,255]
    float w0, w1, u0tmp, u1tmp, u0, u1, u, deltaTmp, deltaMax = 0;
    w0 = w1 = u0tmp = u1tmp = u0 = u1 = u = deltaTmp = 0;
    for (j = Pixel_Min; j < Pixel_Max; j++) {
        w0 += pixelPro[j];  // ��������ÿ���Ҷ�ֵ�����ص���ռ����֮��   ���������ֵı���
        u0tmp += j * pixelPro[j];  // �������� ÿ���Ҷ�ֵ�ĵ�ı��� *�Ҷ�ֵ
        w1 = 1 - w0;
        u1tmp = gray_sum / pixelSum - u0tmp;
        u0 = u0tmp / w0;              // ����ƽ���Ҷ�
        u1 = u1tmp / w1;              // ǰ��ƽ���Ҷ�
        u = u0tmp + u1tmp;            // ȫ��ƽ���Ҷ�
        deltaTmp = (float)(w0 * w1 * (u0 - u1) * (u0 - u1));
        if (deltaTmp > deltaMax) {
            deltaMax = deltaTmp;
            threshold = j;
        }
        if (deltaTmp < deltaMax) {
            break;
        }
    }
    //�̶���ֵ
//if (threshold > 90 && threshold < 130)
//    last_threshold = threshold;
//else
//    threshold = last_threshold;

    return threshold;
}
// ��򷨶�ֵ��
uint16_t white_point = 0;
void turn_to_bin(void) {
    uint8_t image_threshold = otsuThreshold(image);
    printf("��ֵ=%d\n", image_threshold);
    white_point = 0;    //������
    // ����ԭʼͼ�����飬���ݶ�̬��ֵ���ж�ֵ��
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // ���ԭʼͼ��Ҷ�ֵ������ֵ������Ϊ��ɫ����
                white_point++;
            }
            else {
                image[i][j] = 0; // ���ԭʼͼ��Ҷ�ֵС�ڵ�����ֵ������Ϊ��ɫ����
            }
        }
    }
    white_point /= 2;
    printf("white_point = %d\n", white_point);
}
//ͼ��ѹ��һ��
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
//�˲�
/***�˲�***/
#define threshold_max	(255 * 5) // �˲����ɸ����Լ����������
#define threshold_min	(255 * 2) // �˲����ɸ����Լ����������
void image_filter(void) {
    uint16_t i, j;
    uint32_t num = 0;

    for (i = 1; i < image_h - 1; i++)
    {
        for (j = 1; j < (image_w - 1); j++)
        {
            // ͳ�ư˸����������ֵ
            num = bin_image[i - 1][j - 1] + bin_image[i - 1][j] + bin_image[i - 1][j + 1]
                + bin_image[i][j - 1] + bin_image[i][j + 1]
                + bin_image[i + 1][j - 1] + bin_image[i + 1][j] + bin_image[i + 1][j + 1];

            // �����ؽ��д���
            if (num >= threshold_max && bin_image[i][j] == 0)
            {
                bin_image[i][j] = 255; // ��
            }
            if (num <= threshold_min && bin_image[i][j] == 255)
            {
                bin_image[i][j] = 0; // ��
            }
        }
    }
}
//�����ڿ�
void image_draw_rectan(uint8_t image[image_h][image_w]) {
    for (uint8_t i = 0; i < image_h; i++) { //�ѵ�һ�к����һ�л��ڿ�
        image[i][0] = 0;
        image[i][1] = 0;
        image[i][image_w - 1] = 0;
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //��ͼ�������滭�ڿ�
        image[0][i] = 0;
        //image[1][i] = 0;
    }
}
uint8_t start_left[2];    //row col
uint8_t start_right[2];   //row col
uint8_t get_start_point(uint8_t start_row) {
    //���м����ұ������
    for (uint8_t col = image_w / 2; col < image_w; col++) {
        if (bin_image[start_row][col] == 255 && bin_image[start_row][col + 1] == 0
            && bin_image[start_row][col + 2] == 0) {
            start_right[0] = start_row;
            start_right[1] = col;
            break;
        }
    }
    //���м�����������
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
//����������
#define search_amount 100   //��������������ĵ�
uint8_t L_search_count, R_search_count; //�������ĵ�
uint8_t left_points[search_amount][2],right_points[search_amount][2];    //row col
//{-1,-1},{0,-1},{+1,-1},
//{-1, 0},	     {+1, 0},
//{-1,+1},{0,+1},{+1,+1},
//�����˳ʱ��
int8_t L_dir_row[8] = {-1,-1,-1,0,1,1,1,0};
int8_t L_dir_col[8] = {-1,0,1,1,1,0,-1,-1};
//{-1,-1},{0,-1},{+1,-1},
//{-1, 0},	     {+1, 0},
//{-1,+1},{0,+1},{+1,+1},
//�������ʱ��
int8_t R_dir_row[8] = { -1,-1,-1,0,1,1,1,0 };
int8_t R_dir_col[8] = { 1,0,-1,-1,-1,0,1,1 };
uint8_t dir_left[search_amount], dir_right[search_amount];
void search(void) {
    L_search_count = R_search_count = 0;    //�ѵ��ĵ�������
    left_points[0][0] = start_left[0];  //y row
    left_points[0][1] = start_left[1];  //X col
    int8_t L_curr_row = start_left[0];   //��ʼ��row����
    int8_t L_curr_col = start_left[1];  //��ʼ��col����
    for(int j = 1; j < search_amount; j++){
        for (int i = 0; i < 8; i++) {
            if ((bin_image[L_curr_row + L_dir_row[i]][L_curr_col + L_dir_col[i]] == 0)
                && (bin_image[L_curr_row + L_dir_row[(i + 1) & 7]][L_curr_col + L_dir_col[(i + 1) & 7]] == 255)){  
                //�ɺڵ���
                L_search_count++;
                L_curr_row += L_dir_row[(i + 1) & 7];
                L_curr_col += L_dir_col[(i + 1) & 7]; 
                left_points[L_search_count][0] = L_curr_row;   //�洢�׵�
                left_points[L_search_count][1] = L_curr_col;   //�洢�׵� 
                dir_left[L_search_count] = i;
                printf("��%d�Σ�L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
                break;
            }
        }
    }
    right_points[0][0] = start_right[0];  //y row
    right_points[0][1] = start_right[1];  //X col
    int8_t R_curr_row = start_right[0];   //��ʼ��row����
    int8_t R_curr_col = start_right[1];  //��ʼ��col����
    for (int j = 1; j < search_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[R_curr_row + R_dir_row[i]][R_curr_col + R_dir_col[i]] == 0)
                && (bin_image[R_curr_row + R_dir_row[(i + 1) & 7]][R_curr_col + R_dir_col[(i + 1) & 7]] == 255)) {
                R_search_count++;
                R_curr_row += R_dir_row[(i + 1) & 7];
                R_curr_col += R_dir_col[(i + 1) & 7];
                right_points[R_search_count][0] = R_curr_row;   //�洢�׵�
                right_points[R_search_count][1] = R_curr_col;   //�洢�׵�
                dir_right[R_search_count] = i;
                printf("��%d�Σ�R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
//��ȡ����
uint8_t left_line[image_h][2];  //�������� row col
uint8_t right_line[image_h][2];    //�������� row col
uint8_t center_line[image_h][2];   //�������� row col
void get_left(uint8_t total_L){
    uint8_t height = start_line;
    //��ʼ����������
    for (uint8_t i = image_h - 1; i > 0; i--) {
        left_line[i][0] = i;    //row
        left_line[i][1] = 0;    //col
    }
    //��ȡ����,һ��ֻ����һ��
    for (uint8_t i = 0; i < total_L; i++) { //�������еĵ�
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
    //��ʼ����������
    for (uint8_t i = image_h - 1; i > 0; i--) {
        right_line[i][0] = i;    //row
        right_line[i][1] = image_w - 1;    //col
    }
    //��ȡ����,һ��ֻ����һ��
    for (uint8_t i = 0; i < total_L; i++) { //�������еĵ�
        if (height == right_points[i][0]) {
            right_line[height][1] = right_points[i][1];    //col
            height--;
            if (height < 0) break;
            continue;
        }
    }
}
//��С���˷�
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
    if ((end - begin) * x2sum - xsum * xsum) //�жϳ����Ƿ�Ϊ��
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
//����б�ʺͽؾ�
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
    //�������ƽ����
    if (num)
    {
        x_average = (float)(xsum / num);
        y_average = (float)(ysum / num);

    }
    /*����б��*/
    *slope_rate = Slope_Calculate(start, end, border);//б��
    *intercept = y_average - (*slope_rate) * x_average;//�ؾ�
}
uint8_t element_start_line = 5; //���������
uint8_t element_end_line = 40; //���������
//void get_turning_point(void) {
//    //for (int i = element_start_line; i <= 60; i++) {
//    // /*   if(right_points[i][0] == )*/
//    //    printf("��%d�Σ���%d��\n",i, right_points[i][0]);
//    //}
//    for (int i = element_start_line; i < element_end_line; i++) {
//        if(dir_right[i] == 6 && dir_right[i + 1] == 6 && dir_right[i + 2] == 6)  //����
//    }
//    
//}
void image_process(void)
{
    turn_to_bin();  //ͼ���ֵ��
    compressimage(); //ͼ��ѹ��
    image_filter();
    image_draw_rectan(bin_image);
    get_start_point(start_line);
    search();
    get_left(L_search_count);
    get_right(R_search_count);
    //get_turning_point();
}
int main(){
    //��Ƶ����
    VideoCapture cap("/home/anysets/Documents/Video/1.mp4");
    // �����Ƶ�Ƿ�ɹ���
    if (!cap.isOpened()) {
        std::cerr << "Error opening video file" << std::endl;
        return -1;
    }
    Mat frame;
    int frame_number = 0;
    // ѭ����ȡÿһ֡
    while (cap.read(frame)) {
         //����ÿһ֡ΪͼƬ
        string filename = "frame_" + to_string(frame_number) + ".png";
        imwrite(filename, frame);
        waitKey(10);
        // �������ÿһ֡���д���
        dataTrans(frame, image);
        image_process();
        // ����һ����ɫͼ�����ֵ������Ĵ�С��ͬ
       Mat Binimage_color(image_h, image_w, CV_8UC3, Scalar(0, 0, 0)); // ��ʼ��Ϊȫ��ɫ
       // ����ֵ�������еİ�ɫ���ص����Ϊָ����ɫ
       for (int i = 0; i < image_h; i++) {
           for (int j = 0; j < image_w; j++) {
               if (bin_image[i][j] == 255) { // ��ɫ���ص�
                   Binimage_color.at<Vec3b>(i, j) = Vec3b(255, 255, 255); // ����Ϊ��ɫ��BGR��ɫ��
               }
           }
       }
       for (int i = 0; i < image_h; i++) {
           draw_point(left_line[i][1], left_line[i][0],GREEN, Binimage_color);
           draw_point(right_line[i][1], right_line[i][0], GREEN, Binimage_color);
           //draw_point(left_points[i][1], left_points[i][0], PINK, Binimage_color);
           //draw_point(right_points[i][1], right_points[i][0], PINK, Binimage_color);
       }
       // ��ʾͼ��
       write_to_csv(bin_image);
       namedWindow("bin_image", WINDOW_NORMAL);
       imshow("bin_image", Binimage_color);
       resizeWindow("bin_imagel", 360, 240);
        namedWindow("output", WINDOW_NORMAL);
        imshow("output", frame);
        resizeWindow("output", 360, 240);
        // �ȴ�һ��ʱ�䣬����ESC���˳�
        if (waitKey(3) == 27) {
            break;
        }
        frame_number++;
    }

 //   //ͼƬ����
 //   const string ROAD = "E:\\AR\\image\\����ͼ.BMP";
 //   Mat roadImage = imread(ROAD);
 //   if (roadImage.empty())
 //   {
 //       cerr << "Failed to load image!" << endl;
 //       return -1;
 //   }
 //   dataTrans(roadImage, image);
 //   image_process();
 //   // ����һ����ɫͼ�����ֵ������Ĵ�С��ͬ
 //   Mat Binimage_color(image_h, image_w, CV_8UC3, Scalar(0, 0, 0)); // ��ʼ��Ϊȫ��ɫ
 //   // ����ֵ�������еİ�ɫ���ص����Ϊָ����ɫ
 //   for (int i = 0; i < image_h; i++) {
 //       for (int j = 0; j < image_w; j++) {
 //           if (bin_image[i][j] == 255) { // ��ɫ���ص�
 //               Binimage_color.at<Vec3b>(i, j) = Vec3b(255, 255, 255); // ����Ϊ��ɫ��BGR��ɫ��
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
 //   // ��ʾͼ��
 //   write_to_csv(bin_image);
 //   namedWindow("��ֵ��", WINDOW_NORMAL);
 //   imshow("��ֵ��", Binimage_color);
 //   resizeWindow("��ֵ��", 360, 240);
	//namedWindow("output", WINDOW_NORMAL);
	//imshow("output", roadImage);
	//resizeWindow("output", 360, 240);


	waitKey(0);
	destroyAllWindows();
	return 0;
}