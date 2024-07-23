#include "image_process.h"
#include <fstream>
#include <vector>
#include <math.h>
// #include <stdlib.h>
// 58 71
#define _CRT_SECURE_NO_WARNINGS
/*
��������㣨0.0��************  col  *����>*************xֵ���sss
************************************************************
************************************************************
************************************************************
************************************************************
**********************��������һ��ͼ��*************************
row  *******************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
yֵ���***********************************************(94.60)
*/
#define ROAD  "/home/anysets/Documents/AIMC/frames/frame_84.png"
#define VIDEO  "/home/anysets/Documents/AIMC/Video/82.mp4"
#define PVFLAG 1
// #define ROAD "E:\\Files\\OpenCV++\\test.BMP"
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 94
//ע�����border_length����̫�࣬�ᵼ������Խ��
#define border_length 60
int image_center = image_w / 2;
int search_start_set = 25;
int search_stop_set = 35;

//��ͼ�������ӡ��csv�ļ�--------------------------------------------------------
void printToCSV(uint8_t matrix[][image_w], const std::string& filename) {
	std::ofstream outputFile(filename);

	if (!outputFile.is_open()) {
		std::cerr << "Cannot open file!" << filename << std::endl;
		return;
	}

	for (int i = 0; i < image_h; ++i) {
		for (int j = 0; j < image_w; ++j) {
			outputFile << static_cast<int>(matrix[i][j]);
			if (j != image_w - 1)
				outputFile << ",";
		}
		outputFile << std::endl;
	}

	outputFile.close();
	// printf("print to csv successed\n");
}

uint8_t image[IMAGE_H][IMAGE_W];    //ԭʼͼ��
uint8_t bin_image[image_h][image_w]; //���ʹ�õĶ�ֵ��ͼ��
int isStartFlag = 1;
void dataTrans(const Mat& image, uint8_t imageArray[][IMAGE_W]) {
	Mat grayImage;
	cvtColor(image, grayImage, COLOR_BGR2GRAY); // ����ɫͼ��ת��Ϊ�Ҷ�ͼ��

	for (int i = 0; i < IMAGE_H; ++i) {
		for (int j = 0; j < IMAGE_W; ++j) {
			imageArray[i][j] = grayImage.at<uchar>(i, j);
		}
	}
}
#define GrayScale 256

#define BLOCK_SIZE 15
#define C 5

void adaptive_threshold(unsigned char image[][IMAGE_W], int width, int height, int block_size, int c) {
    int half_block = block_size / 2;
    int sum, count, x, y, i, j;

    // ����һ������������ԭʼͼ������
    unsigned char original[IMAGE_H][IMAGE_W];
    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            original[y][x] = image[y][x];
        }
    }

    for (y = 0; y < height; ++y) {
        for (x = 0; x < width; ++x) {
            sum = 0;
            count = 0;

            // ����ֲ���ֵ
            for (j = -half_block; j <= half_block; ++j) {
                for (i = -half_block; i <= half_block; ++i) {
                    int yy = y + j;
                    int xx = x + i;
                    if (yy >= 0 && yy < height && xx >= 0 && xx < width) {
                        sum += original[yy][xx];
                        count++;
                    }
                }
            }

            // ����ֲ���ֵ
            int mean = sum / count;
            int threshold = mean - c;

            // Ӧ����ֵ
            if (original[y][x] > threshold) {
                image[y][x] = 255;
            }
            else {
                image[y][x] = 0;
            }
        }
    }
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
void turn_to_bin(void) {
    uint8_t image_threshold = otsuThreshold(image);
    // uint8_t image_threshold = 176;
    // printf("OTSU=%d\n", image_threshold);

    // ����ԭʼͼ�����飬���ݶ�̬��ֵ���ж�ֵ��
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // ���ԭʼͼ��Ҷ�ֵ������ֵ������Ϊ��ɫ����
            }
            else {
                image[i][j] = 0; // ���ԭʼͼ��Ҷ�ֵС�ڵ�����ֵ������Ϊ��ɫ����
            }
        }
    }
}

//ͼ��ѹ��һ��
void image_compress(uint8_t image[IMAGE_H][IMAGE_W], uint8_t iamge_zip[image_h][image_w]) {
	for (uint8_t row = 0; row < IMAGE_H; row+=2) {
		for (uint8_t col = 0; col < IMAGE_W; col+=2) {
			iamge_zip[row/2][col/2] = image[row][col];
		}
	}
}
/***�˲�***/
#define threshold_max	(255 * 7) // �˲����ɸ����Լ����������
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
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //��ͼ�������滭�ڿ�
        image[0][i] = 0;
        image[1][i] = 0;
        image[image_h - 1][i] = 0;  //��ͼ�������滭һ�кڿ��Խ������Խ������
    }
}

//�����ڿ�
void image_draw_rectan(uint8_t image[image_h][image_w]) {
    for (uint8_t i = 0; i < image_h; i++) { //��ǰ���к�������л��ڿ�
        image[i][0] = 0;
        image[i][1] = 0;
        image[i][image_w - 1] = 0;
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //��ͼ�������滭�ڿ�
        image[0][i] = 0;
        image[1][i] = 0;
        image[image_h - 1][i] = 0;  //��ͼ�������滭һ�кڿ���ͼ�������Խ������
    }
}

//�ҵ���㣬����ײ������Ŀ�ʼѰ��
// 6.14 ���Ը�Ϊ���ǰ���
int start_point_l[2] = { 0 };
int start_point_r[2] = { 0 };
int start_line = image_h - 2;  //���һ�б������ˣ�����һ��
int getStartPoint()
{
    start_line = image_h - 2;
    image_center = image_w / 2;
    int left_flag = 0;
    int right_flag = 0;
    //��ʼ��
    start_point_l[0] = 0;
    start_point_l[1] = 0;
    start_point_r[0] = 0;
    start_point_r[1] = 0;
    // bin_image[start_line][image_center] == 255;
    if (bin_image[start_line][image_center] == 0)
    {
        image_center = image_w / 4;
        if (bin_image[start_line][image_center] == 0)
        {
            image_center = image_w / 4 * 3;
            if (bin_image[start_line][image_center] == 0)
            {
                isStartFlag = 0;
                return 0;
            }
        }
    }
    //���м�����Ѱ��
    for (int i = image_center; i > 1; i--)  //i>1��ԭ����0��1�����Ѿ���ȫ����0�ˣ�����i=0��1ʱ���������Խ��
    {
        if (bin_image[start_line][i] == 255 && bin_image[start_line][i - 1] == 0 && bin_image[start_line][i - 2] == 0)
        {
            start_point_l[0] = start_line; //��
            start_point_l[1] = i; //��
            left_flag = 1;
            break;
        }
    }
    //���м�����Ѱ��
    for (int i = image_center; i < image_w - 2; i++)  //i < image_w-2��ԭ��ͬ��
    {
        if (bin_image[start_line][i] == 255 && bin_image[start_line][i + 1] == 0 && bin_image[start_line][i + 2] == 0)
        {
            start_point_r[0] = start_line; // image_h
            start_point_r[1] = i;  // �˴���¼���ǰ�ɫ��
            right_flag = 1;
            break;
        }
    }
    if (left_flag && right_flag) { isStartFlag = 1;return 1; }
    else { isStartFlag = 1;return 1; }
}

//������rewrite
//��¼��߽������
int border_location_left[border_length][2];  //��¼�߽�λ��
int border_count_left = 0;
int growth_direction_left[border_length];  //��¼��������
int growth_count_left = 0;
//��¼�ұ߽������
int border_location_right[border_length][2];  //��¼�߽�λ��
int border_count_right = 0;
int growth_direction_right[border_length];  //��¼��������
int growth_count_right = 0;
void neighborSearch()
{
    memset(border_location_left, 0, sizeof(border_location_left));
    memset(border_location_right, 0, sizeof(border_location_right));

    int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
    int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
    border_count_left = 1;
    growth_count_left = 1;
    border_count_right = 1;
    growth_count_right = 1;
    //�����ӽ�ȥ
    border_location_left[0][0] = start_point_l[0];
    border_location_left[0][1] = start_point_l[1];
    border_location_right[0][0] = start_point_r[0];
    border_location_right[0][1] = start_point_r[1];
    for (int j = 0; j < border_length-1; j++)
    {
        //��߽粿��
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        // printf("%d\n", j);
        for (int i = 0; i < 8; i++)
        {
            // printf("ciallo:%d\n", bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]]);
            // printf("ciallo+1:%d\n", bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]]);
            if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
            {
                // printf("h=%d\n", bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]]);
                // printf("h=%d\n", bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]]);
                growth_direction_left[growth_count_left] = i;
                growth_count_left++;

                // printf("start_point_l[0] + neighbor_left[i+1][0]=%d\n", start_point_l[0] + neighbor_left[i+1][0]);
                // border_location_left[border_count_left][0] = start_point_l[0] + neighbor_left[i+1][0];
                // border_location_left[border_count_left][1] = start_point_l[1] + neighbor_left[i+1][1];
                // printf("%d : %d, %d\n", j, start_point_l[0] + neighbor_left[i+1][0], start_point_l[1] + neighbor_left[i+1][1]);
                // printf("%d:border_location_left[7]: %d, %d\n", j, border_location_left[7][0], border_location_left[7][1]);
                // printf("%d\n", border_count_left);
                // bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] = 151;
                start_point_l[0] = start_point_l[0] + neighbor_left[i][0];  // ����ǰ�׵�����Ϊ��һ�����
                start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                border_location_left[border_count_left][0] = start_point_l[0];
                border_location_left[border_count_left][1] = start_point_l[1];
                border_count_left++;
                // printf("New Left Start Point: %d, %d\n", start_point_l[0], start_point_l[1]);
                // printf("finished\n");
                // printf("%d\n", j);
                break;
            }
        }
        //�ұ߽粿��
        /*
            [0/8]  [1]  [2]

            [7]  [ ]  [3]

            [6]  [5]  [4]
        */
        for (int i = 0; i < 8; i++)
        {
            if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
            {
                growth_direction_right[growth_count_right] = i;
                growth_count_right++;

                border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                border_count_right++;
                break;
            }
        }
        //�����߽���ʱֹͣ����
        if ((j != 0 && j != border_length-1) &&
        ((border_location_left[j][0] == border_location_right[j][0] && border_location_left[j][1] == border_location_right[j][1])
        || (border_location_left[j+1][0] == border_location_right[j-1][0] && border_location_left[j+1][1] == border_location_right[j-1][1])
        || (border_location_left[j-1][0] == border_location_right[j+1][0] && border_location_left[j-1][1] == border_location_right[j+1][1])))
        {
            break;
        }
    }
    // ���ڷ�����������
    // printf("growth_count_left: %d", growth_count_left);
    // for (int i = 0; i< 200; i++)
    // {
    //  printf("%d\n", growth_direction_left[i]);
    // }

    // ���border_location_left
    // for (int i = 0; i < border_length; i++)
    // {
    //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
    // }
}

// �洢�Ż������߽磬ÿ�����ض��зݣ�ע��˴���ȡֵ��ΧΪ-127��127��������ѹ�����ͼ��Χ��
// ����һ��Ҫע�⣬�˴��Ƕ�ά���飬���ǵ�0��ֵ��image_w����1��ֵ�ǵ�ǰ����ԭ�߽����������ǰֵ����������image_h
int border_location_left_optimized[image_h][2] = { 0 };
int border_location_right_optimized[image_h][2] = { 0 };
void optimizeBorder()
{
    // ------------------------------------------------------------��߽�Ĳ���
    // ѭ���Ѿ���ȡ���ı߽磬ȥ������Ҫ�ĵ�
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //һ��ֻҪһ�������ˣ��е�������[0]����߽���Ҫ����ÿ�����Ҳ�ĵ�
        if (border_location_left[i][1] > border_location_left_optimized[border_location_left[i][0]][0])
        {
            border_location_left_optimized[border_location_left[i][0]][0] = border_location_left[i][1];
            border_location_left_optimized[border_location_left[i][0]][1] = i;
        }
    }
    // for (uint8_t i = 0; i < image_h; i++)
    // {
        // bin_image[i][border_location_left_optimized[i][0]] = 152;  // ���߽��ǳ���
        // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i][0]);
    // }
    // �ұ߽�Ĳ���
    // ѭ���Ѿ���ȡ���ı߽磬ȥ������Ҫ�ĵ�
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //һ��ֻҪһ�������ˣ��е�������[0]���ұ߽���Ҫ����ÿ�������ĵ�
        if (border_location_right_optimized[border_location_right[i][0]][0] == 0 || border_location_right[i][1] < border_location_right_optimized[border_location_right[i][0]][0])
        {
            border_location_right_optimized[border_location_right[i][0]][0] = border_location_right[i][1];
            border_location_right_optimized[border_location_right[i][0]][1] = i;
        }

    }
    // for (int i = 0; i < image_h; i++)
    // {
    //     // bin_image[i][border_location_right_optimized[i][0]] = 152;  // ���߽��ǳ���
    //     printf("border_optimized[%d],:%d\n", i, border_location_right_optimized[i][0]);
    // }
}

// ��ȡ����
// int32_t middle_line[border_length] = { 0 };
int middle_line[image_w] = { 0 };
void getMiddleLine()
{
    memset(middle_line, 0, sizeof(middle_line));
    for (int i = 0; i < image_h; i++)
    {
        middle_line[i] = (border_location_left_optimized[i][0] + border_location_right_optimized[i][0]) / 2;
    }
}

// ��ȡ�Ƕ�
// double angle_deg = 0;
// void getAngle(double_t opposite, uint8_t adjacent)
// {
//     // opposite = abs(opposite);
//     // printf("opposite=%d\n", opposite);
//     // printf("adjacent=%d\n", adjacent);
//     double angle_rad = atan2(opposite, adjacent);
//     angle_deg = angle_rad * (180.0 / 3.14);
//     // printf("angle_deg=%f\n", angle_deg);
// }

int highestPoint = 0;
void getHighestPoint()
{
    highestPoint = 0;

    for (int i = 0; i < image_h; i++)
    {
        // printf("middle_line:%d, %d\n", middle_line[i] ,(0 + image_w-1) / 2);
        // printf("border_location_left_optimized:%d %d\n", border_location_left_optimized[i][0], 0);
        // printf("border_location_right_optimized:%d %d\n", border_location_right_optimized[i][0], image_w-1);
        if ((middle_line[i] != (0 + image_h-1) / 2 && border_location_left_optimized[i][0] != 0 && border_location_right_optimized[i][0] != image_w-1))
        {
            // printf("Ciallo~ %d\n", i);
            highestPoint = i;
            break;
        }
    }
}

// ��ȡƫ��ֵ
// int32_t deviation_number;
// int32_t deviation_number_p = 0;  // previous
double deviation_number;
double deviation_number_p = 0;
int linear_flag = 0;

int search_start = 0;
int search_stop = 0;
void getDeviation()
{
    search_start = search_start_set;
    search_stop = search_start+10;
    int count = 0;
    deviation_number = 0;
    if (highestPoint > search_start)
    {
        if(highestPoint < search_stop)
        {
            search_start = highestPoint;
            // ��ʱ��ʩ���ǵøģ���Ҫ��ע������в�Ҫ����߽磡������������������������������������������������������������������������������������������������������������������������������
            search_stop = search_stop + (highestPoint - search_start);
        }
        else if(highestPoint < 50)
        {
            search_start = highestPoint;
            search_stop = 55;
        }
        else
        {
            search_start = 50;
            search_stop = 57;
        }
    }
    for (int i = search_start; i < search_stop; i++)
    {
        deviation_number += middle_line[i] - (image_w/2);
        count++;
    }
    deviation_number = deviation_number / count;
    if (deviation_number == -47.00)
    {
        deviation_number = deviation_number_p;
    }

    // getAngle(deviation_number, 15);

    // �ڴ��������Ĳ���
    // �����ֱ�ߣ���ƫ��ֵ��Ȩ��
    // if (linear_flag != 0)
    // {
    //     if (fabs(deviation_number_p) < fabs(deviation_number))
    //     {
    //         deviation_number = 0.5 * deviation_number_p + 0.5 * deviation_number;
    //     }
    // }
    deviation_number = deviation_number *0.6 + deviation_number_p * 0.4;
    // ��¼ƫ��ֵ
    deviation_number_p = deviation_number;
}


//Ѱ�ҹյ�
int diffPoint_left[2] = {0}; // ��¼�յ�
int diffPoint_right[2] = {0}; // ��¼�յ�
int findTurningPoint()
{
    memset(diffPoint_left, 0, sizeof(diffPoint_left));
    memset(diffPoint_right, 0, sizeof(diffPoint_right));
    int flag = 0;
    int left_flag = 0;
    int right_flag = 0;
    int diffPoint_temp[2] = {0};
    // ----------------------------------------------------------------------------------�������¹յ�
    // 1������ʹ���Ż����ı��߽����жϡ������жϷ���Ϊ������������������ͻ�����ж�Ϊ�յ�
    // ��������ɨ��
    int diff = 0;  // ��¼��ֵ

    for (int i = image_h-2; i > 1; i--)
    {
        diff = border_location_left_optimized[i-1][0] - border_location_left_optimized[i][0];
        if (diff > 5 || diff < -5)
        {
            diffPoint_temp[0] = i;
            diffPoint_temp[1] = border_location_left_optimized[i][0];
            break;
        }
    }
    // printf("left diff point 1: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    // 2�����������������ҵ��Ĺյ��Ƿ���ȷ
    // ���յ㸽���ĵ����������, ���ǰ���������
    // border_location_right_optimized[diffPoint_temp[0]][1]��ָ�յ���ԭ�߽��е�λ�ã����ҵ��߽�
    if (border_location_left_optimized[diffPoint_temp[0]][1] > 1 && border_location_left_optimized[diffPoint_temp[0]][1] < border_length-2)  // ��ֹ����Խ��
    {
        // printf("[border_location_left_optimized[diffPoint_temp[0]][1]+2]:%d\n", border_location_left_optimized[diffPoint_temp[0]][1]+2);
        if ((growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 4 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // �յ��ڶ�����
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 2 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 3 ||  // �յ���һ���
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 2 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 0 ||  // �յ�ǰ��һ���
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 7 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 0 ||  // �յ�ǰ�ڶ����
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 7 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 1 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 2) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] !=
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2])
            )
            {
                printf("+2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2]);
                printf("+1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1]);
                printf("-1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1]);
                printf("-2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2]);
                left_flag = 1;
            }  //�޶���
        else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //������������������
        }
    }
    else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //������������������
        }
    // printf("left diff point 2: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    // growth_direction_left[border_location_right_optimized[diffPoint_temp[0]][1]];

    // �������ֵ
    diffPoint_left[0] = diffPoint_temp[0];
    diffPoint_left[1] = diffPoint_temp[1];

    // ----------------------------------------------------------------------------------�������¹յ�
    // 1������ʹ���Ż����ı��߽����жϡ������жϷ���Ϊ������������������ͻ�����ж�Ϊ�յ�
    // ��������ɨ��
    diff = 0;  // ��¼��ֵ

    for (int i = image_h-2; i > 1; i--)
    {
        diff = border_location_right_optimized[i-1][0] - border_location_right_optimized[i][0];
        if (diff > 5 || diff < -5)
        {
            diffPoint_temp[0] = i;
            diffPoint_temp[1] = border_location_right_optimized[i][0];
            break;
        }
    }
    // printf("right diff point 1: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    // 2�����������������ҵ��Ĺյ��Ƿ���ȷ
    // ���յ㸽���ĵ����������, ���ǰ���������
    // border_location_right_optimized[diffPoint_temp[0]][1]��ָ�յ���ԭ�߽��е�λ�ã����ҵ��߽�
    // printf("[border_location_right_optimized[diffPoint_temp[0]][1]+2]:%d\n", border_location_right_optimized[diffPoint_temp[0]][1]+2);
    // printf("+2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2]);
    // printf("+1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1]);
    // printf("-1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1]);
    // printf("-2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2]);
    // printf("");
    if (border_location_right_optimized[diffPoint_temp[0]][1] > 1 && border_location_right_optimized[diffPoint_temp[0]][1] < image_h-2)  // ��ֹ����Խ��
    {
        if ((growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 4 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // �յ��ڶ�����
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 2 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 3 ||  // �յ���һ���
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 2 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 0 ||  // �յ�ǰ��һ���
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 7 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 0 ||  // �յ�ǰ�ڶ����
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 7 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 1 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 2) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] !=
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2])
            )
            {
                // printf("+2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2]);
                // printf("+1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1]);
                // printf("-1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1]);
                // printf("-2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2]);
                right_flag = 2;
            }  //�޶���
        else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //������������������
        }
    }
    else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //������������������
        }
    // printf("right diff point 2: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    // growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]];

    // �������ֵ
    diffPoint_right[0] = diffPoint_temp[0];
    diffPoint_right[1] = diffPoint_temp[1];

    // ����flag
    flag = left_flag + right_flag;
    // 3��������ҹյ�߶�������ȡ����·��ĵ���Ϊ�յ����
    if (flag ==  3)
    {
        if (diffPoint_right[0] - diffPoint_left[0] > 25)
        {
            left_flag = 0;
        }
        else if (diffPoint_left[0] - diffPoint_right[0] > 25)
        {
            right_flag = 0;
        }
    }
    // printf("diffPoint_right:%d\n", diffPoint_right[0]);
    // 4���˴�Ϊ��ֹ����Զ�������д����Ӱ������ɾ��
    // printf("diffPoint_left[0]: %d\n", diffPoint_left[0]);
    // printf("diffPoint_right[0]: %d\n", diffPoint_right[0]);
    if (diffPoint_left[0] < 25)
    {
        left_flag = 0;
    }
    if (diffPoint_right[0] < 25)
    {
        right_flag = 0;
    }

    flag = left_flag + right_flag;
    return flag;
}



// ---------------------------------------------------------------------------------------���߲���
// ����������б��
void two_slope(double x1, double y1, double x2, double y2, double* slope, double* intercept) {
    *slope = (y2 - y1) / (x2 - x1);
    *intercept = y1 - (*slope * x1);
}
// ����һ��Ҫע�⣬�˴��Ƕ�ά���飬���ǵ�0��ֵ��image_w����1��ֵ�ǵ�ǰ����ԭ�߽����������ǰֵ����������image_h
// �������㲹��
void FillLine(int line[image_h][2], int x1, int y1, int x2, int y2) {
    double slope, intercept;
    // printf("x1:%d y1:%d x2:%d y2:%d\n", x1, y1, x2, y2);
    // ������벻����Ҫ�����Զ�����
    if (y1 > y2) {
        int temp = x1;
        x1 = x2;
        x2 = temp;
        temp = y1;
        y1 = y2;
        y2 = temp;
    }
    // ���������
    if (x1 == x2)
    {
        for (int i = y1; i <= y2; i++) {
            line[i][0] = x1;
            // printf("line[%d][0]:%d\n", i, line[i][0]);
        }
        return;
    }
    two_slope(x1,y1,x2,y2,&slope,&intercept);
    for (int i = y1; i <= y2; i++) {
        line[i][0] = (i - intercept) / slope;
        // printf("line[%d][0]:%d\n", i, line[i][0]);
    }
}
// ---------------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------
// ֱ���ж�
// ��С���˷�����ֱ����ϵ�б�ʺͽؾ�
void LeastSquaresSlope(int32_t line[image_h], uint8_t start, uint8_t end, double* slope, double* intercept) {
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    uint8_t size = end - start + 1;
    for (int i = start; i <= end; i++) {
        sumX += i;  // �ۼӵ�Ӧ���� x ֵ���� i
        sumY += line[i];  // �ۼӵ��� y ֵ���� line[i]
        sumXY += i * line[i];  // ������ x * y
        sumX2 += i * i;  // ������ x ��ƽ��
    }
    double xMean = sumX / size;
    double yMean = sumY / size;
    *slope = (sumXY - size * xMean * yMean) / (sumX2 - size * xMean * xMean);
    *intercept = (yMean - *slope * xMean);
}
// ������ϵķ���
double variance(int32_t line[image_h], uint8_t start, uint8_t end, uint8_t step, double* slope, double* intercept) {
    LeastSquaresSlope(line, start, end, slope, intercept);
    double sum = 0;
    int size = end - start + 1;
    for (int i = start; i <= end; i += step) {
        double yFit = *slope * i + *intercept;  // �������ֱ���ϵ� y ֵ
        sum += pow(line[i] - yFit, 2);  // ����в��ƽ����
    }
    return sum * step / size;  // ���ر�׼��
}
// void LeastSquaresSlope(uint8_t line[image_h], uint8_t start, uint8_t end, double* slope, double* intercept) {
//     double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
//     uint8_t size = end - start + 1;
//     for (int i = start; i <= end; i++) {
//         sumX += i;  // �ۼӵ�Ӧ���� x ֵ���� i
//         sumY += line[i];  // �ۼӵ��� y ֵ���� line[i]
//         sumXY += i * line[i];  // ������ x * y
//         sumX2 += i * i;  // ������ x ��ƽ��
//     }
//     double xMean = sumX / size;
//     double yMean = sumY / size;
//     *slope = (sumXY - size * xMean * yMean) / (sumX2 - size * xMean * xMean);
//     *intercept = (yMean - *slope * xMean);
// }
// double variance(uint8_t line[image_h], uint8_t start, uint8_t end, uint8_t step, double* slope, double* intercept) {
//     LeastSquaresSlope(line, start, end, slope, intercept);
//     double sum = 0;
//     int size = end - start + 1;
//     for (int i = start; i <= end; i += step) {
//         double yFit = *slope * i + *intercept;  // �������ֱ���ϵ� y ֵ
//         sum += pow(line[i] - yFit, 2);  // ����в��ƽ����
//     }
//     return sum * step / size;  // ���ر�׼��
// }
// ---------------------------------------------------------------------------------------------

// ��¼��ǰ״̬�ı�־
// bool isInRightCircle = false;
// bool isInCrossing = false;
int inRightCircleFlag = 0;
int inCrossingFlag = 0;
int inLeftCircleFlag = 0;

int status3_bwline[2] = {0};  // black white line
long frame_number = 0;
int amount_fps[2] = {0, 0};  // ���ڼ�¼״̬0ʱ�����ҵ�5���ҹյ�

int right_point_row1[2] = {0};  // �����¼���ǳ�̬�µ����
// �ҳ��һ����������
int LeftPoint_RightCircle[2];

// �ҳ��󻷵������ҵ�
int RightPoint_LeftCircle[2];

void FindRightCircleLeftPoint(int flag)
{
    memset(LeftPoint_RightCircle, 0, sizeof(LeftPoint_RightCircle));
    if (flag == 0)
    {
        diffPoint_right[0] = right_point_row1[0];
        diffPoint_right[1] = right_point_row1[1];
    }
    else if(flag == 2)
    {
        diffPoint_right[0] = start_point_l[0] - 1;
        diffPoint_right[1] = start_point_l[1] - 2;
    }
    // �ӵ�ǰ�ҹյ�������Բ
    // ��ǰ�ҹյ㣺diffPoint_right
    for (int i = diffPoint_right[0]-1; i > 1; i--)
    {
        if (inRightCircleFlag == 1)
        {
            if (i < (image_h/4)*3 && bin_image[i][diffPoint_right[1]] == 0)  // �������ǰ��һ���ж�������״̬1ʱ�����һ����ͼ���ϵĲ���
            {
                LeftPoint_RightCircle[0] = i;
                LeftPoint_RightCircle[1] = diffPoint_right[1];
                // printf("Found Circle\n");
                break;  // �ҵ��˼�¼�㲢����ѭ��
            }
        }
        else
        {
            if (bin_image[i][diffPoint_right[1]] == 0)
            {
                LeftPoint_RightCircle[0] = i;
                LeftPoint_RightCircle[1] = diffPoint_right[1];
                // printf("Found Circle\n");
                break;  // �ҵ��˼�¼�㲢����ѭ��
            }
        }

    }
    // ���ҵ��İ�����������ߵĵ�
    /*
        [2]  [3]  [4]

        [1]  [ ]  [5]

      [0/8]  [7]  [6]
    */
    int neighbor[9][2] = {{1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}};
    for (int j = diffPoint_right[0]-1; j > 1; j--)
    {
        for (int i = 0; i < 8; i++)
        {
            if (bin_image[LeftPoint_RightCircle[0] + neighbor[i][0]][LeftPoint_RightCircle[1] + neighbor[i][1]] == 255 && bin_image[LeftPoint_RightCircle[0] + neighbor[i+1][0]][LeftPoint_RightCircle[1] + neighbor[i+1][1]] == 0)
            {
                if (LeftPoint_RightCircle[1] >= LeftPoint_RightCircle[1] + neighbor[i][1])
                {
                    LeftPoint_RightCircle[0] = LeftPoint_RightCircle[0] + neighbor[i][0];
                    LeftPoint_RightCircle[1] = LeftPoint_RightCircle[1] + neighbor[i][1];
                    // printf("LeftPoint_RightCircle:%d %d\n", LeftPoint_RightCircle[0], LeftPoint_RightCircle[1]);
                    break;  // ���ԭ���ĵ��ڰ��������ұߣ���Ѹð�������¼����
                }
                else
                {
                    return;
                }
            }
        }
    }
}



// uint8_t start_point = image_h - 2;  // ��¼��ʼ��
int diffPoint_rightUp[2] = {0, 90};
void findS3TurningPoint_right(int start_point)
{
    int found_flag = 0;  // ָʾ�Ƿ��ҵ��յ㣬����Ƕ��ѭ��
    for (int i = start_point; i > 1; i--)
    {
        for (int j = border_location_left_optimized[i][0]+1; j < image_w - 15; j++)
        {
            // printf("diffPoint_rightUp: %d, %d\n", i, j);
            if (bin_image[i][j] == 0)
            {
                if (border_location_right_optimized[i][0] == j)
                {
                    printf("break!\n");
                    break;
                }
                printf("diffPoint_rightUp: %d, %d\n", i, j);
                // printf("border_location_left_optimized[%d][0]: %d\n",i , border_location_left_optimized[i][0]);
                found_flag = 1;
                diffPoint_rightUp[0] = i;
                diffPoint_rightUp[1] = j;
                break;
            }
        }
        if (found_flag)
        {
            if(amount_fps[0] == 0)
            {
                amount_fps[0]++;
                continue;
            }
            else
            {
                amount_fps[0] = 0;
                break;
            }
        }
    }
}
// �ж���Բ��
void rightCircleHandle(int TurningPointFlag)
{
    if (inRightCircleFlag == 0 && TurningPointFlag == 2 && inCrossingFlag == 0 && inLeftCircleFlag == 0)
    {
        // �Ϸ����������ҵ�5�ιյ�
        // if (amount_fps[0] == 0)
        // {

        //     amount_fps[0]++;
        //     amount_fps[1] = frame_number;
        // }
        // else
        // {
        //     if (frame_number - amount_fps[1] > 1)
        //     {
        //         amount_fps[0] = 0;
        //     }
        //     else if (amount_fps[0] >= 5)
        //     {
        //         amount_fps[0] = 0;
        //         inRightCircleFlag = 1;  // ����״̬1
        //         getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        //         right_point_row1[0] = start_point_r[0];
        //         right_point_row1[1] = start_point_r[1];
        //     }
        //     else
        //     {
        //         amount_fps[0]++;
        //         amount_fps[1] = frame_number;
        //     }
        // }

        // �·������ҵ��Ҳ�յ�󣬼������Ƿ���
        int diff = 0;
        // printf("ah\n");
        // ����Ƿ���ִ���ȶ���
        for (int i = image_h / 5; i < image_h-2; i++)
        {
            diff = border_location_left_optimized[i][0] - border_location_left_optimized[i+1][0];
            printf("diff=%d\n", diff);
            if (diff > 5 || diff < -5)
            {
                printf("What a pity!\n");
                return;
            }
        }
        // ������յ��Ӧλ���Ϸ��Ƿ�����
        // bin_image[diffPoint_right[0]-2][border_location_left_optimized[diffPoint_right[0]][0]]
        if (bin_image[diffPoint_right[0]-2][border_location_left_optimized[diffPoint_right[0]][0]-1] == 255)
        {
            printf("%d, %d\n", diffPoint_right[0]-2, border_location_left_optimized[diffPoint_right[0]][0]);
            printf("What a pity!\n");
            return;
        }
        // �������Ƿ�Ϊһ��ֱ��
        // double var = 0;
        // double slope, intercept;
        // var = variance(border_location_left_optimized, 15, image_h - 4, 1, &slope, &intercept);
        inRightCircleFlag = 1;
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        right_point_row1[0] = start_point_r[0];
        right_point_row1[1] = start_point_r[1];

        FindRightCircleLeftPoint(1);
    }
    else if (inRightCircleFlag == 1)
    {
        if(TurningPointFlag == 2 || TurningPointFlag == 3)
        {
            // printf("in if \n");
            FindRightCircleLeftPoint(1);
            // printf("LeftPoint_RightCircle[1]:%d\n", LeftPoint_RightCircle[1]);
            // printf("LeftPoint_RightCircle[0]:%d\n", LeftPoint_RightCircle[0]);
            FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
            // FindRightCircleLeftPoint(3);
            // FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
            // FindRightCircleLeftPoint(1);
        }
        else
        {
            // printf("in else \n");
            FindRightCircleLeftPoint(0);
            FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
            // FindRightCircleLeftPoint(3);
            // FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
            if (amount_fps[0] == 0)
            {
                amount_fps[0]++;
                amount_fps[1] = frame_number;
            }
            else
            {
                if (frame_number - amount_fps[1] > 1)
                {
                    amount_fps[0] = 0;
                }
                else if (amount_fps[0] >= 2 )
                {
                    amount_fps[0] = 0;
                    getStartPoint();  //��ʺɽ����ô�ѳɵġ�
                    if (start_point_r[1] > 90)
                    {
                        inRightCircleFlag = 2;  // ����״̬2
                    }

                }
                else
                {
                    amount_fps[0]++;
                    amount_fps[1] = frame_number;
                }
            }

            // amount_fps[0] = 0;
            // inRightCircleFlag = 2;  // ����״̬2

        }
    }
    else if (inRightCircleFlag == 2)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        // printf("")
        if (start_point_r[1] < 90 || border_location_right_optimized[image_h /2][0] < 90)
        {
            inRightCircleFlag = 3;
            rightCircleHandle(TurningPointFlag);
            return;
        }
        FindRightCircleLeftPoint(0);
        FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
        // �ұ߽粿��
        /*
            [0/8]  [1]  [2]

            [7]  [ ]  [3]

            [6]  [5]  [4]
        */
        int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
        memset(border_location_right, 0, sizeof(border_location_right));
        border_count_right = 0;
        start_point_r[0] = LeftPoint_RightCircle[0];
        start_point_r[1] = LeftPoint_RightCircle[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
                {
                    border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                    border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                    // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                    start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                    start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                    border_count_right++;
                    break;
                }
            }
        }
        optimizeBorder();




        // ####################################################################
        int start_point = image_h - 2;  // ��¼��ʼ��
        start_point = LeftPoint_RightCircle[0];
        printf("start_point = LeftPoint_RightCircle[0]  %d\n", start_point);
        for (int i = start_point; i > 1; i--)
        {
            if (border_location_right_optimized[i][0] > 90)
            {
                start_point = i-2;
                printf("start_point  test1: %d\n", start_point);
                break;
            }
        }
        // ####################################################################
        // FindRightCircleLeftPoint(3);
        // FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);

        // ####################################################################
        // ����������10����ǰ��ֱ����Ϊ��ʼ�㣨��ʱ�Ҳ��Ѿ����ߣ�
        // getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        // printf("start_point_l[1]: %d\n", start_point_l[1]);
        // uint8_t start_point = image_h - 2;  // ��¼��ʼ��
        // diffPoint_rightUp[2] = {0, 90};
        diffPoint_rightUp[0] = 0;
        diffPoint_rightUp[1] = 90;
        // && start_point_r[1] > 90
        // printf("start_point_l[1]:%d\n", start_point_l[1]);
        // if (start_point_l[1] < 10)
        // // if (false)
        // {
        //     printf("In if\n");
        //     // ʹ���Ż�������߽�����Ѱ�ҹյ�
        //     findS3TurningPoint_right(start_point);
        //     findS3TurningPoint_right(diffPoint_rightUp[0]-2);
        //     // �ҵ��յ��ʼ����
        //     // printf("%d\n", start_point);
        //     // FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 9, start_point_l[0]);
        //     FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        // }
        // else
        // {
        // ����ʹ���ұ߽�Ѱ����ʼ��
        // printf("In else\n");

        // ʹ���Ż�������߽�����Ѱ�ҹյ�
        printf("start_point:%d\n", start_point);
        findS3TurningPoint_right(start_point);
        // printf("")
        // findS3TurningPoint_right(diffPoint_rightUp[0]-2);
        // �ҵ��յ��ʼ����
        // printf("%d\n", start_point);
        // FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], start_point_l[1], start_point_l[0]);
        // FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        if (diffPoint_rightUp[0]+18 < image_h-1)
        {
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], border_location_left_optimized[diffPoint_rightUp[0]+18][0], diffPoint_rightUp[0]+18);
        }
        else{
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        }
        // }
        // �ҵ��յ��������������
        // ������߽磬ֱ���������Ҳ�ĵ�
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        memset(border_location_left, 0, sizeof(border_location_left));
        border_count_left = 0;
        start_point_l[0] = diffPoint_rightUp[0];
        start_point_l[1] = diffPoint_rightUp[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                {
                    start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                    start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                    border_location_left[border_count_left][0] = start_point_l[0];
                    border_location_left[border_count_left][1] = start_point_l[1];
                    border_count_left++;
                    break;
                }
            }
            if (start_point_l[1] > 90)
            {
                printf("break!!!\n");
                break;
            }
        }
        optimizeBorder();
        // ####################################################################

    }
    else if (inRightCircleFlag == 3)
    {
        // FillLine(border_location_left_optimized, image_w / 4, image_h - 1, image_w - 1, 0);
        // ���ж���ǰ�������Ҫ�л��������һ֡����ֹ���ֲ��ÿ�������
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (status3_bwline[0] == 0)
        {
            if (start_point_l[1] > 5)
            {
                status3_bwline[0] = 1;
                status3_bwline[1] = frame_number;
            }
        }
        else if(status3_bwline[0] == 1)
        {
            if (start_point_l[1] <= 5)
            {
                status3_bwline[0] = 2;
            }
        }
        else if(status3_bwline[0] == 2)
        {
            if (start_point_l[1] > 5)
            {
                status3_bwline[0] = 0;
                printf("to status4\n");
                inRightCircleFlag = 4;  // ����״̬4
                return;  // ֱ���˳�
            }
        }

        int start_point = image_h - 2;  // ��¼��ʼ��
        start_point = LeftPoint_RightCircle[0];
        for (int i = start_point; i > 1; i--)
        {
            if (border_location_right_optimized[i][0] >= 85)
            {
                start_point = i-2;
                // printf("start_point: %d\n", start_point);
                break;
            }
        }

        diffPoint_rightUp[0] = 0;
        diffPoint_rightUp[1] = 90;

        findS3TurningPoint_right(start_point);
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (diffPoint_rightUp[0]+35 < image_h-1 && start_point_l[1] >= 5)
        {
            printf("---if 1\n");
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], border_location_left_optimized[diffPoint_rightUp[0]+35][0], diffPoint_rightUp[0]+35);
        }
        else if(start_point_l[1] < 5 && diffPoint_rightUp[1] > 3)
        {
            printf("---else if 2\n");
            printf("diffPoint_rightUp:%d\n", diffPoint_rightUp[0]);
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w/5, image_h-1);
        }
        else if (diffPoint_rightUp[1] > 3)
        {
            printf("---else if 3\n");
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        }
        else{
            printf("---else 4\n");
            for (int i = image_h - 2; i > 5; i--)
            {
                if(bin_image[i][image_w/5] == 0)
                {
                    diffPoint_rightUp[0] = i;
                    diffPoint_rightUp[1] = image_w/5;
                    printf("diffPoint_rightUp0:%d\n", diffPoint_rightUp[0]);
                    printf("diffPoint_rightUp1:%d\n", diffPoint_rightUp[1]);
                    break;
                }
            }
        }
        // }
        // �ҵ��յ��������������
        // ������߽磬ֱ���������Ҳ�ĵ�
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        memset(border_location_left, 0, sizeof(border_location_left));
        border_count_left = 0;
        start_point_l[0] = diffPoint_rightUp[0];
        start_point_l[1] = diffPoint_rightUp[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                {
                    start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                    start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                    border_location_left[border_count_left][0] = start_point_l[0];
                    border_location_left[border_count_left][1] = start_point_l[1];
                    border_count_left++;
                    break;
                }
            }
            if (start_point_l[1] > 90)
            {
                printf("break!!!\n");
                break;
            }
        }
        optimizeBorder();
        // �ֳ����������һ�������½����ߣ�һ�������½�û�ߵ�����һ�����绷��һ�������½�û������һ��Ҳû�绷
        // #########################################################################################################################################################
        // uint8_t start_point = image_h - 2;  // ��¼��ʼ��
        // diffPoint_rightUp[0] = 0;
        // diffPoint_rightUp[1] = 90;
        // // && start_point_r[1] > 90
        // printf("start_point_l[1]:%d\n", start_point_l[1]);
        // if (start_point_l[1] < 10)
        // {
        //     printf("In if\n");
        //     // ʹ���Ż�������߽�����Ѱ�ҹյ�
        //     findS3TurningPoint_right(start_point);
        //     // �ҵ��յ��ʼ����
        //     // printf("%d\n", start_point);
        //     // FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 9, start_point_l[0]);
        //     FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        // }
        // else
        // {
        //     // ����ʹ���ұ߽�Ѱ����ʼ��
        //     printf("In else\n");
        //     for (uint8_t i = image_h - 2; i > 1; i--)
        //     {
        //         if (border_location_right_optimized[i][0] >= 90)
        //         {
        //             start_point = i;
        //             // printf("start_point: %d\n", start_point);
        //             break;
        //         }
        //     }
        //     // ʹ���Ż�������߽�����Ѱ�ҹյ�
        //     findS3TurningPoint_right(start_point);
        //     // �ҵ��յ��ʼ����
        //     // printf("%d\n", start_point);
        //     // FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], start_point_l[1], start_point_l[0]);
        //     FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 3, start_point_l[0]);
        // }
        // // �ҵ��յ��������������
        // // ������߽磬ֱ���������Ҳ�ĵ�
        // /*
        //     [2]  [1]  [0/8]

        //     [3]  [ ]  [7]

        //     [4]  [5]  [6]
        // */
        // int32_t neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        // memset(border_location_left, 0, sizeof(border_location_left));
        // border_count_left = 0;
        // start_point_l[0] = diffPoint_rightUp[0];
        // start_point_l[1] = diffPoint_rightUp[1];
        // for (uint8_t j = 0; j < border_length-1; j++)
        // {
        //     for (uint8_t i = 0; i < 8; i++)
        //     {
        //         if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
        //         {
        //             start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
        //             start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
        //             border_location_left[border_count_left][0] = start_point_l[0];
        //             border_location_left[border_count_left][1] = start_point_l[1];
        //             border_count_left++;
        //             break;
        //         }
        //     }
        //     if (start_point_l[1] > 90)
        //     {
        //         printf("break!!!\n");
        //         // printf("start_point_l[1]: %d\n", start_point_l[1]);
        //         // printf("start_point_l[0]: %d\n", start_point_l[0]);
        //         break;
        //     }
        // }
        // // ���border_location_left
        // // for (int i = 0; i < border_length; i++)
        // // {
        // //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
        // // }
        // optimizeBorder();
        // #########################################################################################################################################################

    }
    else if (inRightCircleFlag == 4)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
            //  || border_location_left_optimized[image_h-15][0] < 5
            // || border_location_right_optimized[image_h-15][0] > 89)
        if ((border_location_left_optimized[image_h-10][0] < 5) &&
            (border_location_right_optimized[image_h-10][0] > 89))
        {
            FillLine(border_location_left_optimized, image_w - 1, 0, image_w / 5, image_h - 1);
            inRightCircleFlag = 5;
            return;
        }

        int right_point = 0;

        // ����4  ѧ����code
        for (int i = 5; i < border_length; i++)
        {
            if (border_location_left[i][1] >= border_location_left[i - 1][1] && border_location_left[i][1] > border_location_left[i + 1][1]
                && border_location_left[i][1] > border_location_left[i - 10][1]) {   //������� һ��Ҫ�е��ں�
                int row = border_location_left[i][0];
                int col = border_location_left[i][1];
                // printf("status4:  row: %d, col: %d\n", row, col);
                //��һ����Ϊ���½ǹ̶��� �ڶ�����Ϊ���Ϲ̶�
                FillLine(border_location_left_optimized, image_w - 5, 5, col, row);  //ѹ��
                break;
            }
        }
    }
    else if (inRightCircleFlag == 5)
    {
        // ��������
        double slope_stick = 0;
        double slope_find = 0;
        double intercept_stick = 0;
        double intercept_find = 0;
        two_slope(image_w - 1, 0, 0, image_h - 1, &slope_stick, &intercept_stick);
        // printf("border_location_left_optimized1:%d,%d\n", border_location_left_optimized[(image_h-5)*4][0], (image_h-5)*4);
        // printf("border_location_left_optimized2:%d,%d\n", border_location_left_optimized[image_h/2][0], image_h/2);
        // two_slope(border_location_left_optimized[(image_h-5)*4][0], border_location_left_optimized[(image_h-5)*4][0], border_location_left_optimized[image_h/2][1], border_location_left_optimized[image_h/2][0], &slope_find, &intercept_find);
        two_slope(border_location_left_optimized[(image_h/3) *2][0], (image_h/3)*2, border_location_left_optimized[image_h/3][0], image_h/3, &slope_find, &intercept_find);
        // printf("slope_find:%f\n", slope_find);
        // printf("slope_stick:%f\n", slope_stick);
        if ((!isinf(slope_find)) && (slope_find < -0.8))
        {
            inRightCircleFlag = 6;
            // printf("slope_find - slope_stick: %f\n", slope_find - slope_stick);
            return;
        }
        // for (uint8_t i = image_h-20; i > 20; i++)
        // {
        //     diff = border_location_left_optimized[i-1][0] - border_location_left_optimized[i][0];
        //     if (diff > 10 || diff < -10)
        //     {
        //         FillLine(border_location_left_optimized, i, border_location_left_optimized[i][0], 0, image_w/2);
        //         break;
        //     }
        //     else
        //     {
        //         inRightCircleFlag = 6;
        //     }
        // }
        // printf("TurningPointFlag: %d\n", TurningPointFlag);
        // if (TurningPointFlag == 1 || TurningPointFlag == 3)
        // {
        //     FillLine(border_location_left_optimized, diffPoint_left[1], diffPoint_left[0], image_w/2, 3);
        //     // break;
        // }
        // else
        // {
        //     inRightCircleFlag = 6;
        // }

        // д����б�Ż���
        // ע�ⲻ���Ǵ������½ǿ�ʼд���������м俿����Ȼ����ȥ
        // FillLine(border_location_left_optimized, image_w - 1, 0, image_w / 5, image_h - 1);
        FillLine(border_location_left_optimized, image_w - 1, 0, 0, image_h - 1);


    }
    else if (inRightCircleFlag == 6)
    {
        if (border_location_right_optimized[image_h-5][0] > 90)
        {
            printf("Enter status 7\n");
            inRightCircleFlag = 7;
        }
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        int start_point = image_h - 2;  // ��¼��ʼ��
        for (int i = image_h - 2; i > 1; i--)
        {
            if (border_location_right_optimized[i][0] >= 90)
            {
                start_point = i;
                // printf("start_point: %d\n", start_point);
                break;
            }
        }
        diffPoint_rightUp[0] = 0;
        diffPoint_rightUp[1] = 90;
        findS3TurningPoint_right(start_point);
        FillLine(border_location_right_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], start_point_r[1], start_point_r[0]);
        // FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
        // FillLine(border_location_left_optimized, image_h-5, 0, 0, image_w/2);
        // FillLine(border_location_left_optimized, image_w / 3 * 2, 0, image_w / 9, image_h - 1);
        // printf("border_location_left_optimized[image_h-2][0]: %d\n", border_location_left_optimized[image_h-2][0]);
        // if (start_point_l[1] > 10)

    }
    else if (inRightCircleFlag == 7)
    {
        // ���½�����㣬�˳�����
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_r[1] < 90 || border_location_right_optimized[image_h/2][0] < 90)
        {
            inRightCircleFlag = 0;
            return;
            // printf("start_point_r[1]: %d\n", start_point_r[1]);
        }
        // FindRightCircleLeftPoint(0);
        int start_point = image_h - 2;  // ��¼��ʼ��
        diffPoint_rightUp[0] = 0;
        diffPoint_rightUp[1] = 90;
        findS3TurningPoint_right(start_point);
        FillLine(border_location_right_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], right_point_row1[1], right_point_row1[0]);
    }
}

// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------���Բ��֣���Ҫֱ����ֲ-----------------------------------------------------------------------------------
// ��¼��ǰ״̬�ı�־
// bool isInRightCircle = false;
// bool isInCrossing = false;


// uint8_t status3_bwline[2] = {0};  // black white line
// long frame_number = 0;
// uint32_t amount_fps[2] = {0, 0};  // ���ڼ�¼״̬0ʱ�����ҵ�5���ҹյ�

int left_point_row1[2] = {0};  // �����¼���ǳ�̬�µ����
// �ҳ��󻷵������ҵ�
// uint8_t RightPoint_LeftCircle[2];

void FindLeftCircleRightPoint(uint8_t flag)
{
    memset(RightPoint_LeftCircle, 0, sizeof(RightPoint_LeftCircle));
    if (flag == 0)
    {
        // diffPoint_left[0] = left_point_row1[0]-1;
        // diffPoint_left[1] = left_point_row1[1]-2;

        diffPoint_left[0] = left_point_row1[0]-1;
        diffPoint_left[1] = image_w/8;
    }
    else if(flag == 2)
    {
        diffPoint_left[0] = start_point_l[0] - 1;
        diffPoint_left[1] = start_point_l[1] - 2;
    }
    printf("flag:%d\n", flag);
    printf("diffPoint_left[0]:%d\n", diffPoint_left[0]);
    printf("diffPoint_left[1]:%d\n", diffPoint_left[1]);
    // �ӵ�ǰ��յ�������Բ
    // ��ǰ��յ㣺diffPoint_left
    for (int i = diffPoint_left[0]-1; i > 1; i--)
    {
        if (inLeftCircleFlag == 1)
        {
            // �ص��ע����ж�����-------------------------------------------
            if (i < (image_h/4)*3 && bin_image[i][diffPoint_left[1]] == 0)  // �������ǰ��һ���ж�������״̬1ʱ�����һ����ͼ���ϵĲ���
            {
                RightPoint_LeftCircle[0] = i;
                RightPoint_LeftCircle[1] = diffPoint_left[1];
                // printf("Found Circle\n");
                break;  // �ҵ��˼�¼�㲢����ѭ��
            }
        }
        else
        {
            if (bin_image[i][diffPoint_left[1]] == 0)
            {
                RightPoint_LeftCircle[0] = i;
                RightPoint_LeftCircle[1] = diffPoint_left[1];
                // printf("Found Circle\n");
                break;  // �ҵ��˼�¼�㲢����ѭ��
            }
        }

    }
    // ���ҵ��İ����������ұߵĵ�
    /*
        [4]  [3]  [2]

        [5]  [ ]  [1]

        [6]  [7]  [0/8]
    */
    //                           0        1        2        3        4        5        6        7        8
    // int32_t neighbor[9][2] = {{1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}};
    int neighbor[9][2] = {{1, 1}, {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}};
    for (int j = diffPoint_left[0]-1; j > 1; j--)
    {
        for (int i = 0; i < 8; i++)
        {
            if (bin_image[RightPoint_LeftCircle[0] + neighbor[i][0]][RightPoint_LeftCircle[1] + neighbor[i][1]] == 255 && bin_image[RightPoint_LeftCircle[0] + neighbor[i+1][0]][RightPoint_LeftCircle[1] + neighbor[i+1][1]] == 0)
            {
                if (RightPoint_LeftCircle[1] <= RightPoint_LeftCircle[1] + neighbor[i][1])
                {
                    RightPoint_LeftCircle[0] = RightPoint_LeftCircle[0] + neighbor[i][0];
                    RightPoint_LeftCircle[1] = RightPoint_LeftCircle[1] + neighbor[i][1];
                    // printf("RightPoint_LeftCircle:%d %d\n", RightPoint_LeftCircle[0], RightPoint_LeftCircle[1]);
                    break;  // ���ԭ���ĵ��ڰ���������ߣ���Ѹð�������¼����
                }
                else
                {
                    printf("RightPoint_LeftCircle:%d %d\n", RightPoint_LeftCircle[0], RightPoint_LeftCircle[1]);
                    return;
                }
            }
        }
    }
}



// uint8_t start_point = image_h - 2;  // ��¼��ʼ��
int diffPoint_leftUp[2] = {0, 90};
void findS3TurningPoint_left(int start_point)
{
    int found_flag = 0;  // ָʾ�Ƿ��ҵ��յ㣬����Ƕ��ѭ��
    for (int i = start_point; i > 1; i--)
    {
        printf("border_location_right_optimized[i][0]=%d\n", border_location_right_optimized[i][0]);
        for (int j = border_location_right_optimized[i][0]-1; j > 15; j--)
        {
            // printf("diffPoint_rightUp: %d, %d\n", i, j);
            if (bin_image[i][j] == 0)
            {
                if (border_location_left_optimized[i][0] == j)
                {
                    printf("break!\n");
                    break;
                }
                printf("diffPoint_leftUp: %d, %d\n", i, j);
                // printf("border_location_left_optimized[%d][0]: %d\n",i , border_location_left_optimized[i][0]);
                found_flag = 1;
                diffPoint_leftUp[0] = i;
                diffPoint_leftUp[1] = j;
                break;
            }
        }
        if (found_flag)
        {
            if(amount_fps[0] == 0)
            {
                amount_fps[0]++;
                continue;
            }
            else
            {
                amount_fps[0] = 0;
                break;
            }
        }
    }
}
// �ж���Բ��
void leftCircleHandle(int TurningPointFlag)
{
    if (inLeftCircleFlag == 0 && TurningPointFlag == 1 && inCrossingFlag == 0 && inRightCircleFlag == 0)
    {
        // �·������ҵ��Ҳ�յ�󣬼������Ƿ���
        int diff = 0;
        for (int i = image_h / 5; i < image_h-2; i++)
        {
            diff = border_location_right_optimized[i][0] - border_location_right_optimized[i+1][0];
            if (diff > 5 || diff < -5)
            {
                // printf("What a pity!\n");
                return;
            }
        }
        // ����Ҳ�յ��Ӧλ���Ϸ��Ƿ�����
        if (bin_image[diffPoint_left[0]-2][border_location_right_optimized[diffPoint_left[0]][0]] == 255)
        {
            return;
        }
        inLeftCircleFlag = 1;
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        left_point_row1[0] = start_point_l[0];
        left_point_row1[1] = start_point_l[1];

        FindLeftCircleRightPoint(1);
    }
    else if (inLeftCircleFlag == 1)
    {
        if(TurningPointFlag == 1 || TurningPointFlag == 3)
        {
            FindLeftCircleRightPoint(1);
            // FillLine(border_location_left_optimized, diffPoint_left[1], diffPoint_left[0], RightPoint_LeftCircle[1], RightPoint_LeftCircle[0]);
            FillLine(border_location_left_optimized, RightPoint_LeftCircle[1], RightPoint_LeftCircle[0], diffPoint_left[1], diffPoint_left[0]);
        }
        else
        {
            FindLeftCircleRightPoint(0);
            // FillLine(border_location_left_optimized, diffPoint_left[1], diffPoint_left[0], RightPoint_LeftCircle[1], RightPoint_LeftCircle[0]);
            FillLine(border_location_left_optimized, RightPoint_LeftCircle[1], RightPoint_LeftCircle[0], diffPoint_left[1], diffPoint_left[0]);
            printf("RightPoint_LeftCircle[1]:%d\n",RightPoint_LeftCircle[1]);
            printf("RightPoint_LeftCircle[0]:%d\n",RightPoint_LeftCircle[0]);
            printf("diffPoint_left[1]:%d\n",diffPoint_left[1]);
            printf("diffPoint_left[0]:%d\n",diffPoint_left[0]);
            // amount_fps[0] = 0;
            // inLeftCircleFlag = 2;  // ����״̬2
            if (amount_fps[0] == 0)
            {
                amount_fps[0]++;
                amount_fps[1] = frame_number;
            }
            else
            {
                if (frame_number - amount_fps[1] > 1)
                {
                    amount_fps[0] = 0;
                }
                else if (amount_fps[0] >= 2)
                {
                    amount_fps[0] = 0;
                    getStartPoint();  //��ʺɽ����ô�ѳɵġ�
                    // printf("start_point_l[1] %d\n", start_point_l[1]);
                    // printf("border_location_left_optimized[image_h-5][0] %d\n", border_location_left_optimized[image_h-5][0]);
                    if (start_point_l[1] < 4)
                    {
                        // printf("border_location_left_optimized[image_h-5 < 3]:%d\n",border_location_left_optimized[image_h-5][0]);
                        inLeftCircleFlag = 2;  // ����״̬2
                    }
                }
                else
                {
                    amount_fps[0]++;
                    amount_fps[1] = frame_number;
                }
            }
        }
    }
    else if (inLeftCircleFlag == 2)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_l[1] > 3 || border_location_left_optimized[image_h /2][0] > 3)
        {
            inLeftCircleFlag = 3;
            leftCircleHandle(TurningPointFlag);
            return;
        }
        FindLeftCircleRightPoint(0);
        FillLine(border_location_left_optimized, RightPoint_LeftCircle[1], RightPoint_LeftCircle[0], diffPoint_left[1], diffPoint_left[0]);

        // ������߽磬ֱ���������Ҳ�ĵ�
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        memset(border_location_left, 0, sizeof(border_location_left));  // ��ʼ������
        border_count_left = 0;
        start_point_l[0] = RightPoint_LeftCircle[0];
        start_point_l[1] = RightPoint_LeftCircle[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                {
                    start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                    start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                    border_location_left[border_count_left][0] = start_point_l[0];
                    border_location_left[border_count_left][1] = start_point_l[1];
                    border_count_left++;
                    break;
                }
            }
            if (start_point_l[1] > 90)
            {
                printf("break!!!\n");
                break;
            }
        }
        optimizeBorder();

        int start_point = image_h - 2;  // ��¼��ʼ��
        start_point = RightPoint_LeftCircle[0];
        for (int i = start_point; i > 1; i--)
        {
            if (border_location_left_optimized[i][0] < 4)
            {
                start_point = i-2;
                // printf("start_point: %d\n", start_point);
                break;
            }
        }

        diffPoint_leftUp[0] = 0;
        diffPoint_leftUp[1] = 0;

        findS3TurningPoint_left(start_point);

        printf("diffPoint_leftUp[0]:%d\n", diffPoint_leftUp[0]);
        printf("diffPoint_leftUp[1]:%d\n", diffPoint_leftUp[1]);

        if (diffPoint_leftUp[0]+18 < image_h-1)
        {
            FillLine(border_location_right_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], border_location_right_optimized[diffPoint_leftUp[0]+18][0], diffPoint_leftUp[0]+18);
        }
        else{
            FillLine(border_location_right_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], image_w / 3 * 2, start_point_r[0]);
        }

        // �����ұ߽磬ֱ�����������ĵ�

        // int32_t neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
        memset(border_location_right, 0, sizeof(border_location_right));
        border_count_right = 0;
        start_point_r[0] = diffPoint_leftUp[0];
        start_point_r[1] = diffPoint_leftUp[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255 && bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0)
                {
                    start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                    start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                    border_location_right[border_count_right][0] = start_point_r[0];
                    border_location_right[border_count_right][1] = start_point_r[1];
                    border_count_right++;
                    break;
                }
            }
            if (start_point_r[1] < 4)
            {
                printf("break!!!\n");
                break;
            }
        }
        optimizeBorder();
    }
    else if (inLeftCircleFlag == 3)
    {
        // FillLine(border_location_right_optimized, image_w / 4*3, image_h - 1, 0, 0);
        // ���ж���ǰ�������Ҫ�л��������һ֡����ֹ���ֲ��ÿ�������
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (status3_bwline[0] == 0)
        {
            // if (start_point_r[1] < 89)
            if (border_location_right_optimized[image_h /2][0] < 89 || start_point_r[1] < 89)
            {
                status3_bwline[0] = 1;
                status3_bwline[1] = frame_number;
            }
        }
        else if(status3_bwline[0] == 1)
        {
            if (start_point_r[1] >= 89)
            {
                status3_bwline[0] = 2;
            }
        }
        else if(status3_bwline[0] == 2)
        {
            if (start_point_r[1] < 89)
            {
                status3_bwline[0] = 0;
                printf("to status4\n");
                inLeftCircleFlag = 4;  // ����״̬4
                return;  // ֱ���˳�
            }
        }
        int start_point = image_h - 2;  // ��¼��ʼ��
        start_point = RightPoint_LeftCircle[0];
        for (int i = start_point; i > 1; i--)
        {
            if (border_location_left_optimized[i][0] <= 9)
            {
                start_point = i-2;
                // printf("start_point: %d\n", start_point);
                break;
            }
        }

        diffPoint_leftUp[0] = 0;
        diffPoint_leftUp[1] = 0;

        findS3TurningPoint_left(start_point);

        printf("diffPoint_leftUp[0]:%d\n", diffPoint_leftUp[0]);
        printf("diffPoint_leftUp[1]:%d\n", diffPoint_leftUp[1]);
        if (diffPoint_leftUp[0]+35 < image_h-1 && start_point_r[1] <= 89)
        // if (diffPoint_leftUp[0]+35< image_h-1)
        {
            FillLine(border_location_right_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], border_location_right_optimized[diffPoint_leftUp[0]+35][0], diffPoint_leftUp[0]+35);
        }
        else if(start_point_r[1] >= 89 && diffPoint_leftUp[1] < 91)
        {
            printf("else if 2\n");
            FillLine(border_location_right_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], image_w/5*4, image_h-1);
        }
        else if (diffPoint_leftUp[1] < 91)
        {
            printf("---else if 3\n");
            FillLine(border_location_right_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], image_w / 3 * 2, start_point_r[0]);
        }
        else{
            printf("---else 4\n");
            for (int i = image_h - 2; i > 5; i--)
            {
                if(bin_image[i][image_w/5*4] == 0)
                {
                    diffPoint_leftUp[0] = i;
                    diffPoint_leftUp[1] = image_w/5*4;
                    printf("diffPoint_leftUp0:%d\n", diffPoint_leftUp[0]);
                    printf("diffPoint_leftUp1:%d\n", diffPoint_leftUp[1]);
                    break;
                }
            }
        }
        // �ұ߽粿��
        /*
            [0/8]  [1]  [2]

            [7]  [ ]  [3]

            [6]  [5]  [4]
        */
        int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
        memset(border_location_right, 0, sizeof(border_location_right));
        border_count_right = 0;
        start_point_r[0] = diffPoint_leftUp[0];
        start_point_r[1] = diffPoint_leftUp[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
                {
                    border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                    border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                    // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                    start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                    start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                    border_count_right++;
                    break;
                }
            }
            if (start_point_r[1] < 4)
            {
                printf("break!!!\n");
                break;
            }
        }
        optimizeBorder();
    }
    else if (inLeftCircleFlag == 4)
        {
            getStartPoint();  //��ʺɽ����ô�ѳɵġ�
            // printf("border_location_right_optimized:%d\n", border_location_right_optimized[image_h-10][0]);
            if ((border_location_left_optimized[image_h-15][0] < 5) &&
                (border_location_right_optimized[image_h-15][0] > 89))
                // (start_point_l[1]< 10 && start_point_r[1]>84))
            {
                FillLine(border_location_right_optimized, 0, 0, image_w / 5 * 4, image_h - 1);
                inLeftCircleFlag = 5;
                return;
            }

            int left_point = 0;

            // ����4  ѧ����code
            for (int i = 5; i < border_length-1; i++)
            {
                if (border_location_right[i][1] <= border_location_right[i - 1][1] && border_location_right[i][1] < border_location_right[i + 1][1]
                    && border_location_right[i][1] < border_location_right[i - 10][1]) {   //������� һ��Ҫ�е��ں�
                    int row = border_location_right[i][0];
                    int col = border_location_right[i][1];
                    // printf("status4:  row: %d, col: %d\n", row, col);
                    //��һ����Ϊ���½ǹ̶��� �ڶ�����Ϊ���Ϲ̶�
                    FillLine(border_location_right_optimized, 4, 20, col, row);  //ѹ��
                    break;
                }
            }
        }
    else if (inLeftCircleFlag == 5)
    {
        // ��������
        double slope_stick = 0;
        double slope_find = 0;
        double intercept_stick = 0;
        double intercept_find = 0;
        two_slope(0, 0, image_w - 1, image_h - 1, &slope_stick, &intercept_stick);
        two_slope(border_location_right_optimized[(image_h/3) *2][0], (image_h/3)*2, border_location_right_optimized[image_h/3][0], image_h/3, &slope_find, &intercept_find);
        printf("slope_find:%f\n", slope_find);
        printf("slope_stick:%f\n", slope_stick);
        // if (!isnan(slope_find) && slope_find > slope_stick)
        if ((!isnan(slope_find)) && (slope_find > 0.8))
        {
            inLeftCircleFlag = 6;
            // printf("slope_find - slope_stick: %f\n", slope_find - slope_stick);
            return;
        }
        // д����б�Ż���
        // ע�ⲻ���Ǵ������½ǿ�ʼд���������м俿����Ȼ����ȥ
        FillLine(border_location_right_optimized, 0, 0, image_w -1, image_h - 1);
    }
    else if (inLeftCircleFlag == 6)
    {
        if (border_location_left_optimized[image_h-5][0] < 3)
        {
            printf("Enter status 7\n");
            inLeftCircleFlag = 7;
        }
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        int start_point = image_h - 2;  // ��¼��ʼ��
        for (int i = image_h - 2; i > 1; i--)
        {
            if (border_location_left_optimized[i][0] <= 4)
            {
                start_point = i;
                // printf("start_point: %d\n", start_point);
                break;
            }
        }
        diffPoint_leftUp[0] = 0;
        diffPoint_leftUp[1] = 0;
        findS3TurningPoint_left(start_point);
        // printf("diffPoint_leftUp[0]=%d, diffPoint_leftUp[1]=%d\n", diffPoint_leftUp[0], diffPoint_leftUp[1]);
        FillLine(border_location_left_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], start_point_l[1], start_point_l[0]);

    }
    else if (inLeftCircleFlag == 7)
    {
        // ���½�����㣬�˳�����
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_l[1] > 4)
        {
            inLeftCircleFlag = 0;
            return;
            // printf("start_point_r[1]: %d\n", start_point_r[1]);
        }
        // FindRightCircleLeftPoint(0);
        int start_point = image_h - 2;  // ��¼��ʼ��
        diffPoint_leftUp[0] = 0;
        diffPoint_leftUp[1] = 0;
        findS3TurningPoint_left(start_point);
        FillLine(border_location_left_optimized, diffPoint_leftUp[1], diffPoint_leftUp[0], left_point_row1[1], left_point_row1[0]);
    }
}
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------/���Բ��֣���Ҫֱ����ֲ---------------------------------------------------------------------------------------

// ʶ��ʮ�ֺ�Բ������Ҫ��������ʮ�����ֻʶ��һ���յ�ʱ����һ���Ƕ��ߵ�
// ��Բ���뻷��������ֻʶ��һ���յ㣬����һ���ǲ����ߵ�
int up_right_point[2] = {0};
int up_left_point[2] = {0};
void findUpTurningPoint(int start_l, int x_l, int start_r, int x_r)
{
    // printf("start_l=%d\n", start_l);
    // printf("start_r=%d\n", start_r);
    // �����Ϲյ�
    int flag = 0;
    for (int i = start_l; i > 2; i--)
    {
        // �ɷ���
        // if ((border_location_left_optimized[i-1][0] - border_location_left_optimized[i][0] < 5) &&
        //     (border_location_left_optimized[i][0] - border_location_left_optimized[i+2][0] > 10))
        // {
        //     // printf("Found Left Up Turning Point\n");
        //     up_left_point[0] = i;
        //     up_left_point[1] = border_location_left_optimized[i][0];
        //     printf("up_left_point[0]=%d\n", up_left_point[0]);
        //     printf("up_left_point[1]=%d\n", up_left_point[1]);
        //     break;
        // }
        // �·���
        if (bin_image[i][x_l] == 0 && bin_image[i-1][x_l] == 0 && bin_image[i-2][x_l] == 0)
        {
            up_left_point[0] = i-2;
            up_left_point[1] = x_l;
            for (int j = x_l; j < image_w - 2; j++)
            {
                if (bin_image[i-2][j] == 0 && bin_image[i-2][j+1] == 255 && bin_image[i-2][j+2] == 255)
                {
                    up_left_point[1] = j+1;
                    flag = 1;
                    break;
                }
            }
            if (flag)
            {
                flag = 0;
                break;
            }
        }
    }
    flag = 0;
    // �����Ϲյ�
    for (int i = start_r; i > 1; i--)
    {
        // if ((border_location_right_optimized[i][0] - border_location_right_optimized[i-1][0] < 5) &&
        //     (border_location_right_optimized[i+2][0] - border_location_right_optimized[i][0] > 10))
        // {
        //     // printf("Found right Up Turning Point\n");
        //     up_right_point[0] = i;
        //     up_right_point[1] = border_location_right_optimized[i][0];
        //     printf("up_right_point[0]=%d\n", up_right_point[0]);
        //     printf("up_right_point[1]=%d\n", up_right_point[1]);
        //     break;
        // }

        if (bin_image[i][x_r] == 0 && bin_image[i-1][x_r] == 0 && bin_image[i-2][x_r] == 0)
        {
            up_right_point[0] = i-2;
            up_right_point[1] = x_r;
            for (int j = x_r; j > 1; j--)
            {
                // printf("j:%d\n",j);
                if (bin_image[i-2][j] == 0 && bin_image[i-2][j-1] == 255 && bin_image[i-2][j-2] == 255)
                {
                    up_right_point[1] = j-1;
                    flag = 1;
                    break;
                }
            }
            if (flag)
            {
                flag = 0;
                break;
            }
        }
    }
        // printf("up_right_point[0]=%d\n", up_right_point[0]);
        // printf("up_right_point[1]=%d\n", up_right_point[1]);
}

// ʶ��ʮ��
int restore_point[4] = {image_h - 1, image_w / 5, image_h - 1, image_w / 5 * 4};

void crossingHandle(int TurningPointFlag)
{
    // if (inRightCircleFlag == 0 && inLeftCircleFlag == 0 && (TurningPointFlag == 1 || TurningPointFlag == 3) && inCrossingFlag == 0)
    if (inRightCircleFlag == 0 && inLeftCircleFlag == 0 && TurningPointFlag != 0 && inCrossingFlag == 0)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        // �����һЩ�ж�
        if(TurningPointFlag == 3)  // �����յ㶼ʶ�𵽣�����״̬1
        {
            inCrossingFlag = 1;
        }
        // else if ((start_point_r[1] > 90 && TurningPointFlag == 1) || (start_point_l[1] < 4 && TurningPointFlag == 2))  // ֻʶ����յ�ʱ������Ƿ�����㣬�����ˣ�����״̬1
        // {
        //     inCrossingFlag = 1;
        // }
    }
    else if (inCrossingFlag == 1)
    {
        // ��ʹ������ʮ�֣�Ҳ�п���ֻʶ����һ���յ㣬����Ҫ�����ж�
        // ��������յ㶼�����ҵ��ģ������ҵ��Ϸ��������յ�
        if (TurningPointFlag == 3)
        {
            restore_point[0] = diffPoint_left[0];
            restore_point[1] = diffPoint_left[1];
            restore_point[2] = diffPoint_right[0];
            restore_point[3] = diffPoint_right[1];
            findUpTurningPoint(diffPoint_left[0], diffPoint_left[1], diffPoint_right[0], diffPoint_right[1]);
            // ����
            FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], diffPoint_left[1], diffPoint_left[0]);
            FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], diffPoint_right[1], diffPoint_right[0]);
            // ����һ�������
            // ������߽磬ֱ���������Ҳ�ĵ�
            /*
                [2]  [1]  [0/8]

                [3]  [ ]  [7]

                [4]  [5]  [6]
            */
            int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
            memset(border_location_left, 0, sizeof(border_location_left));  // ��ʼ������
            border_count_left = 0;
            start_point_l[0] = up_left_point[0];
            start_point_l[1] = up_left_point[1];
            for (int j = 0; j < border_length-1; j++)
            {
                for (int i = 0; i < 8; i++)
                {
                    if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                    {
                        start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                        start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                        border_location_left[border_count_left][0] = start_point_l[0];
                        border_location_left[border_count_left][1] = start_point_l[1];
                        border_count_left++;
                        break;
                    }
                }
            }
            // �ұ߽粿��
            /*
                [0/8]  [1]  [2]

                [7]  [ ]  [3]

                [6]  [5]  [4]
            */
            int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
            memset(border_location_right, 0, sizeof(border_location_right));
            border_count_right = 0;
            start_point_r[0] = up_right_point[0];
            start_point_r[1] = up_right_point[1];
            for (int j = 0; j < border_length-1; j++)
            {
                for (int i = 0; i < 8; i++)
                {
                    if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
                    {
                        border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                        border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                        // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                        start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                        start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                        border_count_right++;
                        break;
                    }
                }
            }
            optimizeBorder();
        }
        else
        {
            if (restore_point[1] != 0 && restore_point[3] != 0)
            {
                findUpTurningPoint(restore_point[0], restore_point[1], restore_point[2], restore_point[3]);
                FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], restore_point[1], restore_point[0]);
                FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], restore_point[3], restore_point[2]);
                // ����һ�������
                // ������߽磬ֱ���������Ҳ�ĵ�
                /*
                    [2]  [1]  [0/8]

                    [3]  [ ]  [7]

                    [4]  [5]  [6]
                */
                int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
                memset(border_location_left, 0, sizeof(border_location_left));
                border_count_left = 0;
                start_point_l[0] = up_left_point[0];
                start_point_l[1] = up_left_point[1];
                for (int j = 0; j < border_length-1; j++)
                {
                    for (int i = 0; i < 8; i++)
                    {
                        if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                        {
                            start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                            start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                            border_location_left[border_count_left][0] = start_point_l[0];
                            border_location_left[border_count_left][1] = start_point_l[1];
                            border_count_left++;
                            break;
                        }
                    }
                    // if (start_point_l[1] > 90)
                    // {
                    //     printf("break!!!\n");
                    //     // printf("start_point_l[1]: %d\n", start_point_l[1]);
                    //     // printf("start_point_l[0]: %d\n", start_point_l[0]);
                    //     break;
                    // }
                }
                // ���border_location_left
                // for (int i = 0; i < border_length; i++)
                // {
                //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
                // }
                //�ұ߽粿��
                /*
                    [0/8]  [1]  [2]

                    [7]  [ ]  [3]

                    [6]  [5]  [4]
                */
                int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
                memset(border_location_right, 0, sizeof(border_location_right));
                border_count_right = 0;
                start_point_r[0] = up_right_point[0];
                start_point_r[1] = up_right_point[1];
                printf("start_point_r[0]: %d\n", start_point_r[0]);
                printf("start_point_r[1]: %d\n", start_point_r[1]);
                for (int j = 0; j < border_length-1; j++)
                {
                    for (int i = 0; i < 8; i++)
                    {
                        if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
                        {
                            border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                            border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                            // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                            start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                            start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                            border_count_right++;
                            break;
                        }
                    }
                }
                optimizeBorder();
            }
        }

        // ���������㣬����״̬1
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_l[1] < 4 && start_point_r[1] > 90)
        {
            inCrossingFlag = 2;
            memset(up_left_point, 0, sizeof(up_left_point));
            memset(up_right_point, 0, sizeof(up_right_point));
        }
    }
    else if (inCrossingFlag == 2)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        // ������಻�����ˣ��˳�ʮ��·��
        if (start_point_l[1] > 4 && start_point_r[1] < 90)
        {
            inCrossingFlag = 0;
            return;
        }
        findUpTurningPoint(start_point_l[0]-1, image_w/7, start_point_r[0]-1, image_w/7*6);
        // ����  ע������ģ�25, 58����75, 58�������ʵ��������
        FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], 25, 58);
        FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], 75, 58);
        // ����һ�������
        // ������߽磬ֱ���������Ҳ�ĵ�
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        int neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        memset(border_location_left, 0, sizeof(border_location_left));
        border_count_left = 0;
        start_point_l[0] = up_left_point[0];
        start_point_l[1] = up_left_point[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                {
                    start_point_l[0] = start_point_l[0] + neighbor_left[i][0];
                    start_point_l[1] = start_point_l[1] + neighbor_left[i][1];
                    border_location_left[border_count_left][0] = start_point_l[0];
                    border_location_left[border_count_left][1] = start_point_l[1];
                    border_count_left++;
                    break;
                }
            }
        }
        // �ұ߽粿��
        /*
            [0/8]  [1]  [2]

            [7]  [ ]  [3]

            [6]  [5]  [4]
        */
        int neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
        memset(border_location_right, 0, sizeof(border_location_right));
        border_count_right = 0;
        start_point_r[0] = up_right_point[0];
        start_point_r[1] = up_right_point[1];
        for (int j = 0; j < border_length-1; j++)
        {
            for (int i = 0; i < 8; i++)
            {
                if (bin_image[start_point_r[0] + neighbor_right[i][0]][start_point_r[1] + neighbor_right[i][1]] == 255  && (bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 0 || bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] == 151))
                {
                    border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i][0];
                    border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i][1];
                    // bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                    start_point_r[0] = start_point_r[0] + neighbor_right[i][0];
                    start_point_r[1] = start_point_r[1] + neighbor_right[i][1];
                    border_count_right++;
                    break;
                }
            }
        }
        // ���border_location_right
        // for (int i = 0; i < border_length; i++)
        // {
        //     printf("border_location_right[%d]: %d, %d\n", i, border_location_right[i][0], border_location_right[i][1]);
        // }

        // printf("optimizeBorder\n");
        optimizeBorder();


    }
}

// �������ж�
int zebraJudge()
{
    int amount = 0;  // ��¼�ڰ���������
    int threshold = 10;  // ÿ���������ֵ
    int over_threshold = 0;  // ��¼������ֵ������
    // ��һ����Χ�ڵ�ͼ������ж�
    for (int h = image_h - 3; h > 40; h--)
    {
        amount = 0;
        for (int w = 0; w < image_w-2; w++)
        {
            if ((bin_image[h][w] == 255 && bin_image[h][w+1] == 0) ||
                (bin_image[h][w] == 0 && bin_image[h][w+1] == 255))
            {
                // printf("Ciallo~\n");
                amount++;
            }
        }
        if(amount > threshold)
        {
            over_threshold++;
        }
    }
    // printf("over_threshold:%d\n", over_threshold);
    if (over_threshold > 4)
    {
        printf("return 1\n");
        return 1;
    }
    else
    {
        return 0;
    }
}



// ֱ���ж�
// uint8_t linear_flag = 0; �ڻ�ȡƫ��ֵ����
void linearJudge() {
    double var = 0;
    double slope, intercept;
    var = variance(middle_line, 15, image_h - 4, 1, &slope, &intercept);
    // printf("var:%f\n", var);
    if (var < 1) 
    {
        linear_flag = 1;
    }
    else if (var < 10)
    {
        linear_flag = 2;
    }
    else 
    {
        linear_flag = 0;
    }
}

// int isPodaoFlag = 0;
// void podaoJudge()
// {
//     if (inLeftCircleFlag == 0 && inRightCircleFlag == 0 && inCrossingFlag == 0)
//     {
//         int diff = 0;
//         int count = 0;
//         for (int i = 20; i < 26; i++)
//         {
//             diff = border_location_right_optimized[i][0] - border_location_left_optimized[i][0];
//             // printf("%d: %d\n", i, diff);
//             if (diff > 20 && diff < 40)
//             {
//                 count += 1;
//             }
//         }
//         if (count >= 5)
//         {
//             isPodaoFlag = 1;
//         }
//         else
//         {
//             isPodaoFlag = 0;
//         }
//     }
//     else
//     {
//         isPodaoFlag = 0;
//     }
// }
// int zhang_flag = 0;
// int zhang2_flag = 0;
int zhang2_flag_r = 0;
int zhang2_flag_l = 0;
void zhangai2_right()
{
    if (inLeftCircleFlag == 0 && inRightCircleFlag == 0 && inCrossingFlag == 0 && zhang2_flag_l == 0 && zhang2_flag_r == 0)
    {
        // �����ж϶������
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_right_optimized[i-1][0] - border_location_right_optimized[i][0];
            if (diff < -7)
            {
                zhangai_point[0] = i-1;
                zhangai_point[1] = border_location_right_optimized[i-1][0];
                // printf("zhangai_point[0]:%d\n", zhangai_point[0]);
                // printf("zhangai_point[1]:%d\n", zhangai_point[1]);
                // printf("diff=%d\n", diff);
                if (bin_image[zhangai_point[0]-4][border_location_left_optimized[zhangai_point[0]][0]-2] == 255)
                {
                    return;
                }
                if (zhangai_point[0]+2 < image_h-1 && bin_image[zhangai_point[0]+2][border_location_left_optimized[zhangai_point[0]][0]-1] == 255)
                {
                    return;
                }
                zhang2_flag_r = 1;  // ����״̬1
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_right_optimized, zhangai_point[1], zhangai_point[0], border_location_right_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
    }
    else if (zhang2_flag_r == 1)
    {
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_right_optimized[i-1][0] - border_location_right_optimized[i][0];
            if (diff < -7)
            {
                zhangai_point[0] = i-1;
                zhangai_point[1] = border_location_right_optimized[i-1][0];
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_right_optimized, zhangai_point[1], zhangai_point[0], border_location_right_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
        
        if (inRightCircleFlag >= 1 || inLeftCircleFlag >= 1 || deviation_number >= 15 || deviation_number <= -15)
        {
            zhang2_flag_r = 2;  // ����״̬2
        }
    }
    else if (zhang2_flag_r == 2)
    {
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_right_optimized[i-1][0] - border_location_right_optimized[i][0];
            if (diff < -7)
            {
                zhangai_point[0] = i-1;
                zhangai_point[1] = border_location_right_optimized[i-1][0];
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_right_optimized, zhangai_point[1], zhangai_point[0], border_location_right_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
        
        // ׼����״̬
        int flag = findTurningPoint();
        if (flag == 0 || flag == 1)
        {
            zhang2_flag_r = 0;
            inRightCircleFlag = 0;
            inLeftCircleFlag = 0;
        }
    }
}

void zhangai2_left()
{
    if (inLeftCircleFlag == 0 && inRightCircleFlag == 0 && inCrossingFlag == 0 && zhang2_flag_l == 0 && zhang2_flag_r == 0)
    {
        // �����ж϶������
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_left_optimized[i+1][0] - border_location_left_optimized[i][0];
            if (diff > 7)
            {
                zhangai_point[0] = i+1;
                zhangai_point[1] = border_location_left_optimized[i+1][0];
                if (bin_image[zhangai_point[0]-4][border_location_right_optimized[zhangai_point[0]][0]+2] == 255)
                {
                    return;
                }
                if (zhangai_point[0]+2 < image_h-1 && bin_image[zhangai_point[0]+2][border_location_left_optimized[zhangai_point[0]][0]+1] == 255)
                {
                    return;
                }
                zhang2_flag_l = 1;  // ����״̬1
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_left_optimized, zhangai_point[1], zhangai_point[0], border_location_left_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
    }
    else if (zhang2_flag_l == 1)
    {
        // �����ж϶������
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_left_optimized[i+1][0] - border_location_left_optimized[i][0];
            if (diff > 7)
            {
                zhangai_point[0] = i+1;
                zhangai_point[1] = border_location_left_optimized[i+1][0];
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_left_optimized, zhangai_point[1], zhangai_point[0], border_location_left_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
        if (inRightCircleFlag >= 1 || inLeftCircleFlag >= 1 || deviation_number >= 15 || deviation_number <= -15)
        {
            zhang2_flag_l = 2;  // ����״̬2
        }
    }
    else if (zhang2_flag_l == 2)
    {
        // �����ж϶������
        int diff = 0;
        int zhangai_point[2] = {0};
        for (int i = image_h-2; i > 5; i--)
        {
            diff = border_location_left_optimized[i+1][0] - border_location_left_optimized[i][0];
            if (diff > 7)
            {
                zhangai_point[0] = i+1;
                zhangai_point[1] = border_location_left_optimized[i+1][0];
                break;
            }
        }
        if (zhangai_point[0] < image_h -17)
        {
            FillLine(border_location_left_optimized, zhangai_point[1], zhangai_point[0], border_location_left_optimized[zhangai_point[0]+15][0], zhangai_point[0]+15);
        }
        
        // ׼����״̬
        int flag = findTurningPoint();
        if (flag == 0 || flag == 2)
        {
            zhang2_flag_l = 0;
            inRightCircleFlag = 0;
            inLeftCircleFlag = 0;
        }
    }
}
void zhangai2()
{
    zhangai2_left();
    zhangai2_right();
}
// int zhang_flag = 0;
// int zhang2_flag = 0;
// void zhangaiJudge()
// {
//     if(zhang_flag == 1)
//     {
//         podaoJudge();
//         if (isPodaoFlag == 0)
//         {
//             zhangai2();
//             zhang2_flag = 1;
//         }else{
//             zhang2_flag = 0;
//         }
//     }
//     else{
//         isPodaoFlag = 0;
//         zhang2_flag = 0;
//     }
// }
// Ԫ���ж�������
int zebra_flag = 0;  // ��¼�Ƿ�ʶ�𵽰�����
void elementsHandle()
{
    zebra_flag = zebraJudge();
    int TurningPointFlag = 0;
    TurningPointFlag = findTurningPoint();
    // printf("deviation_number=%f\n", deviation_number);
    printf("TurningPointFlag:%d\n", TurningPointFlag);
    printf("inRightCircleFlag:%d\n", inRightCircleFlag);
    printf("inLeftCircleFlag:%d\n", inLeftCircleFlag);
    printf("inCrossingFlag:%d\n", inCrossingFlag);
    // printf("linear_flag:%d\n", linear_flag);
    // printf("zebra_flag:%d\n", zebra_flag);
    printf("zebra_flag:%d\n", zebra_flag);
    crossingHandle(TurningPointFlag);
    rightCircleHandle(TurningPointFlag);
    leftCircleHandle(TurningPointFlag);
    // printf("up_right_point[0]: %d\n", up_right_point[0]);
    // printf("up_right_point[1]: %d\n", up_right_point[1]);
    
}

void image_process(void)
{
    // ----------��̬��ͼ����----------
    turn_to_bin();  //ͼ���ֵ��
    // adaptive_threshold(image, IMAGE_W, IMAGE_H, BLOCK_SIZE, C);  // ����Ӧ��ֵ��
    image_compress(image,bin_image); //ͼ��ѹ��
    // image_filter();
    image_draw_rectan(bin_image);
    // ������
    if (getStartPoint())
	{
		// printf("Start point found.\n");
        // printf("Found Left Start Point: %d,%d\n", start_point_l[0], start_point_l[1]);
        // printf("Found Right Start Point: %d,%d\n", start_point_r[0], start_point_r[1]);
        neighborSearch();
        // memset(border_location_left_optimized, 0, sizeof(border_location_left_optimized));
        // memset(border_location_right_optimized, 0, sizeof(border_location_right_optimized));
        memset(border_location_left_optimized, 0, sizeof(border_location_left_optimized));
        for(int i = 0; i < image_h; i++) {
            border_location_right_optimized[i][0] = image_w-1;
        }
	    optimizeBorder();
        elementsHandle();
        getMiddleLine();
        if (inRightCircleFlag == 3)
        {
            highestPoint = 0;
        }
        else
        {
            getHighestPoint();
        }

        printf("Highest Point:%d\n", highestPoint);
        getDeviation();

        // ��ʱ��ʩ
        if (inRightCircleFlag == 4)
        {
            // if (deviation_number > 28)
            // {
            //     deviation_number = 28;
            // }
        }
        
        if (inCrossingFlag == 0 && inRightCircleFlag == 0 && inLeftCircleFlag == 0)
        {
            linearJudge();
        }
        else
        {
            linear_flag = 0;
        }
        // zhangaiJudge();

        zhangai2();
    }
    else
    {
        printf("Find start point failed.\n");
        deviation_number = deviation_number_p;
    }
    printf("deviation_number=%f\n", deviation_number);
    frame_number++;
    // printToCSV(bin_image, "output.csv");
    // usleep(20000);
}

int main(){
    if (PVFLAG == 0)
    {
        inRightCircleFlag = 0;
        // inCrossingFlag = 2;
        // inLeftCircleFlag = 6;
        Mat roadImage = imread(ROAD);
        if (roadImage.empty())
        {
            cerr << "Failed to load image!" << endl;
            return -1;
        }
        dataTrans(roadImage, image);

        image_process();

        //��opencv��ʾ
        Mat Binimage(image_h, image_w, CV_8UC1);
        for (int i = 0; i < image_h; i++) {
            for (int j = 0; j < image_w; j++) {
                Binimage.at<uint8_t>(i, j) = bin_image[i][j];
            }
        }

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
            draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
            draw_point(border_location_right_optimized[i][0], i, GREEN, Binimage_color);
            // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
            // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
            draw_point(middle_line[i], i, BLUE, Binimage_color);
            draw_point(diffPoint_left[1], diffPoint_left[0], RED, Binimage_color);
            draw_point(diffPoint_right[1], diffPoint_right[0], RED, Binimage_color);
            // draw_point(LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], RED, Binimage_color);
        }

        // ��ʾͼ��
        namedWindow("Binimage_color", WINDOW_NORMAL);
        namedWindow("Binimage_color");
        imshow("Binimage_color", Binimage_color);

        namedWindow("output", WINDOW_NORMAL);
        imshow("output", roadImage);
        resizeWindow("output", 360, 240);
        waitKey(0);
        destroyAllWindows();
    }
	
    if(PVFLAG == 1)
    {
        //��Ƶ����
        VideoCapture cap(VIDEO);
        // �����Ƶ�Ƿ�ɹ���
        if (!cap.isOpened()) {
            std::cerr << "Error opening video file" << std::endl;
            return -1;
        }
        double fps = cap.get(CAP_PROP_FPS);
        double delay = 1000 / fps; 
        Mat frame;
        // ѭ����ȡÿһ֡
        while (cap.read(frame)) {
            printf("------------------------------------------frame:%d------------------------------------------\n", frame_number);
            //����ÿһ֡ΪͼƬ
            string filename = "./frames/frame_" + to_string(frame_number) + ".png";
            imwrite(filename, frame);
            // waitKey(1000);
            // �������ÿһ֡���д���
            dataTrans(frame, image);
            // inRightCircleFlag = 2;
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
                draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
                draw_point(border_location_right_optimized[i][0], i, PINK, Binimage_color);
                // draw_point(border_location_left[i][1], border_location_left[i][0], GREEN, Binimage_color);
                // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
                draw_point(middle_line[i], i, BLUE, Binimage_color);
                draw_point(diffPoint_left[1], diffPoint_left[0], RED, Binimage_color);
                draw_point(diffPoint_right[1], diffPoint_right[0], RED, Binimage_color);
                // draw_point(LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], RED, Binimage_color);
            }
            // ��ʾͼ��
            namedWindow("bin_image");  //, WINDOW_NORMAL
            imshow("bin_image", Binimage_color);
            // resizeWindow("bin_image", 360, 240);
            // namedWindow("output", WINDOW_NORMAL);
            // imshow("output", frame);
            // resizeWindow("output", 360, 240);
            // �ȴ�һ��ʱ�䣬����ESC���˳�
            // if (waitKey(100) == 27) {
            //     break;
            // }
            // �����̰���
            int key = waitKey((int)delay) & 0xFF;
            if (key == 27) { // �� 'ESC' �˳�
                break;
            } else if (key == 'd' || key == 'D') { // �� 'D' ���
                delay /= 2;
                if (delay < 1) delay = 1; // ��С�ӳٲ���С��1����
            } else if (key == 'a' || key == 'A') { // �� 'A' ����
                delay *= 2;
            }
            
        }
    }    
	return 0;
}
