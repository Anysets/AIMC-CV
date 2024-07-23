#include "image_process.h"
#include <fstream>
#include <vector>
// #include <stdlib.h>
#define _CRT_SECURE_NO_WARNINGS
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
#define ROAD  "E:\\Files\\OpenCV++\\test\\test\\163.png"
// #define VIDEO  "C:\\Users\\wqg20\\Downloads\\23.mp4"
#define VIDEO "/home/anysets/Documents/Video/rc.mp4"
#define PVFLAG 1
// #define ROAD "E:\\Files\\OpenCV++\\test.BMP"
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 94
//ע�����border_length����̫�࣬�ᵼ������Խ��
#define border_length 60

uint8_t image_center = image_w / 2;

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
int32_t start_point_l[2] = { 0 };
int32_t start_point_r[2] = { 0 };
int32_t start_line = image_h - 2;  //���һ�б������ˣ�����һ��
int8_t getStartPoint()
{
	uint8_t left_flag = 0;
	uint8_t right_flag = 0;
	//��ʼ��
	start_point_l[0] = 0;
	start_point_l[1] = 0;
	start_point_r[0] = 0;
	start_point_r[1] = 0;
	//���м�����Ѱ��
	for (uint32_t i = image_center; i > 1; i--)  //i>1��ԭ����0��1�����Ѿ���ȫ����0�ˣ�����i=0��1ʱ���������Խ��
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i - 1] == 0 && bin_image[start_line][i - 2] == 0)
		{
			start_point_l[0] = start_line; //��
			start_point_l[1] = i-1; //��
			left_flag = 1;
			break;
		}
	}
	//���м�����Ѱ��
	for (uint32_t i = image_center; i < image_w - 2; i++)  //i < image_w-2��ԭ��ͬ��
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i + 1] == 0 && bin_image[start_line][i + 2] == 0)
		{
			start_point_r[0] = start_line; //x����
			start_point_r[1] = i+1;
			right_flag = 1;
			break;
		}
	}
	if (left_flag && right_flag) { return 1; }
	else { return 1; }
}

//������rewrite
//��¼��߽������
int32_t border_location_left[border_length][2];  //��¼�߽�λ��
int32_t border_count_left = 0;
int32_t growth_direction_left[border_length];  //��¼��������
int32_t growth_count_left = 0;
//��¼�ұ߽������
int32_t border_location_right[border_length][2];  //��¼�߽�λ��
int32_t border_count_right = 0;
int32_t growth_direction_right[border_length];  //��¼��������
int32_t growth_count_right = 0;
void neighborSearch()
{
    memset(border_location_left, 0, sizeof(border_location_left));
    memset(border_location_right, 0, sizeof(border_location_right));

    int32_t neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
    int32_t neighbor_right[9][2] = {{-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}};
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
                start_point_l[0] = start_point_l[0] + neighbor_left[i+1][0];
                start_point_l[1] = start_point_l[1] + neighbor_left[i+1][1];
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

                border_location_right[border_count_right][0] = start_point_r[0] + neighbor_right[i+1][0];
                border_location_right[border_count_right][1] = start_point_r[1] + neighbor_right[i+1][1];
				// bin_image[start_point_r[0] + neighbor_right[i+1][0]][start_point_r[1] + neighbor_right[i+1][1]] = 151;
                start_point_r[0] = start_point_r[0] + neighbor_right[i+1][0];
                start_point_r[1] = start_point_r[1] + neighbor_right[i+1][1];
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
	// 	printf("%d\n", growth_direction_left[i]);
	// }

    // ���border_location_left
    // for (int i = 0; i < border_length; i++)
    // {
    //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
    // }
}

// �洢�Ż������߽磬ÿ�����ض��зݣ�ע��˴���ȡֵ��ΧΪ-127��127��������ѹ�����ͼ��Χ��
// ����һ��Ҫע�⣬�˴��Ƕ�ά���飬���ǵ�0��ֵ��image_w����1��ֵ�ǵ�ǰ����ԭ�߽����������ǰֵ����������image_h
uint8_t border_location_left_optimized[image_h][2] = { 0 };
uint8_t border_location_right_optimized[image_h][2] = { 0 };
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
    // for (int i = 0; i < image_h; i++)
    // {
    //     // bin_image[i][border_location_left_optimized[i][0]] = 152;  // ���߽��ǳ���
    //     // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
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
    //     // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
    // }   
}

// ��ȡ����
int32_t middle_line[border_length] = { 0 }; 
void getMiddleLine()
{
    memset(middle_line, 0, sizeof(middle_line));
    for (int i = 0; i < image_h; i++)
    {
        middle_line[i] = (border_location_left_optimized[i][0] + border_location_right_optimized[i][0]) / 2;
    }
}

// ��ȡ�Ƕ�
double angle_deg = 0;
void getAngle(int8_t opposite, uint8_t adjacent)
{
    // opposite = abs(opposite);
    // printf("opposite=%d\n", opposite);
    // printf("adjacent=%d\n", adjacent);
    double angle_rad = atan2(opposite, adjacent);
    angle_deg = angle_rad * (180.0 / 3.14);
    // printf("angle_deg=%f\n", angle_deg);
}

// ��ȡƫ��ֵ
int8_t deviation_number;
void getDeviation()
{
    int count = 0;
    deviation_number = 0;
    for (int i = 40; i < 55; i++)
    {
        deviation_number += middle_line[i] - (image_w/2);
        count++;
    }
    deviation_number = deviation_number / count;
    printf("deviation_number=%d\n", deviation_number);
    getAngle(deviation_number, 15);
}




//Ѱ�ҹյ�
uint8_t diffPoint_left[2] = {0}; // ��¼�յ�
uint8_t diffPoint_right[2] = {0}; // ��¼�յ�
uint8_t findTurningPoint()
{
    memset(diffPoint_left, 0, sizeof(diffPoint_left));
    memset(diffPoint_right, 0, sizeof(diffPoint_right));
    uint8_t flag = 0;
    uint8_t left_flag = 0;
    uint8_t right_flag = 0;
    uint8_t diffPoint_temp[2] = {};
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
        if ((growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // �յ��ڶ�����
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
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 1))
            {
                // printf("+2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2]);
                // printf("+1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1]);
                // printf("-1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1]);
                // printf("-2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2]);
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
    growth_direction_left[border_location_right_optimized[diffPoint_temp[0]][1]];

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
    if (border_location_right_optimized[diffPoint_temp[0]][1] > 1 && border_location_right_optimized[diffPoint_temp[0]][1] < border_length-2)  // ��ֹ����Խ��
    {
        if ((growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // �յ��ڶ�����
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
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 1))
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
    // 3��������ҹյ�߶�������ȡ����·��ĵ���Ϊ�յ����
    if (diffPoint_right[0] - diffPoint_left[0] > 25)
    {
        left_flag = 0;
    }
    else if (diffPoint_left[0] - diffPoint_right[0] > 25)
    {
        right_flag = 0;
    }

    // 4���˴�Ϊ��ֹ���ж�д����Ӱ������ɾ��
    printf("diffPoint_left[0]: %d\n", diffPoint_left[0]);
    printf("diffPoint_right[0]: %d\n", diffPoint_right[0]);
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
void FillLine(uint8_t line[image_h][2], uint8_t x1, uint8_t y1, uint8_t x2, uint8_t y2) {
    double slope, intercept;
    // printf("x1:%d y1:%d x2:%d y2:%d\n", x1, y1, x2, y2);
    two_slope(x1,y1,x2,y2,&slope,&intercept);
    for (int i = y1; i <= y2; i++) {
        line[i][0] = (i - intercept) / slope;
        // printf("line[%d][0]:%d\n", i, line[i][0]);
    }
}
// ---------------------------------------------------------------------------------------------

// ��¼��ǰ״̬�ı�־
// bool isInRightCircle = false;
// bool isInCrossing = false;
uint8_t inRightCircleFlag = 0;
uint8_t inCrossingFlag = 0;

uint8_t status3_bwline[2] = {0};  // black white line
long frame1 = 0;
uint32_t amount_fps[2] = {0, 0};  // ���ڼ�¼״̬0ʱ�����ҵ�5���ҹյ�

uint8_t right_point_row1[2] = {0};  // �����¼���ǳ�̬�µ����
// �ҳ��һ����������
uint8_t LeftPoint_RightCircle[2];
void FindRightCircleLeftPoint(uint8_t flag)
{
    memset(LeftPoint_RightCircle, 0, sizeof(LeftPoint_RightCircle));
    if (flag == 0)
    {
        diffPoint_right[0] = right_point_row1[0];
        diffPoint_right[1] = right_point_row1[1];
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
    int32_t neighbor[9][2] = {{1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 1}, {1, 1}, {1, 0}, {1, -1}};
    for (int j = diffPoint_right[0]-1; j > 1; j--)
    {
        for (int i = 0; i < 8; i++)
        {
            if (bin_image[LeftPoint_RightCircle[0] + neighbor[i][0]][LeftPoint_RightCircle[1] + neighbor[i][1]] == 255 && bin_image[LeftPoint_RightCircle[0] + neighbor[i+1][0]][LeftPoint_RightCircle[1] + neighbor[i+1][1]] == 0)
            {       
                if (LeftPoint_RightCircle[1] >= LeftPoint_RightCircle[1] + neighbor[i+1][1])
                {
                    LeftPoint_RightCircle[0] = LeftPoint_RightCircle[0] + neighbor[i+1][0];
                    LeftPoint_RightCircle[1] = LeftPoint_RightCircle[1] + neighbor[i+1][1];
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

// �ж���Բ��
void rightCircleHandle(uint8_t TurningPointFlag)
{
    printf("amount_fps[0]: %d\n", amount_fps[0]);
    printf("amount_fps[1]: %d\n", amount_fps[1]);
    if (inRightCircleFlag == 0 && TurningPointFlag == 2 && inCrossingFlag == 0)
    {
        if (amount_fps[0] == 0)
        {
            
            amount_fps[0]++;
            amount_fps[1] = frame1;
        }
        else
        {
            if (frame1 - amount_fps[1] > 1)
            {
                amount_fps[0] = 0;
            }
            else if (amount_fps[0] >= 5)
            {
                amount_fps[0] = 0;
                inRightCircleFlag = 1;  // ����״̬1
                getStartPoint();  //��ʺɽ����ô�ѳɵġ�
                right_point_row1[0] = start_point_r[0];
                right_point_row1[1] = start_point_r[1];
            }
            else
            {
                amount_fps[0]++;
                amount_fps[1] = frame1;
            }
        }
    }
    else if (inRightCircleFlag == 1)
    {
        if(TurningPointFlag == 2 || TurningPointFlag == 3)
        {
            FindRightCircleLeftPoint(1);
            printf("LeftPoint_RightCircle[1]:%d\n", LeftPoint_RightCircle[1]);
            printf("LeftPoint_RightCircle[0]:%d\n", LeftPoint_RightCircle[0]);
            FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
        }
        else
        {
            FindRightCircleLeftPoint(0);
            FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
            if (amount_fps[0] == 0)
            {
                amount_fps[0]++;
                amount_fps[1] = frame1;
            }
            else
            {
                if (frame1 - amount_fps[1] > 1)
                {
                    amount_fps[0] = 0;
                }
                else if (amount_fps[0] >= 3)
                {
                    amount_fps[0] = 0;
                    inRightCircleFlag = 2;  // ����״̬2
                }
                else
                {
                    amount_fps[0]++;
                    amount_fps[1] = frame1;
                }
            }
        }
    }
    else if (inRightCircleFlag == 2)
    {
        FindRightCircleLeftPoint(0);
        FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_r[1] < 86)
        {
            inRightCircleFlag = 3;
        }
    }
    else if (inRightCircleFlag == 3)
    {
        // ����������10����ǰ��ֱ����Ϊ��ʼ�㣨��ʱ�Ҳ��Ѿ����ߣ�
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        // printf("start_point_l[1]: %d\n", start_point_l[1]);
        uint8_t start_point = image_h - 2;  // ��¼��ʼ��
        uint8_t diffPoint_rightUp[2] = {0, 90};
        // && start_point_r[1] > 90
        if (start_point_l[1] < 10)
        // if (false)
        {
            printf("In if\n");
            // ʹ���Ż�������߽�����Ѱ�ҹյ�
            bool found_flag = false;  // ָʾ�Ƿ��ҵ��յ㣬����Ƕ��ѭ��
            
            for (uint8_t i = start_point; i > 1; i--)
            {
                for (uint8_t j = border_location_left_optimized[i][0]+1; j < image_w - 15; j++)
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
                        printf("border_location_left_optimized[%d][0]: %d\n",i , border_location_left_optimized[i][0]);
                        
                        found_flag = true;
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
            // �ҵ��յ��ʼ����
            // printf("%d\n", start_point);
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 10, start_point_l[0]);
        }
        else
        {
            // ����ʹ���ұ߽�Ѱ����ʼ��
            printf("In else\n");
            for (uint8_t i = image_h - 2; i > 1; i--)
            {
                if (border_location_right_optimized[i][0] >= 90)
                {
                    start_point = i;
                    // printf("start_point: %d\n", start_point);
                    break;
                }
            }
            // ʹ���Ż�������߽�����Ѱ�ҹյ�
            bool found_flag = false;  // ָʾ�Ƿ��ҵ��յ㣬����Ƕ��ѭ��
            for (uint8_t i = start_point; i > 1; i--)
            {
                for (uint8_t j = border_location_left_optimized[i][0]+1; j < image_w - 15; j++)
                {
                    // printf("diffPoint_rightUp: %d, %d\n", i, j);
                    if (bin_image[i][j] == 0)
                    {
                        printf("diffPoint_rightUp: %d, %d\n", i, j);
                        // printf("border_location_left_optimized[i][0]: %d\n", border_location_left_optimized[i][0]);
                        found_flag = true;
                        diffPoint_rightUp[0] = i;
                        diffPoint_rightUp[1] = j;
                        break;
                    }
                }
                // if (found_flag)
                // {
                    
                //     break;
                // }
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
            // �ҵ��յ��ʼ����
            // printf("%d\n", start_point);
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], start_point_l[1], start_point_l[0]);
        }
        // �ҵ��յ��������������
        // ������߽磬ֱ���������Ҳ�ĵ�
        /*
            [2]  [1]  [0/8]

            [3]  [ ]  [7]

            [4]  [5]  [6]
        */
        int32_t neighbor_left[9][2] = {{-1, 1}, {-1, 0}, {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}};
        memset(border_location_left, 0, sizeof(border_location_left));
        border_count_left = 0;
        start_point_l[0] = diffPoint_rightUp[0];
        start_point_l[1] = diffPoint_rightUp[1];
        for (uint8_t j = 0; j < border_length-1; j++)
        {
            for (uint8_t i = 0; i < 8; i++)
            {
                if (bin_image[start_point_l[0] + neighbor_left[i][0]][start_point_l[1] + neighbor_left[i][1]] == 255 && bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] == 0)
                {
                    start_point_l[0] = start_point_l[0] + neighbor_left[i+1][0];
                    start_point_l[1] = start_point_l[1] + neighbor_left[i+1][1];
                    border_location_left[border_count_left][0] = start_point_l[0];
                    border_location_left[border_count_left][1] = start_point_l[1];
                    border_count_left++;
                    break;
                }
            }
            if (start_point_l[1] > 91)
            {
                printf("break!!!\n");
                // printf("start_point_l[1]: %d\n", start_point_l[1]);
                // printf("start_point_l[0]: %d\n", start_point_l[0]);
                break;
            }
        }
        // ���border_location_left
        // for (int i = 0; i < border_length; i++)
        // {
        //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
        // }
        optimizeBorder();
        

        // ������һ��״̬
        // // ����1:���Ϲյ���ʧ
        // printf("diffPoint_rightUp[0]: %d\n", diffPoint_rightUp[0]);
        // printf("diffPoint_rightUp[1]: %d\n", diffPoint_rightUp[1]);
        // if (diffPoint_rightUp[0] > 50 && diffPoint_rightUp[1] < image_w/2)  // ���Ϲյ���45������
        // {
        //     inRightCircleFlag = 4;
        //     // ��¼��ǰ֡��
        //     amount_fps[1] = frame1;
        //     printf("Status 4\n");
        // }
        
        // ����2�����������ʧ�ٳ���  С��״̬��
        // 
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (status3_bwline[0] == 0)
        {
            if (start_point_l[1] > 5)
            {
                status3_bwline[0] = 1;
                status3_bwline[1] = frame1;
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
                inRightCircleFlag = 4;  // ����״̬4
            }
        }
    }
    else if (inRightCircleFlag == 4)
    {
        uint8_t right_point = 0;
        uint8_t diff = 0;
        // û�дﵽָ��֡�������������һ״̬
        if (frame1 - amount_fps[1] > 10)
        // if (true)
        {
            // ����1
            // for (uint8_t i = image_h-10; i > 10; i++)
            // {
            //     diff = border_location_left_optimized[i-1][0] - border_location_left_optimized[i][0];
            //     if (diff > 10 || diff < -10)
            //     {
            //         FillLine(border_location_left_optimized, i, border_location_left_optimized[i][0], 0, image_w/2);
            //         inRightCircleFlag = 5;
            //         break;
            //     }
            // }

            // ����2
            // // �ùյ��ж�
            // if (TurningPointFlag == 1 || TurningPointFlag == 3)
            // {
            //     FillLine(border_location_left_optimized, diffPoint_left[1], diffPoint_left[0], image_w/2, 0);
            //     inRightCircleFlag = 5;
            //     // break;
            // }

            // ����3
            // ����߽����ҵ��ж�
            // uint8_t rightest_w = 0;
            // for (uint8_t i = image_h-2; i > 1; i--)
            // {
            //     printf("border_location_left_optimized[%d][0]: %d\n", i, border_location_left_optimized[i][0]);
            //     if (border_location_left_optimized[i][0] > rightest_w)
            //     {
                    
            //         rightest_w = border_location_left_optimized[i][0];
            //     }
            // }
            // printf("rightest_w: %d\n", rightest_w);
            // if (rightest_w < image_w - 3)
            // {
            //     inRightCircleFlag = 5;
            // }

            // ����4  ѧ����code
            for (int i = 5; i < 100; i++) 
            {
                if (border_location_left[i][1] >= border_location_left[i - 1][1] && border_location_left[i][1] > border_location_left[i + 1][1]
                    && border_location_left[i][1] > border_location_left[i - 10][1]) {   //������� һ��Ҫ�е��ں�
                    uint8_t row = border_location_left[i][0];
                    uint8_t col = border_location_left[i][1];
                    // printf("status4:  row: %d, col: %d\n", row, col);
                    //��һ����Ϊ���½ǹ̶��� �ڶ�����Ϊ���Ϲ̶�
                    FillLine(border_location_left_optimized, image_w - 5, 5, col, row);  //ѹ��
                    break;
                }
            }

            // �������㶪�ˣ�������һ״̬
            getStartPoint();  //��ʺɽ����ô�ѳɵġ�
            if (start_point_l[1] < 3)
            {
                inRightCircleFlag = 5;
            }
        }
        
    }
    else if (inRightCircleFlag == 5)
    {
        // for (uint8_t i = image_h-20; i > 20; i++)
        // {
        //     uint8_t diff = 0;
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
        FillLine(border_location_left_optimized, image_w - 1, image_h / 3, image_w / 3, image_h - 1);

        // ��������
        double slope_stick = 0;
        double slope_find = 0;
        double intercept_stick = 0;
        double intercept_find = 0;
        two_slope(image_w - 1, 0, 0, image_h - 1, &slope_stick, &intercept_stick);
        two_slope(border_location_left_optimized[image_h-5][1], border_location_left_optimized[image_h-5][0], border_location_left_optimized[image_h/2][1], border_location_left_optimized[image_h/2][0], &slope_find, &intercept_find);    
        if (slope_find - slope_stick> 2)
        {
            inRightCircleFlag = 6;
            // printf("slope_find - slope_stick: %f\n", slope_find - slope_stick);
        }
    }
    else if (inRightCircleFlag == 6)
    {
        FillLine(border_location_left_optimized, image_h-5, 0, 0, image_w/2);
        // printf("border_location_left_optimized[image_h-2][0]: %d\n", border_location_left_optimized[image_h-2][0]);
        if (border_location_left_optimized[image_h-2][0] > 10)
        {
            inRightCircleFlag = 7;
        }
    }
    else if (inRightCircleFlag == 7)
    {
        FindRightCircleLeftPoint(0);
        FillLine(border_location_right_optimized, LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], diffPoint_right[1], diffPoint_right[0]);
        // ���½�����㣬�˳�����
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_r[1] < 90)
        {
            inRightCircleFlag = 0;
            // printf("start_point_r[1]: %d\n", start_point_r[1]);
        }
    }
}

// ʶ��ʮ�ֺ�Բ������Ҫ��������ʮ�����ֻʶ��һ���յ�ʱ����һ���Ƕ��ߵ�
// ��Բ���뻷��������ֻʶ��һ���յ㣬����һ���ǲ����ߵ�
uint8_t up_right_point[2] = {0};
uint8_t up_left_point[2] = {0};
void findUpTurningPoint(uint8_t start_l, uint8_t start_r)
{
    
    // �����Ϲյ�
    for (uint8_t i = start_l; i > 1; i--)
    {
        if ((border_location_left_optimized[i-1][0] - border_location_left_optimized[i][0] < 5) &&
            (border_location_left_optimized[i][0] - border_location_left_optimized[i+2][0] > 10))
        {
            // printf("Found Left Up Turning Point\n");
            up_left_point[0] = i;
            up_left_point[1] = border_location_left_optimized[i][0];
            break;
        }
    }
    // �����Ϲյ�
    for (uint8_t i = start_r; i > 1; i--)
    {
        if ((border_location_left_optimized[i][0] - border_location_left_optimized[i-1][0] < 5) &&
            (border_location_right_optimized[i+2][0] - border_location_right_optimized[i][0] > 10))
        {
            // printf("Found right Up Turning Point\n");
            up_right_point[0] = i;
            up_right_point[1] = border_location_right_optimized[i][0];
            break;
        }
    }
}

// ʶ��ʮ��
void crossingHandle(uint8_t TurningPointFlag)
{
    if (inRightCircleFlag == 0 && (TurningPointFlag == 1 || TurningPointFlag == 3) && inCrossingFlag == 0)
    {
        getStartPoint();  //��ʺɽ����ô�ѳɵġ�
        if (start_point_r[1] > 90 && TurningPointFlag == 1)  // ֻʶ����յ�ʱ������Ƿ�����㣬�����ˣ�����״̬1
        {
            inCrossingFlag = 1;
        }
        else if(TurningPointFlag == 3)  // �����յ㶼ʶ�𵽣�����״̬1
        {
            inCrossingFlag = 1;
        }
    }
    else if (inCrossingFlag == 1)
    {
        // ��ʹ������ʮ�֣�Ҳ�п���ֻʶ����һ���յ㣬����Ҫ�����ж�
        // ��������յ㶼�����ҵ��ģ������ҵ��Ϸ��������յ�
        if (TurningPointFlag == 3)
        {
            findUpTurningPoint(diffPoint_left[0], diffPoint_right[0]);

            // ����
            FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], diffPoint_left[1], diffPoint_left[0]);
            FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], diffPoint_right[1], diffPoint_right[0]);
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
        // printf("start_point_l[0]: %d\n", start_point_l[0]);
        // printf("start_point_r[0]: %d\n", start_point_r[0]);
        // printf("up_left_point[1]: %d\n", up_left_point[1]);
        // printf("up_left_point[0]: %d\n", up_left_point[0]);
        // printf("up_right_point[1]: %d\n", up_right_point[1]);
        // printf("up_right_point[0]: %d\n", up_right_point[0]);
        findUpTurningPoint(start_point_l[0]-1, start_point_r[0]-1);
        // ����  ע������ģ�25, 58����75, 58�������ʵ�������� 
        FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], 25, 58);
        FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], 75, 58);
        // ������಻�����ˣ��˳�ʮ��·��
        
        if (start_point_l[1] > 4 && start_point_r[1] < 90)
        {
            inCrossingFlag = 0;
        }
    }
}

// Ԫ���ж�������
void elementsHandle()
{
    uint8_t TurningPointFlag = 0;
    TurningPointFlag = findTurningPoint();
    printf("TurningPointFlag:%d\n", TurningPointFlag);
    printf("inRightCircleFlag:%d\n", inRightCircleFlag);
    printf("inCrossingFlag:%d\n", inCrossingFlag);
    crossingHandle(TurningPointFlag);
    rightCircleHandle(TurningPointFlag);
    
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
        memset(border_location_left_optimized, 0, sizeof(border_location_left_optimized));
        memset(border_location_right_optimized, 0, sizeof(border_location_right_optimized));
	    optimizeBorder();
        elementsHandle();
        getMiddleLine();
        getDeviation();
        
    }
    else
    {
        // printf("Find start point failed.\n");
    }

    frame1++;
    // printToCSV(bin_image, "output.csv");
    // usleep(20000);
}

int main(){
    if (PVFLAG == 0)
    {
        inRightCircleFlag = 3;
        inCrossingFlag = 0;
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
            //draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
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
        Mat frame;
        // ѭ����ȡÿһ֡
        while (cap.read(frame)) {
            printf("------------------------------------------frame:%d------------------------------------------\n", frame1);
            //����ÿһ֡ΪͼƬ
            string filename = to_string(frame1) + ".png";
            imwrite(filename, frame);
            // waitKey(1000);
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
                draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
                draw_point(border_location_right_optimized[i][0], i, GREEN, Binimage_color);
                // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
                //  draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
                draw_point(middle_line[i], i, BLUE, Binimage_color);
                draw_point(diffPoint_left[1], diffPoint_left[0], RED, Binimage_color);
                draw_point(diffPoint_right[1], diffPoint_right[0], RED, Binimage_color);
                // draw_point(LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], RED, Binimage_color);
            }
            // ��ʾͼ��
            namedWindow("bin_image", WINDOW_NORMAL);
            imshow("bin_image", Binimage_color);
            resizeWindow("bin_image", 360, 240);
            // namedWindow("output", WINDOW_NORMAL);
            // imshow("output", frame);
            // resizeWindow("output", 360, 240);
            // �ȴ�һ��ʱ�䣬����ESC���˳�
            if (waitKey(20) == 27) {
                break;
            }
            
        }
    }
    
    
	return 0;
}
