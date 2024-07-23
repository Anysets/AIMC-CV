#define _CRT_SECURE_NO_WARNINGS
#include "image_process.h"

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
#define start_line image_h - 2
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
        //printf("Error opening file!\n");
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
//����ֵ����
double double_abs(double value)
{
    if (value >= 0) return value;
    else return -value;
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
uint32_t black_points;
void turn_to_bin(void) {
    black_points = 0;
    uint8_t image_threshold = otsuThreshold(image);
    // ����ԭʼͼ�����飬���ݶ�̬��ֵ���ж�ֵ��
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // ���ԭʼͼ��Ҷ�ֵ������ֵ������Ϊ��ɫ����
            }
            else {
                image[i][j] = 0; // ���ԭʼͼ��Ҷ�ֵС�ڵ�����ֵ������Ϊ��ɫ����
                black_points++;
            }
        }
    }
    //printf("��ɫ�ĵ�%d\n", black_points);
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
#define threshold_max	(255 * 8) // �˲����ɸ����Լ����������
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
        image[39][i] = 0;
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
//���������
uint8_t get_single_point(uint8_t end_row) {
    for (uint8_t row = 37; row > end_row; row--) {
        if (bin_image[row][47] == 255 && bin_image[row - 1][47] == 0) {
            start_right[0] = row;
            start_right[1] = 47;
            start_left[0] = row;
            start_left[1] = 47;
            break;
        }
    }
    return 0;
}
//����������
#define search_amount 80   //����������������ĵ�
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
                //printf("��%d�Σ�L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
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
                //printf("��%d�Σ�R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
//���߰�����
uint8_t single_amount = 20; //������������ĵ�
void single_search(void) {
    L_search_count = R_search_count = 0;    //�ѵ��ĵ�������
    left_points[0][0] = start_left[0];  //y row
    left_points[0][1] = start_left[1];  //X col
    int8_t L_curr_row = start_left[0];   //��ʼ��row����
    int8_t L_curr_col = start_left[1];  //��ʼ��col����
    for (int j = 1; j < single_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[L_curr_row + L_dir_row[i]][L_curr_col + L_dir_col[i]] == 0)
                && (bin_image[L_curr_row + L_dir_row[(i + 1) & 7]][L_curr_col + L_dir_col[(i + 1) & 7]] == 255)) {
                //�ɺڵ���
                L_search_count++;
                L_curr_row += L_dir_row[(i + 1) & 7];
                L_curr_col += L_dir_col[(i + 1) & 7];
                left_points[L_search_count][0] = L_curr_row;   //�洢�׵�
                left_points[L_search_count][1] = L_curr_col;   //�洢�׵� 
                dir_left[L_search_count] = i;
                //printf("��%d�Σ�L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
                break;
            }
        }
    }
    right_points[0][0] = start_right[0];  //y row
    right_points[0][1] = start_right[1];  //X col
    int8_t R_curr_row = start_right[0];   //��ʼ��row����
    int8_t R_curr_col = start_right[1];  //��ʼ��col����
    for (int j = 1; j < single_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[R_curr_row + R_dir_row[i]][R_curr_col + R_dir_col[i]] == 0)
                && (bin_image[R_curr_row + R_dir_row[(i + 1) & 7]][R_curr_col + R_dir_col[(i + 1) & 7]] == 255)) {
                R_search_count++;
                R_curr_row += R_dir_row[(i + 1) & 7];
                R_curr_col += R_dir_col[(i + 1) & 7];
                right_points[R_search_count][0] = R_curr_row;   //�洢�׵�
                right_points[R_search_count][1] = R_curr_col;   //�洢�׵�
                dir_right[R_search_count] = i;
                //printf("��%d�Σ�R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
uint8_t single_left;    //1�Ƕ���
uint8_t single_right;   //1�Ƕ���
uint8_t single_center;  //1�Ƕ���
void single_element() {
    single_left = single_right = single_center =0;
    uint8_t center_count = 0;
    for (center_count = image_h - 2; center_count > single_amount; center_count--) {
        if (bin_image[center_count][image_w / 2] == 255 && bin_image[center_count - 1][ image_w / 2] == 0) {    //�׵���
            break;
        }
    }
    if (center_count == single_amount) {
        single_center = 1;
        printf("�м䵥�߶���\n");
    }
    for (uint8_t count = 0; count < single_amount - 1; count++) {
        if (right_points[count][0] > right_points[count + 1][0]) {
            uint8_t row = right_points[count][0];
            uint8_t col = right_points[count][1];
            if(bin_image[row - 5][col - 5] == 255){
                right_points[single_amount - 1][0] = right_points[count][0];    //����
                right_points[single_amount - 1][1] = right_points[count][1];
                single_right = 1;
                printf("right\n");
                break;
            }
        }
    }
    for (uint8_t count = 0; count < single_amount - 1; count++) {
        if (left_points[count][0] > left_points[count + 1][0]) {
            uint8_t row = left_points[count][0];
            uint8_t col = left_points[count][1];
            if (bin_image[row - 5][col + 5] == 255) {
                left_points[single_amount - 1][0] = left_points[count][0];
                left_points[single_amount - 1][1] = left_points[count][1];
                single_left = 1;
                printf("left\n");
                break;
            }
        }
    }
}
//��ȡ����
uint8_t left_line[image_h][2];  //�������� row col
uint8_t right_line[image_h][2];    //�������� row col
uint8_t center_line[image_h][2];   //�������� row col
uint8_t width[image_h]; //�������
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
//��С���˷�����ֱ����ϵ�б�ʺͽؾ�
void LeastSquaresSlope(uint8_t line[image_h][2], uint8_t start, uint8_t end, double* slope, double* intercept) {
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
    uint8_t size = end - start + 1;
    for (int i = start; i <= end; i++) {
        sumX += line[i][1];
        sumY += line[i][0];
        sumXY += line[i][1] * line[i][0];
        sumX2 += line[i][1] * line[i][1];
    }
    double xMean = sumX / size;
    double yMean = sumY / size;
    *slope = (sumXY - size * xMean * yMean) / (sumX2 - size * xMean * xMean);
    *intercept = (sumY - *slope * sumX) / size;
}
void fill_line(uint8_t line[image_h][2], uint8_t start_row, uint8_t start, uint8_t end) {
    double slope, intercept;
    LeastSquaresSlope(line, start_row + start, start_row + end, &slope, &intercept);
    for (int i = start_row; i > 5; i--) {
        line[i][1] = (line[i][0] - intercept) / slope;
    }
}
//����������б��
void two_slope(double x1, double y1, double x2, double y2, double* slope, double* intercept) {
    *slope = (y2 - y1) / (x2 - x1);
    *intercept = y1 - (*slope * x1);
}
//�������㲹��
void two_line(uint8_t line[image_h][2], uint8_t x1, uint8_t y1, uint8_t x2, uint8_t y2) {
    double slope, intercept;
    two_slope(x1,y1,x2,y2,&slope,&intercept);
    for (int i = y1; i >= y2; i--) {
        line[i][1] = (line[i][0] - intercept) / slope;
    }
}
// ������ϵķ���
double variance(uint8_t line[image_h][2], uint8_t start, uint8_t end,uint8_t step, double* slope, double* intercept) {
    LeastSquaresSlope(line, start, end, slope, intercept);
    double sum = 0;
    int size = end - start + 1;
    if (*slope == 1) {
        return 0;
    }
    else{
        for (int i = start; i <= end; i+=step) {
            double yFit = *slope * line[i][1] + *intercept;
            sum += pow(line[i][0] - yFit, 2);
        }
        return sum * step/ size;
    }
}
//���������ж�
uint8_t dir_judge(uint8_t dir[search_amount], uint8_t direction, uint8_t count) {
    uint8_t flag = 0;
    for (uint8_t i = count; i < count + 15; i++) {
        if (dir[i] == direction) {
            flag++;
        }
        if (flag >= 12) {

            return 1;
        }
    }
    return 0;
}
//����ɨ��
uint8_t col_middle = image_w / 2;
uint8_t min_right_col = 188;
uint8_t min_right_row;
uint8_t row_right_search(uint8_t start_row,uint8_t start_col) { 
    min_right_col = 188;
    for (uint8_t row = start_row; row > 5; row--) {
        if (bin_image[row][start_col] == 255 && bin_image[row - 1][start_col] == 0) {             
            for (uint8_t i = row - 1; i >= 5; i--) {
                for(int j = col_middle;j < 90;j++){
                    if (bin_image[i][j] == 255 && bin_image[i - 1][j] == 0) {
                        if (min_right_col > j){
                            min_right_col = j;
                            min_right_row = i;
                        }
                        break;  //ÿһ��ֻ�ѵ�һ����
                    }
                }
            }
            if (min_right_col > 10) return 1;
            else return 0;
        }
    }
}
uint8_t max_left_row;
uint8_t max_left_col = 0;
uint8_t row_left_search(uint8_t start_row, uint8_t start_col) {
    max_left_col = 0;
    for (uint8_t row = start_row; row > 5; row--) { //��ʼ����ɨ�ҵ���һ���ڵ� �پ�׼��λ
        if (bin_image[row][start_col] == 255 && bin_image[row - 1][start_col] == 0) {  //�׺� �ҵ���
            for (uint8_t i = row - 1; i >= 5; i--) {
                for (uint8_t j = col_middle; j > 3; j--) {
                    if (bin_image[i][j] == 255 && bin_image[i - 1][j] == 0) {
                        if (max_left_col < j) {
                            max_left_col = j;
                            max_left_row = i;
                        }
                        break;  //ÿһ��ֻ�ѵ�һ����
                    }
                }
            }
            if (max_left_col > 10) return 1;
            else return 0;
        }
    }
}
uint8_t element_start_line = 38; //���������
uint8_t element_end_line = 10;
uint8_t circle_right_row;
uint8_t circle_left_row;
uint8_t cross_right_row;
uint8_t cross_left_row;
void get_turning_point(void) {
    //��־λ����
    circle_right_row = circle_left_row = cross_right_row = cross_left_row =  0;
    for (int i = 0; i < R_search_count; i++) {  //�ҹյ��ж�
        double slope, intercept;
        uint8_t turn_row = right_points[i][0];
        if (turn_row <= element_end_line) break;
        if (turn_row >= element_start_line) continue;
        //�������߽���Ԫ���ж�
        if (right_line[turn_row - 2][1] == 91 && right_line[turn_row - 4][1] == 91 && right_line[turn_row - 6][1] == 91) {
            double err = variance(left_line, turn_row - 10, image_h - 5,1,&slope,&intercept);
            if (dir_judge(dir_right, 6, i)) {
                if(double_abs(err) <= 0.3 ){
                    circle_right_row = turn_row;
                    break;
                }
                else {
                    cross_right_row = turn_row;
                    break;
                }
            }   
        }   
    }   
    for (int i = 0; i < L_search_count; i++) {
        double slope, intercept;
        uint8_t turn_row = left_points[i][0];
        if (turn_row <= element_end_line) break;
        if (turn_row >= element_start_line) continue;
        //�������߽���Ԫ���ж�
        if (left_line[turn_row - 2][1] == 2 && left_line[turn_row - 4][1] == 2 && left_line[turn_row - 6][1] == 2) {
            double err = variance(right_line, turn_row - 10, image_h - 5, 1,&slope, &intercept);
            //printf("err= %f\n",err);
            if (dir_judge(dir_left, 6, i)) {
                if (double_abs(err) <= 0.3) {
                    circle_left_row = turn_row;
              
                    break;
                }
                else {
                    cross_left_row = turn_row;
                    break;
                }
            }
        }

    }
}
//�ж�����ʮ�ֻ����ұ�ʮ�֣��ұ�Ϊ1�����Ϊ2
uint8_t corss_dir(void) {
    turn_to_bin();  //ͼ���ֵ��
    compressimage(); //ͼ��ѹ��
    image_filter();
    image_draw_rectan(bin_image);
    get_start_point(10);
    search();
    uint32_t sum_center_line = 0;
    uint8_t startline = 9;
    uint8_t endline = 5;
    for (int i = startline; i > endline; i--)
    {
        center_line[i][1] = (left_points [i] [1] + right_points[i][1]) >> 1;//������
        sum_center_line += center_line[i][1];
    }        
    float temp_value = (float)sum_center_line / (startline - endline);
    printf("temp_value = %f", temp_value);
    if (temp_value - image_w / 2 > 0) {
        printf("�ұ�ʮ��\n");
        return 1;
    }
    else {
        printf("���ʮ��\n");
        return 2;
    }
    return  temp_value;
}
uint8_t circle_right_flag;
uint8_t circle_right_last_flag;
uint8_t circle_left_flag;
uint8_t circle_left_last_flag;
uint8_t cross_flag;
uint8_t cross_last_flag;
uint8_t zqzq = 0;
void element(void) {
    if (cross_left_row && cross_right_row) {
        uint8_t col_right = right_line[cross_right_row][1];
        uint8_t result_right = row_right_search(cross_right_row, col_right);
        uint8_t col_left = left_line[cross_left_row][1];
        uint8_t result_left = row_left_search(cross_left_row, col_left);
        //if (result_right && result_left) {
        //    two_line(left_line, col_left, cross_left_row, max_left_col, max_left_row);
        //    two_line(right_line, col_right, cross_right_row, min_right_col, min_right_row);
        //}
        corss_dir();
        cross_flag = cross_last_flag = 1;
    }
    else if (circle_left_row) {
        uint8_t col_left = left_line[circle_left_row][1];
        uint8_t result_left = row_left_search(circle_left_row, col_left);
        if (result_left) {
            two_line(left_line, col_left, circle_left_row, max_left_col, max_left_row);
        }
        circle_left_flag = circle_left_last_flag = 1;
    }
    else if (circle_right_row) {
        uint8_t col_right = right_line[circle_right_row][1];
        uint8_t result_right = row_right_search(circle_right_row, col_right);
        if (result_right) {
            two_line(right_line, col_right, circle_right_row, min_right_col, min_right_row);
        }
        circle_right_flag = circle_right_last_flag = 1;
    }
    else {
        cross_flag = circle_left_flag = circle_right_flag = 0;
        zqzq = 0;
    }

    if (circle_right_last_flag == 1 && circle_right_flag == 0) {
        printf("��Բ��");

        //��Բ����־λ
        circle_right_last_flag = 0;
    }
    if (circle_left_last_flag == 1 && circle_left_flag == 0) {
        printf("��Բ��");       
        zqzq = 1;
        //��Բ����־λ
        circle_left_last_flag = 0;
    }
    if (cross_last_flag == 1 && cross_flag == 0) {
        printf("ʮ��");
        //ʮ�ֱ�־λ
        cross_last_flag = 0;
    }
}

//uint8_t single_line[3][2]; //row col [0]������� [1]�����ұ� [2]�����м�
//uint8_t single_left = 0, single_right = 0,single_middle = 0; //1������ 2����δ����
//uint8_t single_end_line = image_h - 10 ;  //���·�����Ѷ�����
//void single_edge(){
//    int mid_row;
//    single_left = 0, single_right = 0;  //����
//    for (int row = image_h - 1; row >= 1; row--) {
//        if (bin_image[row][47] == 255 && bin_image[row - 1][47] == 0 && bin_image[row - 2][47] == 0){  //�׺ں�
//            single_line[2][0] = row;
//            single_line[2][1] = 47;
//            mid_row = row;  //ͼ���м����
//            break;
//        }
//        if (row <= single_end_line){
//            single_middle = 1;  //����
//            mid_row = row;
//        }
//    }
//    for (int row = image_h - 1; row >= mid_row - 5; row--){ //����ѵ��м���5������
//        if (bin_image[row][67] == 255 && bin_image[row - 1][67] == 0 && bin_image[row - 2][67] == 0) {
//            single_line[1][0] = row;
//            single_line[1][1] = 67;
//            single_right = 2;  //δ����
//            break;
//        }
//        if (row == mid_row - 5){
//            single_right = 1;  //����
//        }
//    }
//    for (int row = image_h - 1; row >= mid_row - 5; row--){
//        if (bin_image[row][27] == 255 && bin_image[row - 1][27] == 0 && bin_image[row - 2][27] == 0){
//            single_line[0][0] = row;
//            single_line[0][1] = 27;
//            single_left = 2;    //δ����
//            break;
//        }
//        if (row == mid_row - 5) {
//            single_left = 1;  //����
//        }
//    }
//}
void image_process(void){
    turn_to_bin();  //ͼ���ֵ��
    compressimage(); //ͼ��ѹ��
    image_filter();
    image_draw_rectan(bin_image);
    get_start_point(start_line);
    search();
    get_left(L_search_count);
    get_right(R_search_count);
    get_turning_point();
    element();
}
void single_process(void) {
    turn_to_bin();  //ͼ���ֵ��
    compressimage(); //ͼ��ѹ��
    image_filter();
    image_draw_rectan(bin_image);
    get_single_point(5);
    single_search();
    single_element();
}
void process(void) {
    image_process();
    //corss_dir();
    //single_process();
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
    for (int i = 0; i < image_h; i++)
    {
        center_line[i][0] = i;
        center_line[i][1] = (left_line[i][1] + right_line[i][1]) >> 1;//������
    } 
    //draw_point(single_line[0][1],single_line[0][0], GREEN, Binimage_color)
    //draw_point(single_line[1][1], single_line[1][0], GREEN, Binimage_color)
    //draw_point(single_line[2][1], single_line[2][0], GREEN, Binimage_color)
    for (int i = 0; i < image_h; i++) {
        draw_point(left_line[i][1], left_line[i][0],GREEN, Binimage_color); 
        draw_point(right_line[i][1], right_line[i][0], GREEN, Binimage_color);
        draw_point(center_line[i][1], center_line[i][0], BLUE, Binimage_color);
    }
    //for (int i = 0; i < L_search_count; i++) {
    //    draw_point(left_points[i][1], left_points[i][0], PINK, Binimage_color);
    //    draw_point(right_points[i][1], right_points[i][0], PINK, Binimage_color);
    //    //draw_point(center_line[i][1], center_line[i][0], BLUE, Binimage_color);
    //}
    draw_point(start_left[1], start_left[0], BLUE, Binimage_color);
    //draw_point(right_points[single_amount - 1][1], right_points[single_amount - 1][0], BLUE, Binimage_color);
    //draw_point(left_points[single_amount - 1][1], left_points[single_amount - 1][0], BLUE, Binimage_color);

    draw_point(min_right_col, min_right_row, BLUE, Binimage_color);
    draw_point(max_left_col, max_left_row, BLUE, Binimage_color);
    if (circle_right_row)
    {
        draw_point(right_line[circle_right_row][1], right_line[circle_right_row][0], RED, Binimage_color);
    }
    if (circle_left_row) {
        draw_point(left_line[circle_left_row][1], left_line[circle_left_row][0], RED, Binimage_color);
    }   
    if (cross_right_row)
    {
        draw_point(right_line[cross_right_row][1], right_line[cross_right_row][0], BLUE, Binimage_color);
    }
    if (cross_left_row)
    {
        draw_point(left_line[cross_left_row][1], left_line[cross_left_row][0], BLUE, Binimage_color);
    }
    // ��ʾͼ��
    write_to_csv(bin_image);
    namedWindow("��ֵ��", WINDOW_NORMAL);
    imshow("��ֵ��", Binimage_color);
    resizeWindow("��ֵ��", 360, 240);
}
int picture(uint8_t type, const char* filepath) {
    if (type == 1) {
        //ͼƬ����
        const string ROAD = filepath;
        Mat roadImage = imread(ROAD);
        if (roadImage.empty())
        {
            cerr << "Failed to load image!" << endl;
            return -1;
        }
        dataTrans(roadImage, image);
        process();

        namedWindow("output", WINDOW_NORMAL);
        imshow("output", roadImage);
        resizeWindow("output", 360, 240);
    }
    else ;
}
int video(uint8_t type, const char* filepath) {
    if (type == 1) {
        //��Ƶ����
        VideoCapture cap(filepath);
        // �����Ƶ�Ƿ�ɹ���
        if (!cap.isOpened()) {
            std::cerr << "Error opening video file" << std::endl;
            return -1;
        }
        Mat frame;
        int frame_number = 0;
        // ѭ����ȡÿһ֡
        while (cap.read(frame)) {
            // �������ÿһ֡���д���
            dataTrans(frame, image);
            process();
            if (zqzq == 1) {
                string filename = "zqzq" + to_string(frame_number) + ".BMP";
                imwrite(filename, frame);
            }
            //if (circle_right_row) {
            //    string filename = "zcircle_right_row" + to_string(frame_number) + ".BMP";
            //    imwrite(filename, frame);
            //}
            //if (circle_left_row) {
            //    string filename = "zcircle_left_row" + to_string(frame_number) + ".BMP";
            //    imwrite(filename, frame);
            //}
            namedWindow("output", WINDOW_NORMAL);
            imshow("output", frame);
            resizeWindow("output", 360, 240);
            //����ÿһ֡ΪͼƬ

            string filename = "zframe_" + to_string(frame_number) + ".BMP";
            imwrite(filename, frame);

            if (waitKey(1) == 27) {
                break;
            }
            frame_number++;
        }
    }
    else;
}
int main() {
    const char* picture_filepath = "/home/anysets/Documents/AIMC/frames/frame_15.png";
    const char* video_filepath = "/home/anysets/Documents/AIMC/frames/frame_15.png";
    //video(1, video_filepath);
    picture(1, picture_filepath);

    waitKey(0);
    destroyAllWindows();
    return 0;
}
