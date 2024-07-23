#define _CRT_SECURE_NO_WARNINGS
#include "image_process.h"

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
#define start_line image_h - 2
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
        //printf("Error opening file!\n");
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
//绝对值函数
double double_abs(double value)
{
    if (value >= 0) return value;
    else return -value;
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
uint32_t black_points;
void turn_to_bin(void) {
    black_points = 0;
    uint8_t image_threshold = otsuThreshold(image);
    // 遍历原始图像数组，根据动态阈值进行二值化
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // 如果原始图像灰度值大于阈值，设置为白色像素
            }
            else {
                image[i][j] = 0; // 如果原始图像灰度值小于等于阈值，设置为黑色像素
                black_points++;
            }
        }
    }
    //printf("黑色的点%d\n", black_points);
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
#define threshold_max	(255 * 8) // 此参数可根据自己的需求调节
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
        image[39][i] = 0;
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
//单边找起点
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
//八邻域爬线
#define search_amount 80   //允许单边最多搜索到的点
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
                //printf("第%d次，L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
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
                //printf("第%d次，R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
//单边八邻域
uint8_t single_amount = 20; //单边最多搜索的点
void single_search(void) {
    L_search_count = R_search_count = 0;    //搜到的点数清零
    left_points[0][0] = start_left[0];  //y row
    left_points[0][1] = start_left[1];  //X col
    int8_t L_curr_row = start_left[0];   //初始化row坐标
    int8_t L_curr_col = start_left[1];  //初始化col坐标
    for (int j = 1; j < single_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[L_curr_row + L_dir_row[i]][L_curr_col + L_dir_col[i]] == 0)
                && (bin_image[L_curr_row + L_dir_row[(i + 1) & 7]][L_curr_col + L_dir_col[(i + 1) & 7]] == 255)) {
                //由黑到白
                L_search_count++;
                L_curr_row += L_dir_row[(i + 1) & 7];
                L_curr_col += L_dir_col[(i + 1) & 7];
                left_points[L_search_count][0] = L_curr_row;   //存储白点
                left_points[L_search_count][1] = L_curr_col;   //存储白点 
                dir_left[L_search_count] = i;
                //printf("第%d次，L_curr_row = %d,L_curr_col = %d\n", i, L_curr_row, L_curr_col);
                break;
            }
        }
    }
    right_points[0][0] = start_right[0];  //y row
    right_points[0][1] = start_right[1];  //X col
    int8_t R_curr_row = start_right[0];   //初始化row坐标
    int8_t R_curr_col = start_right[1];  //初始化col坐标
    for (int j = 1; j < single_amount; j++) {
        for (int i = 0; i < 8; i++) {
            if ((bin_image[R_curr_row + R_dir_row[i]][R_curr_col + R_dir_col[i]] == 0)
                && (bin_image[R_curr_row + R_dir_row[(i + 1) & 7]][R_curr_col + R_dir_col[(i + 1) & 7]] == 255)) {
                R_search_count++;
                R_curr_row += R_dir_row[(i + 1) & 7];
                R_curr_col += R_dir_col[(i + 1) & 7];
                right_points[R_search_count][0] = R_curr_row;   //存储白点
                right_points[R_search_count][1] = R_curr_col;   //存储白点
                dir_right[R_search_count] = i;
                //printf("第%d次，R_curr_row = %d,R_curr_col = %d\n", dir_right[R_search_count], R_curr_row, R_curr_col);
                break;
            }
        }
    }
}
uint8_t single_left;    //1是丢线
uint8_t single_right;   //1是丢线
uint8_t single_center;  //1是丢线
void single_element() {
    single_left = single_right = single_center =0;
    uint8_t center_count = 0;
    for (center_count = image_h - 2; center_count > single_amount; center_count--) {
        if (bin_image[center_count][image_w / 2] == 255 && bin_image[center_count - 1][ image_w / 2] == 0) {    //白到黑
            break;
        }
    }
    if (center_count == single_amount) {
        single_center = 1;
        printf("中间单边丢线\n");
    }
    for (uint8_t count = 0; count < single_amount - 1; count++) {
        if (right_points[count][0] > right_points[count + 1][0]) {
            uint8_t row = right_points[count][0];
            uint8_t col = right_points[count][1];
            if(bin_image[row - 5][col - 5] == 255){
                right_points[single_amount - 1][0] = right_points[count][0];    //丢线
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
//提取边线
uint8_t left_line[image_h][2];  //左线数组 row col
uint8_t right_line[image_h][2];    //右线数组 row col
uint8_t center_line[image_h][2];   //中线数组 row col
uint8_t width[image_h]; //赛道宽度
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
//最小二乘法计算直线拟合的斜率和截距
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
//根据两点求斜率
void two_slope(double x1, double y1, double x2, double y2, double* slope, double* intercept) {
    *slope = (y2 - y1) / (x2 - x1);
    *intercept = y1 - (*slope * x1);
}
//根据两点补线
void two_line(uint8_t line[image_h][2], uint8_t x1, uint8_t y1, uint8_t x2, uint8_t y2) {
    double slope, intercept;
    two_slope(x1,y1,x2,y2,&slope,&intercept);
    for (int i = y1; i >= y2; i--) {
        line[i][1] = (line[i][0] - intercept) / slope;
    }
}
// 计算拟合的方差
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
//生长方向判断
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
//横向扫线
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
                        break;  //每一行只搜第一个点
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
    for (uint8_t row = start_row; row > 5; row--) { //开始网上扫找到第一个黑点 再精准定位
        if (bin_image[row][start_col] == 255 && bin_image[row - 1][start_col] == 0) {  //白黑 找到了
            for (uint8_t i = row - 1; i >= 5; i--) {
                for (uint8_t j = col_middle; j > 3; j--) {
                    if (bin_image[i][j] == 255 && bin_image[i - 1][j] == 0) {
                        if (max_left_col < j) {
                            max_left_col = j;
                            max_left_row = i;
                        }
                        break;  //每一行只搜第一个点
                    }
                }
            }
            if (max_left_col > 10) return 1;
            else return 0;
        }
    }
}
uint8_t element_start_line = 38; //从下面计数
uint8_t element_end_line = 10;
uint8_t circle_right_row;
uint8_t circle_left_row;
uint8_t cross_right_row;
uint8_t cross_left_row;
void get_turning_point(void) {
    //标志位清零
    circle_right_row = circle_left_row = cross_right_row = cross_left_row =  0;
    for (int i = 0; i < R_search_count; i++) {  //右拐点判断
        double slope, intercept;
        uint8_t turn_row = right_points[i][0];
        if (turn_row <= element_end_line) break;
        if (turn_row >= element_start_line) continue;
        //连续丢线进入元素判断
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
        //连续丢线进入元素判断
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
//判断是左十字还是右边十字，右边为1，左边为2
uint8_t corss_dir(void) {
    turn_to_bin();  //图像二值化
    compressimage(); //图像压缩
    image_filter();
    image_draw_rectan(bin_image);
    get_start_point(10);
    search();
    uint32_t sum_center_line = 0;
    uint8_t startline = 9;
    uint8_t endline = 5;
    for (int i = startline; i > endline; i--)
    {
        center_line[i][1] = (left_points [i] [1] + right_points[i][1]) >> 1;//求中线
        sum_center_line += center_line[i][1];
    }        
    float temp_value = (float)sum_center_line / (startline - endline);
    printf("temp_value = %f", temp_value);
    if (temp_value - image_w / 2 > 0) {
        printf("右边十字\n");
        return 1;
    }
    else {
        printf("左边十字\n");
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
        printf("右圆环");

        //右圆环标志位
        circle_right_last_flag = 0;
    }
    if (circle_left_last_flag == 1 && circle_left_flag == 0) {
        printf("左圆环");       
        zqzq = 1;
        //左圆环标志位
        circle_left_last_flag = 0;
    }
    if (cross_last_flag == 1 && cross_flag == 0) {
        printf("十字");
        //十字标志位
        cross_last_flag = 0;
    }
}

//uint8_t single_line[3][2]; //row col [0]代表左边 [1]代表右边 [2]代表中间
//uint8_t single_left = 0, single_right = 0,single_middle = 0; //1代表丢线 2代表未丢线
//uint8_t single_end_line = image_h - 10 ;  //从下方最大搜多少行
//void single_edge(){
//    int mid_row;
//    single_left = 0, single_right = 0;  //清零
//    for (int row = image_h - 1; row >= 1; row--) {
//        if (bin_image[row][47] == 255 && bin_image[row - 1][47] == 0 && bin_image[row - 2][47] == 0){  //白黑黑
//            single_line[2][0] = row;
//            single_line[2][1] = 47;
//            mid_row = row;  //图像中间的行
//            break;
//        }
//        if (row <= single_end_line){
//            single_middle = 1;  //丢线
//            mid_row = row;
//        }
//    }
//    for (int row = image_h - 1; row >= mid_row - 5; row--){ //最多搜到中间行5行以上
//        if (bin_image[row][67] == 255 && bin_image[row - 1][67] == 0 && bin_image[row - 2][67] == 0) {
//            single_line[1][0] = row;
//            single_line[1][1] = 67;
//            single_right = 2;  //未丢线
//            break;
//        }
//        if (row == mid_row - 5){
//            single_right = 1;  //丢线
//        }
//    }
//    for (int row = image_h - 1; row >= mid_row - 5; row--){
//        if (bin_image[row][27] == 255 && bin_image[row - 1][27] == 0 && bin_image[row - 2][27] == 0){
//            single_line[0][0] = row;
//            single_line[0][1] = 27;
//            single_left = 2;    //未丢线
//            break;
//        }
//        if (row == mid_row - 5) {
//            single_left = 1;  //丢线
//        }
//    }
//}
void image_process(void){
    turn_to_bin();  //图像二值化
    compressimage(); //图像压缩
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
    turn_to_bin();  //图像二值化
    compressimage(); //图像压缩
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
    for (int i = 0; i < image_h; i++)
    {
        center_line[i][0] = i;
        center_line[i][1] = (left_line[i][1] + right_line[i][1]) >> 1;//求中线
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
    // 显示图像
    write_to_csv(bin_image);
    namedWindow("二值化", WINDOW_NORMAL);
    imshow("二值化", Binimage_color);
    resizeWindow("二值化", 360, 240);
}
int picture(uint8_t type, const char* filepath) {
    if (type == 1) {
        //图片处理
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
        //视频处理
        VideoCapture cap(filepath);
        // 检查视频是否成功打开
        if (!cap.isOpened()) {
            std::cerr << "Error opening video file" << std::endl;
            return -1;
        }
        Mat frame;
        int frame_number = 0;
        // 循环读取每一帧
        while (cap.read(frame)) {
            // 在这里对每一帧进行处理
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
            //保存每一帧为图片

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
