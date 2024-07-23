#include "image_process.h"
#include <fstream>
#include <vector>
// #include <stdlib.h>
#define _CRT_SECURE_NO_WARNINGS
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
#define ROAD  "E:\\Files\\OpenCV++\\test\\test\\163.png"
// #define VIDEO  "C:\\Users\\wqg20\\Downloads\\23.mp4"
#define VIDEO "/home/anysets/Documents/Video/rc.mp4"
#define PVFLAG 1
// #define ROAD "E:\\Files\\OpenCV++\\test.BMP"
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 94
//注意这个border_length不能太多，会导致数组越界
#define border_length 60

uint8_t image_center = image_w / 2;

//将图像数组打印至csv文件--------------------------------------------------------
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
#define GrayScale 256

#define BLOCK_SIZE 15
#define C 5

void adaptive_threshold(unsigned char image[][IMAGE_W], int width, int height, int block_size, int c) {
    int half_block = block_size / 2;
    int sum, count, x, y, i, j;

    // 创建一个副本来保存原始图像数据
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

            // 计算局部均值
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

            // 计算局部阈值
            int mean = sum / count;
            int threshold = mean - c;

            // 应用阈值
            if (original[y][x] > threshold) {
                image[y][x] = 255;
            }
            else {
                image[y][x] = 0;
            }
        }
    }
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
void turn_to_bin(void) {
    uint8_t image_threshold = otsuThreshold(image);
    // uint8_t image_threshold = 176;
    // printf("OTSU=%d\n", image_threshold);

    // 遍历原始图像数组，根据动态阈值进行二值化
    for (int i = 0; i < IMAGE_H; ++i) {
        for (int j = 0; j < IMAGE_W; ++j) {
            if (image[i][j] > image_threshold) {
                image[i][j] = 255; // 如果原始图像灰度值大于阈值，设置为白色像素
            }
            else {
                image[i][j] = 0; // 如果原始图像灰度值小于等于阈值，设置为黑色像素
            }
        }
    }
}

//图像压缩一倍
void image_compress(uint8_t image[IMAGE_H][IMAGE_W], uint8_t iamge_zip[image_h][image_w]) {
	for (uint8_t row = 0; row < IMAGE_H; row+=2) {
		for (uint8_t col = 0; col < IMAGE_W; col+=2) {
			iamge_zip[row/2][col/2] = image[row][col];
		}
	}
}
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
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //把图像最上面画黑框
        image[0][i] = 0;
        image[1][i] = 0;
        image[image_h - 1][i] = 0;  //在图像最下面画一行黑框，以解决数组越界问题
    }
}

//画个黑框
void image_draw_rectan(uint8_t image[image_h][image_w]) {
    for (uint8_t i = 0; i < image_h; i++) { //把前两列和最后两列画黑框
        image[i][0] = 0;
        image[i][1] = 0;
        image[i][image_w - 1] = 0;
        image[i][image_w - 2] = 0;
    }
    for (uint8_t i = 0; i < image_w; i++) { //把图像最上面画黑框
        image[0][i] = 0;
        image[1][i] = 0;
        image[image_h - 1][i] = 0;  //在图像最下面画一行黑框，试图解决数组越界问题
    }
}

//找到起点，从最底部的中心开始寻找
int32_t start_point_l[2] = { 0 };
int32_t start_point_r[2] = { 0 };
int32_t start_line = image_h - 2;  //最后一行被画黑了，上移一行
int8_t getStartPoint()
{
	uint8_t left_flag = 0;
	uint8_t right_flag = 0;
	//初始化
	start_point_l[0] = 0;
	start_point_l[1] = 0;
	start_point_r[0] = 0;
	start_point_r[1] = 0;
	//从中间向左寻找
	for (uint32_t i = image_center; i > 1; i--)  //i>1的原因是0和1两列已经被全部置0了，且若i=0或1时下面数组会越界
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i - 1] == 0 && bin_image[start_line][i - 2] == 0)
		{
			start_point_l[0] = start_line; //行
			start_point_l[1] = i-1; //列
			left_flag = 1;
			break;
		}
	}
	//从中间向右寻找
	for (uint32_t i = image_center; i < image_w - 2; i++)  //i < image_w-2的原因同上
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i + 1] == 0 && bin_image[start_line][i + 2] == 0)
		{
			start_point_r[0] = start_line; //x坐标
			start_point_r[1] = i+1;
			right_flag = 1;
			break;
		}
	}
	if (left_flag && right_flag) { return 1; }
	else { return 1; }
}

//八邻域rewrite
//记录左边界的数据
int32_t border_location_left[border_length][2];  //记录边界位置
int32_t border_count_left = 0;
int32_t growth_direction_left[border_length];  //记录生长方向
int32_t growth_count_left = 0;
//记录右边界的数据
int32_t border_location_right[border_length][2];  //记录边界位置
int32_t border_count_right = 0;
int32_t growth_direction_right[border_length];  //记录生长方向
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
    //把起点加进去
    border_location_left[0][0] = start_point_l[0];
    border_location_left[0][1] = start_point_l[1];
    border_location_right[0][0] = start_point_r[0];
    border_location_right[0][1] = start_point_r[1];
    for (int j = 0; j < border_length-1; j++)
    {
		//左边界部分
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
        //右边界部分
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
        //当两线交汇时停止搜线
        if ((j != 0 && j != border_length-1) && 
        ((border_location_left[j][0] == border_location_right[j][0] && border_location_left[j][1] == border_location_right[j][1])
        || (border_location_left[j+1][0] == border_location_right[j-1][0] && border_location_left[j+1][1] == border_location_right[j-1][1])
        || (border_location_left[j-1][0] == border_location_right[j+1][0] && border_location_left[j-1][1] == border_location_right[j+1][1])))
        {
            break;
        }
    }
	// 用于分析生长方向
	// printf("growth_count_left: %d", growth_count_left);
	// for (int i = 0; i< 200; i++)
	// {
	// 	printf("%d\n", growth_direction_left[i]);
	// }

    // 检查border_location_left
    // for (int i = 0; i < border_length; i++)
    // {
    //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
    // }
}

// 存储优化后的左边界，每行像素都有份，注意此处的取值范围为-127～127，正好在压缩后的图像范围内
// 这里一定要注意，此处是二维数组，但是第0个值是image_w，第1个值是当前点在原边界的索引，当前值的索引才是image_h
uint8_t border_location_left_optimized[image_h][2] = { 0 };
uint8_t border_location_right_optimized[image_h][2] = { 0 };
void optimizeBorder()
{
    

    // ------------------------------------------------------------左边界的部分
    // 循环已经提取到的边界，去除不需要的点
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //一行只要一个就行了，行的索引是[0]，左边界需要的是每行最右侧的点
        if (border_location_left[i][1] > border_location_left_optimized[border_location_left[i][0]][0])
        {
            border_location_left_optimized[border_location_left[i][0]][0] = border_location_left[i][1];
            border_location_left_optimized[border_location_left[i][0]][1] = i;
        }  
    }
    // for (int i = 0; i < image_h; i++)
    // {
    //     // bin_image[i][border_location_left_optimized[i][0]] = 152;  // 将边界标记出来
    //     // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
    // }
    // 右边界的部分
    // 循环已经提取到的边界，去除不需要的点
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //一行只要一个就行了，行的索引是[0]，右边界需要的是每行最左侧的点
        if (border_location_right_optimized[border_location_right[i][0]][0] == 0 || border_location_right[i][1] < border_location_right_optimized[border_location_right[i][0]][0])
        {
            border_location_right_optimized[border_location_right[i][0]][0] = border_location_right[i][1];
            border_location_right_optimized[border_location_right[i][0]][1] = i;
        }
        
    }
    // for (int i = 0; i < image_h; i++)
    // {
    //     // bin_image[i][border_location_right_optimized[i][0]] = 152;  // 将边界标记出来
    //     // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
    // }   
}

// 获取中线
int32_t middle_line[border_length] = { 0 }; 
void getMiddleLine()
{
    memset(middle_line, 0, sizeof(middle_line));
    for (int i = 0; i < image_h; i++)
    {
        middle_line[i] = (border_location_left_optimized[i][0] + border_location_right_optimized[i][0]) / 2;
    }
}

// 获取角度
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

// 获取偏差值
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




//寻找拐点
uint8_t diffPoint_left[2] = {0}; // 记录拐点
uint8_t diffPoint_right[2] = {0}; // 记录拐点
uint8_t findTurningPoint()
{
    memset(diffPoint_left, 0, sizeof(diffPoint_left));
    memset(diffPoint_right, 0, sizeof(diffPoint_right));
    uint8_t flag = 0;
    uint8_t left_flag = 0;
    uint8_t right_flag = 0;
    uint8_t diffPoint_temp[2] = {};
    // ----------------------------------------------------------------------------------先找左下拐点
    // 1、首先使用优化过的边线进行判断。具体判断方法为如果遇到数组里的数字突变则判断为拐点
    // 从下往上扫描
    int diff = 0;  // 记录差值
    
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
    // 2、根据生长方向检查找到的拐点是否正确
    // 检查拐点附近的点的生长方向, 检查前后的两个点
    // border_location_right_optimized[diffPoint_temp[0]][1]代指拐点在原边界中的位置，以找到边界
    if (border_location_left_optimized[diffPoint_temp[0]][1] > 1 && border_location_left_optimized[diffPoint_temp[0]][1] < border_length-2)  // 防止数组越界
    {
        if ((growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // 拐点后第二个点
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 2 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 3 ||  // 拐点后第一格点
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 2 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 0 ||  // 拐点前第一格点
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 7 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1] == 1) &&
            (growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 0 ||  // 拐点前第二格点
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 7 ||
            growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2] == 1))
            {
                // printf("+2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+2]);
                // printf("+1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]+1]);
                // printf("-1 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-1]);
                // printf("-2 = %d\n",growth_direction_left[border_location_left_optimized[diffPoint_temp[0]][1]-2]);
                left_flag = 1;
            }  //无动作
        else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //不满足条件将点清零
        }
    }
    else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //不满足条件将点清零
        }
    // printf("left diff point 2: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    growth_direction_left[border_location_right_optimized[diffPoint_temp[0]][1]];

    // 将结果赋值
    diffPoint_left[0] = diffPoint_temp[0];
    diffPoint_left[1] = diffPoint_temp[1];

    // ----------------------------------------------------------------------------------再找右下拐点
    // 1、首先使用优化过的边线进行判断。具体判断方法为如果遇到数组里的数字突变则判断为拐点
    // 从下往上扫描
    diff = 0;  // 记录差值
    
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
    // 2、根据生长方向检查找到的拐点是否正确
    // 检查拐点附近的点的生长方向, 检查前后的两个点
    // border_location_right_optimized[diffPoint_temp[0]][1]代指拐点在原边界中的位置，以找到边界
    if (border_location_right_optimized[diffPoint_temp[0]][1] > 1 && border_location_right_optimized[diffPoint_temp[0]][1] < border_length-2)  // 防止数组越界
    {
        if ((growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 3 ||  // 拐点后第二个点
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 2 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 3 ||  // 拐点后第一格点
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 2 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 0 ||  // 拐点前第一格点
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 7 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1] == 1) &&
            (growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 0 ||  // 拐点前第二格点
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 7 ||
            growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2] == 1))
            {
                // printf("+2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+2]);
                // printf("+1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]+1]);
                // printf("-1 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-1]);
                // printf("-2 = %d\n",growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]-2]);
                right_flag = 2;
            }  //无动作
        else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //不满足条件将点清零
        }
    }
    else
        {
            diffPoint_temp[0] = 0;
            diffPoint_temp[1] = 0;  //不满足条件将点清零
        }
    // printf("right diff point 2: %d, %d\n",  diffPoint_temp[0], diffPoint_temp[1]);
    // growth_direction_right[border_location_right_optimized[diffPoint_temp[0]][1]];

    // 将结果赋值
    diffPoint_right[0] = diffPoint_temp[0];
    diffPoint_right[1] = diffPoint_temp[1];

    // 传出flag
    // 3、如果左右拐点高度相差过大，取相对下方的点作为拐点输出
    if (diffPoint_right[0] - diffPoint_left[0] > 25)
    {
        left_flag = 0;
    }
    else if (diffPoint_left[0] - diffPoint_right[0] > 25)
    {
        right_flag = 0;
    }

    // 4、此处为防止误判而写，若影响则考虑删除
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



// ---------------------------------------------------------------------------------------画线部分
// 根据两点求斜率
void two_slope(double x1, double y1, double x2, double y2, double* slope, double* intercept) {
    *slope = (y2 - y1) / (x2 - x1);
    *intercept = y1 - (*slope * x1);
}
// 这里一定要注意，此处是二维数组，但是第0个值是image_w，第1个值是当前点在原边界的索引，当前值的索引才是image_h
// 根据两点补线
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

// 记录当前状态的标志
// bool isInRightCircle = false;
// bool isInCrossing = false;
uint8_t inRightCircleFlag = 0;
uint8_t inCrossingFlag = 0;

uint8_t status3_bwline[2] = {0};  // black white line
long frame1 = 0;
uint32_t amount_fps[2] = {0, 0};  // 用于记录状态0时连续找到5次右拐点

uint8_t right_point_row1[2] = {0};  // 这个记录的是常态下的起点
// 找出右环岛的最左点
uint8_t LeftPoint_RightCircle[2];
void FindRightCircleLeftPoint(uint8_t flag)
{
    memset(LeftPoint_RightCircle, 0, sizeof(LeftPoint_RightCircle));
    if (flag == 0)
    {
        diffPoint_right[0] = right_point_row1[0];
        diffPoint_right[1] = right_point_row1[1];
    }
    // 从当前右拐点向上找圆
    // 当前右拐点：diffPoint_right
    for (int i = diffPoint_right[0]-1; i > 1; i--)
    {
        if (inRightCircleFlag == 1)
        {
            if (i < (image_h/4)*3 && bin_image[i][diffPoint_right[1]] == 0)  // 这里加了前面一个判断条件，状态1时这个点一定在图像靠上的部分
            {
                LeftPoint_RightCircle[0] = i;
                LeftPoint_RightCircle[1] = diffPoint_right[1];
                // printf("Found Circle\n");
                break;  // 找到了记录点并跳出循环
            }
        }
        else
        {
            if (bin_image[i][diffPoint_right[1]] == 0) 
            {
                LeftPoint_RightCircle[0] = i;
                LeftPoint_RightCircle[1] = diffPoint_right[1];
                // printf("Found Circle\n");
                break;  // 找到了记录点并跳出循环
            }
        }
        
    }
    // 用找到的八邻域找最左边的点
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
                    break;  // 如果原来的点在八邻域点的右边，则把该八邻域点记录下来
                }
                else
                {
                    return;
                }
            }
        }
    }
}

// 判断右圆环
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
                inRightCircleFlag = 1;  // 进入状态1
                getStartPoint();  //《屎山是怎么堆成的》
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
                    inRightCircleFlag = 2;  // 进入状态2
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
        getStartPoint();  //《屎山是怎么堆成的》
        if (start_point_r[1] < 86)
        {
            inRightCircleFlag = 3;
        }
    }
    else if (inRightCircleFlag == 3)
    {
        // 如果左起点在10列以前，直接作为起始点（此时右侧已经丢线）
        getStartPoint();  //《屎山是怎么堆成的》
        // printf("start_point_l[1]: %d\n", start_point_l[1]);
        uint8_t start_point = image_h - 2;  // 记录起始点
        uint8_t diffPoint_rightUp[2] = {0, 90};
        // && start_point_r[1] > 90
        if (start_point_l[1] < 10)
        // if (false)
        {
            printf("In if\n");
            // 使用优化过的左边界向右寻找拐点
            bool found_flag = false;  // 指示是否找到拐点，跳出嵌套循环
            
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
            // 找到拐点后开始画线
            // printf("%d\n", start_point);
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], image_w / 10, start_point_l[0]);
        }
        else
        {
            // 首先使用右边界寻找起始点
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
            // 使用优化过的左边界向右寻找拐点
            bool found_flag = false;  // 指示是否找到拐点，跳出嵌套循环
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
            // 找到拐点后开始画线
            // printf("%d\n", start_point);
            FillLine(border_location_left_optimized, diffPoint_rightUp[1], diffPoint_rightUp[0], start_point_l[1], start_point_l[0]);
        }
        // 找到拐点后重新爬八邻域
        // 先爬左边界，直到碰到最右侧的点
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
        // 检查border_location_left
        // for (int i = 0; i < border_length; i++)
        // {
        //     printf("border_location_left[%d]: %d, %d\n", i, border_location_left[i][0], border_location_left[i][1]);
        // }
        optimizeBorder();
        

        // 进入下一个状态
        // // 方案1:右上拐点消失
        // printf("diffPoint_rightUp[0]: %d\n", diffPoint_rightUp[0]);
        // printf("diffPoint_rightUp[1]: %d\n", diffPoint_rightUp[1]);
        // if (diffPoint_rightUp[0] > 50 && diffPoint_rightUp[1] < image_w/2)  // 左上拐点在45行以上
        // {
        //     inRightCircleFlag = 4;
        //     // 记录当前帧数
        //     amount_fps[1] = frame1;
        //     printf("Status 4\n");
        // }
        
        // 方案2：左起点先消失再出现  小型状态机
        // 
        getStartPoint();  //《屎山是怎么堆成的》
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
                inRightCircleFlag = 4;  // 进入状态4
            }
        }
    }
    else if (inRightCircleFlag == 4)
    {
        uint8_t right_point = 0;
        uint8_t diff = 0;
        // 没有达到指定帧数不允许进入下一状态
        if (frame1 - amount_fps[1] > 10)
        // if (true)
        {
            // 方案1
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

            // 方案2
            // // 用拐点判断
            // if (TurningPointFlag == 1 || TurningPointFlag == 3)
            // {
            //     FillLine(border_location_left_optimized, diffPoint_left[1], diffPoint_left[0], image_w/2, 0);
            //     inRightCircleFlag = 5;
            //     // break;
            // }

            // 方案3
            // 用左边界最右点判断
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

            // 方案4  学长的code
            for (int i = 5; i < 100; i++) 
            {
                if (border_location_left[i][1] >= border_location_left[i - 1][1] && border_location_left[i][1] > border_location_left[i + 1][1]
                    && border_location_left[i][1] > border_location_left[i - 10][1]) {   //行数最高 一定要有等于号
                    uint8_t row = border_location_left[i][0];
                    uint8_t col = border_location_left[i][1];
                    // printf("status4:  row: %d, col: %d\n", row, col);
                    //第一个点为左下角固定点 第二个点为右上固定
                    FillLine(border_location_left_optimized, image_w - 5, 5, col, row);  //压弯
                    break;
                }
            }

            // 如果左起点丢了，进入下一状态
            getStartPoint();  //《屎山是怎么堆成的》
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

        // 写死，斜着画线
        // 注意不能是从最左下角开始写死，得往中间靠，不然出不去
        FillLine(border_location_left_optimized, image_w - 1, image_h / 3, image_w / 3, image_h - 1);

        // 出环条件
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
        // 右下角有起点，退出环岛
        getStartPoint();  //《屎山是怎么堆成的》
        if (start_point_r[1] < 90)
        {
            inRightCircleFlag = 0;
            // printf("start_point_r[1]: %d\n", start_point_r[1]);
        }
    }
}

// 识别十字和圆环的重要的区别是十字如果只识别到一个拐点时，另一侧是丢线的
// 而圆环入环条件就是只识别到一个拐点，但另一侧是不丢线的
uint8_t up_right_point[2] = {0};
uint8_t up_left_point[2] = {0};
void findUpTurningPoint(uint8_t start_l, uint8_t start_r)
{
    
    // 找左上拐点
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
    // 找右上拐点
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

// 识别十字
void crossingHandle(uint8_t TurningPointFlag)
{
    if (inRightCircleFlag == 0 && (TurningPointFlag == 1 || TurningPointFlag == 3) && inCrossingFlag == 0)
    {
        getStartPoint();  //《屎山是怎么堆成的》
        if (start_point_r[1] > 90 && TurningPointFlag == 1)  // 只识别到左拐点时，检查是否丢右起点，若丢了，进入状态1
        {
            inCrossingFlag = 1;
        }
        else if(TurningPointFlag == 3)  // 两个拐点都识别到，进入状态1
        {
            inCrossingFlag = 1;
        }
    }
    else if (inCrossingFlag == 1)
    {
        // 即使进入了十字，也有可能只识别到了一个拐点，所以要继续判断
        // 如果两个拐点都是能找到的，向上找到上方的两个拐点
        if (TurningPointFlag == 3)
        {
            findUpTurningPoint(diffPoint_left[0], diffPoint_right[0]);

            // 画线
            FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], diffPoint_left[1], diffPoint_left[0]);
            FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], diffPoint_right[1], diffPoint_right[0]);
        }

        // 两侧均丢起点，结束状态1
        getStartPoint();  //《屎山是怎么堆成的》
        if (start_point_l[1] < 4 && start_point_r[1] > 90)
        {
            inCrossingFlag = 2;
            memset(up_left_point, 0, sizeof(up_left_point));
            memset(up_right_point, 0, sizeof(up_right_point));
        }
    }
    else if (inCrossingFlag == 2)
    {
        getStartPoint();  //《屎山是怎么堆成的》
        // printf("start_point_l[0]: %d\n", start_point_l[0]);
        // printf("start_point_r[0]: %d\n", start_point_r[0]);
        // printf("up_left_point[1]: %d\n", up_left_point[1]);
        // printf("up_left_point[0]: %d\n", up_left_point[0]);
        // printf("up_right_point[1]: %d\n", up_right_point[1]);
        // printf("up_right_point[0]: %d\n", up_right_point[0]);
        findUpTurningPoint(start_point_l[0]-1, start_point_r[0]-1);
        // 画线  注意这里的（25, 58）（75, 58）是普适的左右起点 
        FillLine(border_location_left_optimized, up_left_point[1], up_left_point[0], 25, 58);
        FillLine(border_location_right_optimized, up_right_point[1], up_right_point[0], 75, 58);
        // 如果两侧不丢线了，退出十字路口
        
        if (start_point_l[1] > 4 && start_point_r[1] < 90)
        {
            inCrossingFlag = 0;
        }
    }
}

// 元素判断主函数
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
    // ----------常态化图像处理----------
    turn_to_bin();  //图像二值化
    // adaptive_threshold(image, IMAGE_W, IMAGE_H, BLOCK_SIZE, C);  // 自适应二值化
    image_compress(image,bin_image); //图像压缩
    // image_filter();
    image_draw_rectan(bin_image);
    // 获得起点
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

        //用opencv显示
        Mat Binimage(image_h, image_w, CV_8UC1);
        for (int i = 0; i < image_h; i++) {
            for (int j = 0; j < image_w; j++) {
                Binimage.at<uint8_t>(i, j) = bin_image[i][j];
            }
        }

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
            //draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
            draw_point(border_location_right_optimized[i][0], i, GREEN, Binimage_color);
            // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
            // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
            draw_point(middle_line[i], i, BLUE, Binimage_color);
            draw_point(diffPoint_left[1], diffPoint_left[0], RED, Binimage_color);
            draw_point(diffPoint_right[1], diffPoint_right[0], RED, Binimage_color);
            // draw_point(LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], RED, Binimage_color);
        }

        // 显示图像
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
        //视频处理
        VideoCapture cap(VIDEO);
        // 检查视频是否成功打开
        if (!cap.isOpened()) {
            std::cerr << "Error opening video file" << std::endl;
            return -1;
        }
        Mat frame;
        // 循环读取每一帧
        while (cap.read(frame)) {
            printf("------------------------------------------frame:%d------------------------------------------\n", frame1);
            //保存每一帧为图片
            string filename = to_string(frame1) + ".png";
            imwrite(filename, frame);
            // waitKey(1000);
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
                draw_point(border_location_left_optimized[i][0], i, GREEN, Binimage_color);
                draw_point(border_location_right_optimized[i][0], i, GREEN, Binimage_color);
                // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
                //  draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
                draw_point(middle_line[i], i, BLUE, Binimage_color);
                draw_point(diffPoint_left[1], diffPoint_left[0], RED, Binimage_color);
                draw_point(diffPoint_right[1], diffPoint_right[0], RED, Binimage_color);
                // draw_point(LeftPoint_RightCircle[1], LeftPoint_RightCircle[0], RED, Binimage_color);
            }
            // 显示图像
            namedWindow("bin_image", WINDOW_NORMAL);
            imshow("bin_image", Binimage_color);
            resizeWindow("bin_image", 360, 240);
            // namedWindow("output", WINDOW_NORMAL);
            // imshow("output", frame);
            // resizeWindow("output", 360, 240);
            // 等待一段时间，按下ESC键退出
            if (waitKey(20) == 27) {
                break;
            }
            
        }
    }
    
    
	return 0;
}
