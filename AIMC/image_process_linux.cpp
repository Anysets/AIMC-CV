#include "image_process.h"
#include <fstream>
#include <vector>
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
#define ROAD  "/home/anysets/Documents/AIMC-CV/AIMC/frame_40.png"
#define VIDEO  "/home/anysets/Documents/Video/3.mp4"
// #define ROAD "E:\\Files\\OpenCV++\\test.BMP"
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 94
//注意这个border_length不能太多，会导致数组越界
#define border_length 100
uint8_t image_center = image_w / 2;

//将图像数组打印至csv文件--------------------------------------------------------
void printToCSV(uint8_t matrix[][image_w], const std::string& filename) {
	std::ofstream outputFile(filename);

	if (!outputFile.is_open()) {
		std::cerr << "无法打开文件：" << filename << std::endl;
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
	printf("print to csv successed\n");
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
    printf("OTSU=%d\n", image_threshold);

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
        image[image_h - 1][i] = 0;  //在图像最下面画一行黑框，以解决数组越界问题
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
	// cout << "Found Left Start Point: " << start_point_l[0] << ", " << start_point_l[1] << endl;
    // cout << "Found Right Start Point: " << start_point_r[0] << ", " << start_point_r[1] << endl;
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
	else { return 0; }
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
    border_count_left = 0;
    growth_count_left = 0;
    border_count_right = 0;
    growth_count_right = 0;
    for (int j = 0; j < border_length; j++)
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
                bin_image[start_point_l[0] + neighbor_left[i+1][0]][start_point_l[1] + neighbor_left[i+1][1]] = 151;
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
                // printf("New Right Start Point: %d, %d\n", start_point_r[0], start_point_r[1]);
				// printf("growth_direction_right: %d\n", growth_direction_right[j]);
                // printf("finished\n");
				// printf("%d\n", j);
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
uint8_t border_location_left_optimized[image_h] = { 0 };
uint8_t border_location_right_optimized[image_h] = { 0 };
void optimizeBorder()
{
    memset(border_location_left_optimized, 0, sizeof(border_location_left_optimized));
    memset(border_location_right_optimized, 0, sizeof(border_location_right_optimized));

    // 左边界的部分
    // 循环已经提取到的边界，去除不需要的点
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        // border_location_left[i][0];  //一行只要一个就行了，行的索引是[0]，左边界需要的是每行最右侧的点
        // 开始判断
        // printf("i=%d, border_location_left[i][1]=%d, border_location_left[i][0]=%d\n", i, border_location_left[i][1], border_location_left[i][0]);
        if (border_location_left[i][1] > border_location_left_optimized[border_location_left[i][0]])
        {
            border_location_left_optimized[border_location_left[i][0]] = border_location_left[i][1];
        }  
    }
    for (int i = 0; i < image_h; i++)
    {
        bin_image[i][border_location_left_optimized[i]] = 152;
        // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
    }
    // 右边界的部分
    // 循环已经提取到的边界，去除不需要的点
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //一行只要一个就行了，行的索引是[0]，右边界需要的是每行最左侧的点
        // 开始判断
        // 右边界需要有一步初始化
        // if (border_location_right_optimized[border_location_right[i][0]] == 0)
        // {
        //     border_location_right_optimized[border_location_right[i][0]] = border_location_right[i][1];
        // }
        if (border_location_right_optimized[border_location_right[i][0]] == 0 || border_location_right[i][1] < border_location_right_optimized[border_location_right[i][0]])
        {
            border_location_right_optimized[border_location_right[i][0]] = border_location_right[i][1];
        }
        
    }
    for (int i = 0; i < image_h; i++)
    {
        bin_image[i][border_location_right_optimized[i]] = 152;
        // printf("border_optimized:%d, %d\n", i, border_location_left_optimized[i]);
    }
}

// 获取中线
int32_t middle_line[border_length] = { 0 }; 
void getMiddleLine()
{
    memset(middle_line, 0, sizeof(middle_line));
    for (int i = 0; i < image_h; i++)
    {
        middle_line[i] = (border_location_left_optimized[i] + border_location_right_optimized[i]) / 2;
    }
}

// 获取偏差值
int32_t deviation_number;
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
}

void image_process(void)
{
    turn_to_bin();  //图像二值化
    image_compress(image,bin_image); //图像压缩
    image_filter();
    image_draw_rectan(bin_image);
    if (getStartPoint())  //先获得起点
	{
		printf("Start point found.\n");
        cout << "Found Left Start Point: " << start_point_l[0] << ", " << start_point_l[1] << endl;
        cout << "Found Right Start Point: " << start_point_r[0] << ", " << start_point_r[1] << endl;
        neighborSearch();
		// search((uint16)USE_num, image, bin_image, &data_statics_l, &data_statics_r, start_point_l[0], start_point_l[1], start_point_r[0], start_point_r[1], &hightest);
	    optimizeBorder();
        getMiddleLine();
        getDeviation();
    }
    else
    {
        printf("Find start point failed.\n");
    }

    printToCSV(bin_image, "output.csv");
}

int main(){
    /*
	// const string ROAD = "/home/anysets/Documents/Photos/test.BMP";
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
        draw_point(border_location_left_optimized[i], i, GREEN, Binimage_color);
        draw_point(border_location_right_optimized[i], i, GREEN, Binimage_color);
        // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
        // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
        draw_point(middle_line[i], i, RED, Binimage_color);
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
    */
    
    //视频处理
    VideoCapture cap(VIDEO);
    // 检查视频是否成功打开
    if (!cap.isOpened()) {
        std::cerr << "Error opening video file" << std::endl;
        return -1;
    }
    Mat frame;
    int frame_number = 0;
    // 循环读取每一帧
    while (cap.read(frame)) {
        printf("------------------------------------------frame:%d------------------------------------------\n", frame_number);
         //保存每一帧为图片
        string filename = "frame_" + to_string(frame_number) + ".png";
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
            draw_point(border_location_left_optimized[i], i, GREEN, Binimage_color);
            draw_point(border_location_right_optimized[i], i, GREEN, Binimage_color);
            // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
            // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
            draw_point(middle_line[i], i, RED, Binimage_color);
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
        frame_number++;
    }
    
	return 0;
}
