#include <iostream>
#include <fstream>
#include <vector>
#include <opencv2/opencv.hpp>
using namespace std;
using namespace cv;

//#define ROADBMP  "E:\\Files\\OpenCV++\\14.bmp"
#define ROADBMP  "/run/media/anysets/Files/Files/OpenCV++/test.BMP"

//数据类型
typedef   signed                char int8;
typedef   signed short      int int16;
typedef   signed                int int32;
typedef unsigned              char uint8;
typedef unsigned short     int uint16;
typedef unsigned               int uint32;

//自定义图像大小
#define image_h	120
#define image_w	188
uint32 image_center = image_w / 2;

//读取图片
Mat img = imread("/run/media/anysets/Files/Files/OpenCV++/test.BMP");

//图像数组声明
//为使该段程序能直接移植至单片机中
uint8 image[image_h][image_w];
uint8 bin_image[image_h][image_w];

//将三通道Mat转换为单通道数组，结果存放在image数组中
void Data_trans(Mat image_)
{
    int channels = image_.channels();
    cout << "channels的值为：" << channels << endl;
    //uint16 row = 0, line = 0;

    /*at()函数
    对于单通道图像"picture1"，picture1.at<uchar>(i,j)就表示在第i行第j列的像素值。
    对于多通道图像如RGB图像"picture2"，可以用picture2.at<Vec3b>(i,j)[c]来表示某个通道中在(i,j)位置的像素值。*/
    for (int i = 0; i < image_h; i += 1)
    {
        for (int j = 0; j < image_w; j += 1)
        {
            image[i][j] = image_.at<uchar>(i, j * channels);
        }
    }        
}

//滤波
#define threshold_max	255*5
#define threshold_min	255*2
void image_filter(void)
{
	uint16 i, j;
	uint32 num = 0;
	int channels = img.channels();
	for (i = 1; i < image_h - 1; i++)
	{
		for (j = 1; j < (image_w - 1); j++)
		{
			num =
				image[i - 1][j - 1] + image[i - 1][j] + image[i - 1][j + 1]
				+ image[i][j - 1] + image[i][j + 1]
				+ image[i + 1][j - 1] + image[i + 1][j] + image[i + 1][j + 1];
			if (num >= threshold_max && image[i][j] == 0)
			{
				//printf("黑滤波—image[%d][%d]\n", i, j);
				image[i][j] = 255;
				//三通道赋值
				img.at<Vec3b>(i, j)[0] = 255;
				img.at<Vec3b>(i, j)[1] = 255;
				img.at<Vec3b>(i, j)[2] = 255;
			}
			if (num <= threshold_min && image[i][j] == 255)
			{
				//printf("白点滤波—image[%d][%d]\n", i, j);
				image[i][j] = 0;
				//三通道赋值
				img.at<Vec3b>(i, j)[0] = 0;
				img.at<Vec3b>(i, j)[1] = 0;
				img.at<Vec3b>(i, j)[2] = 0;
			}

		}
	}
}

//otsu法确定阈值
uint8 otsuThreshold(uint8* image, uint16 col, uint16 row)
{
#define GrayScale 256
	uint16 Image_Width = col;
	uint16 Image_Height = row;
	int X; uint16 Y;
	uint8* data = image;
	int HistGram[GrayScale] = { 0 };

	uint32 Amount = 0;
	uint32 PixelBack = 0;
	uint32 PixelIntegralBack = 0;
	uint32 PixelIntegral = 0;
	int32 PixelIntegralFore = 0;
	int32 PixelFore = 0;
	double OmegaBack = 0, OmegaFore = 0, MicroBack = 0, MicroFore = 0, SigmaB = 0, Sigma = 0; // 类间方差;
	uint8 MinValue = 0, MaxValue = 0;
	uint8 Threshold = 0;
	for (Y = 0; Y < Image_Height; Y++) //Y<Image_Height改为Y =Image_Height；以便进行 行二值化
	{
		//Y=Image_Height;
		for (X = 0; X < Image_Width; X++)
		{
			HistGram[(int)data[Y * Image_Width + X]]++; //统计每个灰度值的个数信息
		}
	}

	for (MinValue = 0; MinValue < 256 && HistGram[MinValue] == 0; MinValue++);        //获取最小灰度的值
	for (MaxValue = 255; MaxValue > MinValue && HistGram[MinValue] == 0; MaxValue--); //获取最大灰度的值

	if (MaxValue == MinValue)
	{
		return MaxValue;          // 图像中只有一个颜色
	}
	if (MinValue + 1 == MaxValue)
	{
		return MinValue;      // 图像中只有二个颜色
	}

	for (Y = MinValue; Y <= MaxValue; Y++)
	{
		Amount += HistGram[Y];        //  像素总数
	}

	PixelIntegral = 0;
	for (Y = MinValue; Y <= MaxValue; Y++)
	{
		PixelIntegral += HistGram[Y] * Y;//灰度值总数
	}
	SigmaB = -1;
	for (Y = MinValue; Y < MaxValue; Y++)
	{
		PixelBack = PixelBack + HistGram[Y];    //前景像素点数
		PixelFore = Amount - PixelBack;         //背景像素点数
		OmegaBack = (double)PixelBack / Amount;//前景像素百分比
		OmegaFore = (double)PixelFore / Amount;//背景像素百分比
		PixelIntegralBack += HistGram[Y] * Y;  //前景灰度值
		PixelIntegralFore = PixelIntegral - PixelIntegralBack;//背景灰度值
		MicroBack = (double)PixelIntegralBack / PixelBack;//前景灰度百分比
		MicroFore = (double)PixelIntegralFore / PixelFore;//背景灰度百分比
		Sigma = OmegaBack * OmegaFore * (MicroBack - MicroFore) * (MicroBack - MicroFore);//g
		if (Sigma > SigmaB)//遍历最大的类间方差g
		{
			SigmaB = Sigma;
			Threshold = (uint8)Y;
		}
	}
	return Threshold;
}

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
	printf("输出为csv成功\n");
}

//二值化处理函数
void imageToBin(uint8 image_threshold_)
{
	for (int i = 0; i < image_h; i++)
	{
		for (int j = 0; j < image_w; j++)
		{
			if (image[i][j] >= image_threshold_)
			{
				bin_image[i][j] = 255;
			}
			else
			{
				bin_image[i][j] = 0;
			}
		}
	}
}

//给二值化图像周围加一圈黑边
void addBlackPoints()
{
	//把最上面两行变黑
	for (uint16 i = 0; i < image_w; i++)
	{
		bin_image[0][i] = 0;
		bin_image[1][i] = 0;
	}
	////把最下面两行变黑
	//for (uint16 i = 0; i < image_w; i++)
	//{
	//	bin_image[image_h-1][i] = 0;
	//	bin_image[image_h-2][i] = 0;
	//}
	//把最左边两行变黑
	for (uint16 i = 0; i < image_h; i++)
	{
		bin_image[i][0] = 0;
		bin_image[i][1] = 0;
	}
	//把最右边两行变黑
	for (uint16 i = 0; i < image_h; i++)
	{
		bin_image[i][image_w-1] = 0;
		bin_image[i][image_w-2] = 0;
	}
	printf("Ciallo~~~\n");
}

//找到起点，从最底部的中心开始寻找
uint32 start_point_l[2] = { 0 };
uint32 start_point_r[2] = { 0 };
uint32 start_line = image_h - 1;
uint8 getStartPoint()
{
	uint8 left_flag = 0;
	uint8 right_flag = 0;
	//初始化
	start_point_l[0] = 0;
	start_point_l[1] = 0;
	start_point_r[0] = 0;
	start_point_r[1] = 0;
	//从中间向左寻找
	for (uint32 i = image_center; i > 1; i--)  //i>1的原因是0和1两列已经被全部置0了，且若i=0或1时下面数组会越界
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i - 1] == 0 && bin_image[start_line][i - 2] == 0)
		{
			start_point_l[0] = i; //x坐标
			start_point_l[1] = start_line;
			left_flag = 1;
			break;
		}
	}
	//从中间向右寻找
	for (uint32 i = image_center; i < image_w - 2; i++)  //i < image_w-2的原因同上
	{
		if (bin_image[start_line][i] == 255 && bin_image[start_line][i + 1] == 0 && bin_image[start_line][i + 2] == 0)
		{
			start_point_r[0] = i; //x坐标
			start_point_r[1] = start_line;
			right_flag = 1;
			break;
		}
	}
	if (left_flag && right_flag) { return 1; }
	else { return 0; }
}


//八邻域
//存放点的x，y坐标
#define USE_num	image_h*4
uint16 points_l[(uint16)USE_num][2] = { 0 };//左线
uint16 points_r[(uint16)USE_num][2] = { 0 };//右线
uint16 points_c[(uint16)USE_num][2] = { 0 };//中线
uint16 data_statics_l = 0;//统计左边找到点的个数
uint16 data_statics_r = 0;//统计右边找到点的个数
uint16 dir_r[(uint16)USE_num] = { 0 };
uint16 dir_l[(uint16)USE_num] = { 0 };
#define draw_point(x,y,color,map)	circle(map, Point(x, y), 0, color, -1);  // 画半径为0的圆(画点)
//定义几种颜色
#define BLACK	Scalar(0, 0, 255)
#define WHITE	Scalar(255, 0, 0)
#define RED		Scalar(0, 0, 255)
#define RED1	Scalar(0, 17, 226)
#define BLUE	Scalar(255, 0, 0)
#define GREEN	Scalar(0, 255, 0)
#define YELLOW	Scalar(0, 255, 255)
#define BLUE1	Scalar(255, 255, 0)
#define BLUE2	Scalar(225, 169, 36)
#define BLUE3	Scalar(175, 165, 1)
#define PINK	Scalar(255, 0, 255)
#define PINK1	Scalar(67, 13, 90)
#define PINK2	Scalar(100, 13, 244)
#define PINK3	Scalar(101, 67, 256)
#define COLOR1	Scalar(154, 157, 252)
int my_abs(int value)
{
	if (value >= 0) return value;
	else return -value;
}
void search(uint16 break_flag, uint8(*image)[image_w], uint8(*bin_image)[image_w], uint16* l_stastic, uint16* r_stastic, uint8 l_start_x, uint8 l_start_y, uint8 r_start_x, uint8 r_start_y, uint8* hightest)
{

	uint8 i = 0, j = 0;

	//左边变量
	uint8 search_filds_l[8][2] = { {  0 } };
	uint8 index_l = 0;
	uint8 temp_l[8][2] = { {  0 } };
	uint8 center_point_l[2] = { {  0 } };
	uint16 l_data_statics;//统计左边
	//定义八个邻域  
	static int8 seeds_l[8][2] = { {0,  1},{-1,1},{-1,0},{-1,-1},{0,-1},{1,-1},{1,  0},{1, 1} };
	//{-1,-1},{0,-1},{+1,-1},
	//{-1, 0},	     {+1, 0},
	//{-1,+1},{0,+1},{+1,+1},
	//这个是顺时针
	static int8 square_5[24][2] = {
	{-2,-2},{-1,-2},{0,-2},{+1,-2},{+2,-2},
	{-2,-1},{-1,-1},{0,-1},{+1,-1},{+2,-1},
	{-2,-0},{-1, 0},{+1, 0},{+2,-0},
	{-2,+1},{-1,+1},{0,+1},{+1,+1},{+2,+1},
	{-2,+2},{-1,+2},{0,+2},{+1,+2},{+2,+2}
	};
	//右边变量
	uint8 search_filds_r[8][2] = { {  0 } };
	uint8 center_point_r[2] = { 0 };//中心坐标点
	uint8 index_r = 0;//索引下标
	uint8 temp_r[8][2] = { {  0 } };
	uint16 r_data_statics;//统计右边
	//定义八个邻域
	static int8 seeds_r[8][2] = { {0,  1},{1,1},{1,0}, {1,-1},{0,-1},{-1,-1}, {-1,  0},{-1, 1} };
	//{-1,-1},{0,-1},{+1,-1},
	//{-1, 0},	     {+1, 0},
	//{-1,+1},{0,+1},{+1,+1},
	//这个是逆时针
	int thres_r = 0;
	int thres_l = 0;
	l_data_statics = *l_stastic;//统计找到了多少个点，方便后续把点全部画出来
	r_data_statics = *r_stastic;//统计找到了多少个点，方便后续把点全部画出来

	//第一次更新坐标点  将找到的起点值传进来
	center_point_l[0] = l_start_x;//x
	center_point_l[1] = l_start_y;//y
	center_point_r[0] = r_start_x;//x
	center_point_r[1] = r_start_y;//y

	//开启邻域循环
	while (break_flag--)
	{
		//左边
		for (i = 0; i < 8; i++)//传递8F坐标
		{
			search_filds_l[i][0] = center_point_l[0] + seeds_l[i][0];//x
			search_filds_l[i][1] = center_point_l[1] + seeds_l[i][1];//y
			thres_l = 0;
			for (uint8 j = 0; j < 25; j++)
			{
				thres_l += image[search_filds_l[i][1] + square_5[j][1]][search_filds_l[i][0] + square_5[j][0]];
			}
			thres_l = thres_l + 201 + 200 + 200 + 201 + 200 + 199 + 196 + 194 + 192;
			thres_l /= 33;//计算局部平均阈值
			if (image[search_filds_l[i][1]][search_filds_l[i][0]] > thres_l)
			{
				bin_image[search_filds_l[i][1]][search_filds_l[i][0]] = 255;
			}
			else bin_image[search_filds_l[i][1]][search_filds_l[i][0]] = 0;

		}
		//中心坐标点填充到已经找到的点内
		points_l[l_data_statics][0] = center_point_l[0];//x
		points_l[l_data_statics][1] = center_point_l[1];//y

		uint8 x_l = 0, y_l = 0;

		x_l = center_point_l[0];
		y_l = center_point_l[1];

		l_data_statics++;//索引加一

		//右边
		for (i = 0; i < 8; i++)//传递坐标
		{
			search_filds_r[i][0] = center_point_r[0] + seeds_r[i][0];//x
			search_filds_r[i][1] = center_point_r[1] + seeds_r[i][1];//y
			thres_r = 0;
			for (uint8 j = 0; j < 25; j++)
			{
				thres_r += image[search_filds_r[i][1] + square_5[j][1]][search_filds_r[i][0] + square_5[j][0]];
			}
			thres_r = thres_r + +201 + 200 + 200 + 201 + 200 + 199 + 196 + 194 + 192;
			thres_r /= 33;//计算局部阈值
			if (image[search_filds_r[i][1]][search_filds_r[i][0]] > thres_r)
			{
				bin_image[search_filds_r[i][1]][search_filds_r[i][0]] = 255;
			}
			else bin_image[search_filds_r[i][1]][search_filds_r[i][0]] = 0;
		}
		//中心坐标点填充到已经找到的点内
		points_r[r_data_statics][0] = center_point_r[0];//x
		points_r[r_data_statics][1] = center_point_r[1];//y
		uint8 x_r = 0, y_r = 0;
		x_r = center_point_r[0];
		y_r = center_point_r[1];
		index_l = 0;//先清零，后使用
		for (i = 0; i < 8; i++)
		{
			temp_l[i][0] = 0;//先清零，后使用
			temp_l[i][1] = 0;//先清零，后使用
		}
		//左边判断
		for (i = 0; i < 8; i++)
		{
			if (bin_image[search_filds_l[i][1]][search_filds_l[i][0]] == 0
				&& bin_image[search_filds_l[(i + 1) & 7][1]][search_filds_l[(i + 1) & 7][0]] == 255)
			{
				temp_l[index_l][0] = search_filds_l[(i + 1)][0];
				temp_l[index_l][1] = search_filds_l[(i + 1)][1];
				index_l++;
				dir_l[l_data_statics - 1] = (i);//记录生长方向
			}
			if (index_l)
			{
				//更新坐标点
				center_point_l[0] = temp_l[0][0];//x
				center_point_l[1] = temp_l[0][1];//y
				for (j = 0; j < index_l; j++)
				{
					if (center_point_l[1] > temp_l[j][1])
					{
						center_point_l[0] = temp_l[j][0];//x
						center_point_l[1] = temp_l[j][1];//y
					}
				}
			}
		}
		if (my_abs(points_r[r_data_statics][0] - points_l[l_data_statics - 1][0]) < 2
			&& my_abs(points_r[r_data_statics][1] - points_l[l_data_statics - 1][1] < 2)
			)
		{
			printf("\n左右相遇退出\n");
			*hightest = (points_r[r_data_statics][1] + points_l[l_data_statics - 1][1]) >> 1;//取出最高点
			//printf("\n左右相遇退出\n");
			break;
		}
		if ((points_r[r_data_statics][1] < points_l[l_data_statics - 1][1]))
		{
			printf("\n右边比左边高了,等待一次\n");
			continue;//如果右边比左边高了
		}
		//
		//陷入O环退出
		//if (center_point_l[0] == points_l[l_data_statics - 4][0]
		//	&& center_point_l[1] == points_l[l_data_statics - 4][1])
		//{
		//	printf("o环退出\n");
		//	break;
		//}
		r_data_statics++;//索引加一

		index_r = 0;//先清零，后使用
		for (i = 0; i < 8; i++)
		{
			temp_r[i][0] = 0;//先清零，后使用
			temp_r[i][1] = 0;//先清零，后使用
		}

		//右边判断
		for (i = 0; i < 8; i++)
		{
			if (bin_image[search_filds_r[i][1]][search_filds_r[i][0]] == 0
				&& bin_image[search_filds_r[(i + 1) & 7][1]][search_filds_r[(i + 1) & 7][0]] == 255)
			{
				temp_r[index_r][0] = search_filds_r[(i + 1)][0];
				temp_r[index_r][1] = search_filds_r[(i + 1)][1];
				index_r++;//索引加一
				dir_r[r_data_statics - 1] = (i);//记录生长方向
				//printf("dir[%d]:%d\n", r_data_statics - 1, dir_r[r_data_statics - 1]);
			}
			if (index_r)
			{
				//更新坐标点
				center_point_r[0] = temp_r[0][0];//x
				center_point_r[1] = temp_r[0][1];//y
				for (j = 0; j < index_r; j++)
				{
					if (center_point_r[1] > temp_r[j][1])
					{
						center_point_r[0] = temp_r[j][0];//x
						center_point_r[1] = temp_r[j][1];//y
					}
				}

			}
		}


		draw_point(center_point_r[0], center_point_r[1], PINK2, img);//显示起点
		draw_point(x_r, y_r, GREEN, img);//显示起点
		draw_point(center_point_l[0], center_point_l[1], PINK2, img);//显示起点
		draw_point(x_l, y_l, BLUE, img);//显示起点

		namedWindow("ouput", WINDOW_FREERATIO);
		imshow("ouput", img);
		waitKey(10);

	}
}



//处理图像函数
void processImage()
{
	//处理图像
	Data_trans(img);
	image_filter();
	//使用otsu确定阈值
	uint8 image_threshold = otsuThreshold(image[0], image_h, image_w);
	//二值化处理
	imageToBin(image_threshold);
	//给图像加黑边
	addBlackPoints();

}

int main()
{
    //namedWindow("输入图像");
    imshow("img", img);
 //   //处理图像
 //   Data_trans(img);
	//image_filter();
 //   //使用otsu确定阈值
	//uint8 image_threshold = otsuThreshold(image[0], image_h, image_w); //第一个参数不理解
	//printf("----------------%d----------------", image_threshold);
	////二值化处理
	//imageToBin(image_threshold);
	////给图像加黑边
	//addBlackPoints();
	processImage();
	//打印图像至csv
	printToCSV(image, "output_original.csv");
	printToCSV(bin_image, "output.csv");
	//准备八邻域
	uint8 hightest = 0;
	if (getStartPoint())  //先获得起点
	{
		printf("Start point found.");
		search((uint16)USE_num, image, bin_image, &data_statics_l, &data_statics_r, start_point_l[0], start_point_l[1], start_point_r[0], start_point_r[1], &hightest);
	}


	//画出边界点
	for (uint32 i = 0; i < data_statics_l; i++)
	{
		//printf("dir_l[%d]:%d\n", i,dir_l[i]);
		draw_point(points_l[i][0], points_l[i][1], PINK2, img);//显示起点
	}
	for (uint32 i = 0; i < data_statics_r; i++)
	{
		//printf("dir_r[%d]:%d\n", i,dir_r[i]);
		draw_point(points_r[i][0], points_r[i][1], PINK, img);//显示起点
	}
	namedWindow("ouput2", WINDOW_FREERATIO);
	imshow("ouput2", img);

	//显示输出图像
	for (uint8 h = 0; h < image_h - 1; h++)
	{
		for (uint8 w = 0; w < image_w - 1; w++)
		{
			if (bin_image[h][w] == 255)
			{
				img.at<Vec3b>(h, w)[0] = 255;
				img.at<Vec3b>(h, w)[1] = 255;
				img.at<Vec3b>(h, w)[2] = 255;
			}
			else
			{
				img.at<Vec3b>(h, w)[0] = 0;
				img.at<Vec3b>(h, w)[1] = 0;
				img.at<Vec3b>(h, w)[2] = 0;
			}

		}
	}
	namedWindow("bin", WINDOW_FREERATIO);
	imshow("bin", img);


    waitKey();
    destroyAllWindows();


 
    return 0;
}
