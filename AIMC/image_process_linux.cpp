#include "image_process.h"
#include <fstream>
#include <vector>
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
#define ROAD  "/home/anysets/Documents/AIMC-CV/AIMC/frame_40.png"
#define VIDEO  "/home/anysets/Documents/Video/3.mp4"
// #define ROAD "E:\\Files\\OpenCV++\\test.BMP"
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 94
//ע�����border_length����̫�࣬�ᵼ������Խ��
#define border_length 100
uint8_t image_center = image_w / 2;

//��ͼ�������ӡ��csv�ļ�--------------------------------------------------------
void printToCSV(uint8_t matrix[][image_w], const std::string& filename) {
	std::ofstream outputFile(filename);

	if (!outputFile.is_open()) {
		std::cerr << "�޷����ļ���" << filename << std::endl;
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
    printf("OTSU=%d\n", image_threshold);

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
        image[image_h - 1][i] = 0;  //��ͼ�������滭һ�кڿ��Խ������Խ������
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
	// cout << "Found Left Start Point: " << start_point_l[0] << ", " << start_point_l[1] << endl;
    // cout << "Found Right Start Point: " << start_point_r[0] << ", " << start_point_r[1] << endl;
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
	else { return 0; }
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
    border_count_left = 0;
    growth_count_left = 0;
    border_count_right = 0;
    growth_count_right = 0;
    for (int j = 0; j < border_length; j++)
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
                // printf("New Right Start Point: %d, %d\n", start_point_r[0], start_point_r[1]);
				// printf("growth_direction_right: %d\n", growth_direction_right[j]);
                // printf("finished\n");
				// printf("%d\n", j);
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
uint8_t border_location_left_optimized[image_h] = { 0 };
uint8_t border_location_right_optimized[image_h] = { 0 };
void optimizeBorder()
{
    memset(border_location_left_optimized, 0, sizeof(border_location_left_optimized));
    memset(border_location_right_optimized, 0, sizeof(border_location_right_optimized));

    // ��߽�Ĳ���
    // ѭ���Ѿ���ȡ���ı߽磬ȥ������Ҫ�ĵ�
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        // border_location_left[i][0];  //һ��ֻҪһ�������ˣ��е�������[0]����߽���Ҫ����ÿ�����Ҳ�ĵ�
        // ��ʼ�ж�
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
    // �ұ߽�Ĳ���
    // ѭ���Ѿ���ȡ���ı߽磬ȥ������Ҫ�ĵ�
    for (int i = 0; i < border_length; i++)
    {
        // printf("%d\n", i);
        //һ��ֻҪһ�������ˣ��е�������[0]���ұ߽���Ҫ����ÿ�������ĵ�
        // ��ʼ�ж�
        // �ұ߽���Ҫ��һ����ʼ��
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

// ��ȡ����
int32_t middle_line[border_length] = { 0 }; 
void getMiddleLine()
{
    memset(middle_line, 0, sizeof(middle_line));
    for (int i = 0; i < image_h; i++)
    {
        middle_line[i] = (border_location_left_optimized[i] + border_location_right_optimized[i]) / 2;
    }
}

// ��ȡƫ��ֵ
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
    turn_to_bin();  //ͼ���ֵ��
    image_compress(image,bin_image); //ͼ��ѹ��
    image_filter();
    image_draw_rectan(bin_image);
    if (getStartPoint())  //�Ȼ�����
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
        draw_point(border_location_left_optimized[i], i, GREEN, Binimage_color);
        draw_point(border_location_right_optimized[i], i, GREEN, Binimage_color);
        // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
        // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
        draw_point(middle_line[i], i, RED, Binimage_color);
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
    */
    
    //��Ƶ����
    VideoCapture cap(VIDEO);
    // �����Ƶ�Ƿ�ɹ���
    if (!cap.isOpened()) {
        std::cerr << "Error opening video file" << std::endl;
        return -1;
    }
    Mat frame;
    int frame_number = 0;
    // ѭ����ȡÿһ֡
    while (cap.read(frame)) {
        printf("------------------------------------------frame:%d------------------------------------------\n", frame_number);
         //����ÿһ֡ΪͼƬ
        string filename = "frame_" + to_string(frame_number) + ".png";
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
            draw_point(border_location_left_optimized[i], i, GREEN, Binimage_color);
            draw_point(border_location_right_optimized[i], i, GREEN, Binimage_color);
            // draw_point(border_location_left[i][1], border_location_left[i][0], PINK, Binimage_color);
            // draw_point(border_location_right[i][1], border_location_right[i][0], PINK, Binimage_color);
            draw_point(middle_line[i], i, RED, Binimage_color);
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
        frame_number++;
    }
    
	return 0;
}
