#include "image_process.h"
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
#define IMAGE_H 120
#define IMAGE_W 188
#define image_h 60
#define image_w 120
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
    printf("��ֵ=%d\n", image_threshold);

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
    }
}
void image_process(void)
{
    turn_to_bin();  //ͼ���ֵ��
    image_compress(image,bin_image); //ͼ��ѹ��
    image_draw_rectan(bin_image);
}
int main(){
	const string ROAD = "E:\\AR\\image\\y1.BMP";
	Mat roadImage = imread(ROAD);
	if (roadImage.empty())
    {
        cerr << "Failed to load image!" << endl;
        return -1;
    }
	dataTrans(roadImage, image);

    image_process();
    Mat Binimage(image_h, image_w, CV_8UC1);
    for (int i = 0; i < image_h; i++) {
        for (int j = 0; j < image_w; j++) {
            Binimage.at<uint8_t>(i, j) = bin_image[i][j];
        }
    }
    // ��ʾͼ��
    namedWindow("��ֵ��");
    imshow("��ֵ��", Binimage);

	namedWindow("output", WINDOW_NORMAL);
	imshow("output", roadImage);
	resizeWindow("output", 360, 240);
	waitKey(0);
	destroyAllWindows();
	return 0;
}