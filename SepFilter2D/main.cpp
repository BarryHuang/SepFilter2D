#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))

void sepFilter2D(const float *src, float *tempBuf, float *dst, int ksize, int width, int height, float *kx, float *ky);
void sepFilter2D_v2(const float *src, float *tempBuf, float *dst, int ksize, int width, int height, float *kx, float *ky);
void sepFilter2D_v3(const float *src, float *tempBuf, float *dst, int ksize, int width, int height, float *kx, float *ky);

int main()
{
    printf("test\n");

    const int max_iter = 5000;
    const int width = 120;
    const int height = 80;
    const int ksize = 9;

    float *src = (float *)malloc(width * height * sizeof(float));
    float *dst1 = (float *)malloc(width * height * sizeof(float));
    float *dst2 = (float *)malloc(width * height * sizeof(float));
    float *tmp = (float *)malloc(width * (height + ksize - 1) * sizeof(float));
    float *kx = (float *)malloc(ksize * sizeof(float));
    float *ky = (float *)malloc(ksize * sizeof(float));

    memset(kx, 0, ksize * sizeof(float));
    memset(ky, 0, ksize * sizeof(float));

    kx[0] = -1.f;
    kx[ksize-1] = 1.f;

    ky[0] = 3.f / 16.f;
    ky[ksize/2] = 10.f / 16.f;
    ky[ksize-1] = 3.f / 16.f;

    for (int i = 0; i < width * height; i++)
        src[i] = (float)(rand() % 256) / 255.f;

    float time1, time2, time3;
    {
        clock_t t1 = clock();
        for (int iter = 0; iter < max_iter; iter++)
        {
            sepFilter2D(src, tmp, dst1, ksize, width, height, kx, ky);
            sepFilter2D(src, tmp, dst1, ksize, width, height, ky, kx);
        }
            
        clock_t t2 = clock();
        time1 = (float)(t2 - t1) / CLOCKS_PER_SEC * 1000.f;
    }

    {
        clock_t t1 = clock();
        for (int iter = 0; iter < max_iter; iter++)
        {
            sepFilter2D_v2(src, tmp, dst2, ksize, width, height, kx, ky);
            sepFilter2D_v2(src, tmp, dst2, ksize, width, height, ky, kx);
        }
            
        clock_t t2 = clock();
        time2 = (float)(t2 - t1) / CLOCKS_PER_SEC * 1000.f;
    }

    {
        clock_t t1 = clock();
        for (int iter = 0; iter < max_iter; iter++)
        {
            sepFilter2D_v3(src, tmp, dst2, ksize, width, height, kx, ky);
            sepFilter2D_v3(src, tmp, dst2, ksize, width, height, ky, kx);
        }

        clock_t t2 = clock();
        time3 = (float)(t2 - t1) / CLOCKS_PER_SEC * 1000.f;
    }

    {
        for (int i = 0; i < width * height; i++)
            if (fabsf(dst1[i] - dst2[i]) >= 1.f/255.f)
            {
                printf("Not bit true %f - %f\n", dst1[i], dst2[i]);
                break;
            }
                
    }
    
    printf("sepFilter2D %d times, processing time: %.2f ms \n", max_iter, time1);
    printf("sepFilter2D_v2 %d times, processing time: %.2f ms, speed up %f \n", max_iter, time2, time1 / time2);
    printf("sepFilter2D_v3 %d times, processing time: %.2f ms, speed up %f \n", max_iter, time3, time1 / time3);

    free(src);
    free(dst1);
    free(dst2);
    free(tmp);
    free(kx);
    free(ky);

    system("PAUSE");

    return 0;
}

void sepFilter2D_v2(
    const float *src,
    float *tempBuf,
    float *dst,
    int ksize,
    int width,
    int height,
    float *kx,
    float *ky)
{
    const int iBias = (int)(ksize / 2.f);

    int i = 0, j = 0, k = 0;
    int iTmpX, iTmpY;
    int img_index;
    float fSum;

    // Horizantal
    const float *src_scan = src;
    float *tmp_scan = tempBuf;
    for (i = 0; i < height; i++)
    {
        j = 0;
        // 0 ~ iBias
        for (; j < iBias; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpX = j + k - iBias;
                iTmpX = MAX(iTmpX, 0);

                fSum += src_scan[iTmpX] * kx[k];
            }
            tmp_scan[j] = fSum;
        }

        // iBias ~ width - iBias
        for (; j < width - iBias; j++)
        {
            fSum = 0.f;
            iTmpX = j - iBias;
            for (k = 0; k < ksize; k++)
            {
                fSum += src_scan[iTmpX] * kx[k];
                iTmpX++;
            }
            tmp_scan[j] = fSum;
        }

        // width - iBias ~ width
        for (; j < width; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpX = j + k - iBias;
                iTmpX = MIN(iTmpX, width - 1);

                fSum += src_scan[iTmpX] * kx[k];
            }
            tmp_scan[j] = fSum;
        }

        src_scan += width;
        tmp_scan += width;
    }

    // Vertical
    float *dst_scan = dst;
    for (i = 0; i < iBias; i++)
    {
        for (j = 0; j < width; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpY = i + k - iBias;
                iTmpY = MAX(iTmpY, 0);

                img_index = iTmpY * width + j;

                fSum += tempBuf[img_index] * ky[k];
            }
            dst_scan[j] = fSum;
        }
        dst_scan += width;
    }

    float *tmp_k0_scan = tempBuf;
    for (; i < height - iBias; i++)
    {
        for (j = 0; j < width; j++)
        {
            fSum = 0.f;
            tmp_scan = tmp_k0_scan;
            for (k = 0; k < ksize; k++)
            {
                fSum += tmp_scan[j] * ky[k];
                tmp_scan += width;
            }
            dst_scan[j] = fSum;
        }
        dst_scan += width;
        tmp_k0_scan += width;
    }

    for (; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpY = i + k - iBias;
                iTmpY = MIN(iTmpY, height - 1);

                img_index = iTmpY * width + j;

                fSum += tempBuf[img_index] * ky[k];
            }

            dst_scan[j] = fSum;
        }
        dst_scan += width;
    }
}

void sepFilter2D_v3(
    const float *src,
    float *tempBuf,
    float *dst,
    int ksize,
    int width,
    int height,
    float *kx,
    float *ky)
{
    const int iBias = (int)(ksize / 2.f);

    int i = 0, j = 0, k = 0;
    int iTmpX, iTmpY;
    int img_index;
    float fSum;

    // Horizantal
    const float *src_scan = src;
    float *tmp_scan = tempBuf + iBias * width;
    for (i = 0; i < height; i++)
    {
        j = 0;
        // 0 ~ iBias
        for (; j < iBias; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpX = j + k - iBias;
                iTmpX = MAX(iTmpX, 0);

                fSum += src_scan[iTmpX] * kx[k];
            }
            tmp_scan[j] = fSum;
        }

        // iBias ~ width - iBias
        for (; j < width - iBias; j++)
        {
            fSum = 0.f;
            iTmpX = j - iBias;
            for (k = 0; k < ksize; k++)
            {
                fSum += src_scan[iTmpX] * kx[k];
                iTmpX++;
            }
            tmp_scan[j] = fSum;
        }

        // width - iBias ~ width
        for (; j < width; j++)
        {
            fSum = 0.f;
            for (k = 0; k < ksize; k++)
            {
                iTmpX = j + k - iBias;
                iTmpX = MIN(iTmpX, width - 1);

                fSum += src_scan[iTmpX] * kx[k];
            }
            tmp_scan[j] = fSum;
        }

        src_scan += width;
        tmp_scan += width;
    }

    // Copy padding

    for (i = 0; i < iBias; i++)
    {
        memcpy(tempBuf + i * width, tempBuf + iBias * width, width * sizeof(float));
        memcpy(tempBuf + (height + iBias + i) * width, tempBuf + (height - 1 + iBias) * width, width * sizeof(float));
    }

    // Vertical
    float *dst_scan = dst;
    float *tmp_k0_scan = tempBuf;
    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            fSum = 0.f;
            tmp_scan = tmp_k0_scan;
            for (k = 0; k < ksize; k++)
            {
                fSum += tmp_scan[j] * ky[k];
                tmp_scan += width;
            }
            dst_scan[j] = fSum;
        }
        dst_scan += width;
        tmp_k0_scan += width;
    }
}

void sepFilter2D(
    const float *src,
    float *tempBuf,
    float *dst,
    int ksize,
    int width,
    int height,
    float *kx,
    float *ky)
{
    int i = 0, j = 0, k = 0;
    int iTmpX, iTmpY;
    int k_index;
    int img_index;
    float fSum;
    int iBias = (int)(ksize / 2.f);

    for (i = 0; i < height; i++)
    {
        for (j = 0; j < width; j++)
        {
            fSum = 0;
            for (k = j - iBias; k <= j + iBias; k++)
            {
                k_index = k - (j - iBias);

                iTmpX = MAX(k, 0);
                iTmpX = MIN(iTmpX, width - 1);
                iTmpY = i;
                img_index = iTmpY * width + iTmpX;

                fSum += src[img_index] * kx[k_index];
            }
            tempBuf[i * width + j] = fSum;
        }
    }

    for (j = 0; j < width; j++)
    {
        for (i = 0; i < height; i++)
        {
            fSum = 0;
            for (k = i - iBias; k <= i + iBias; k++)
            {
                k_index = k - (i - iBias);

                iTmpX = j;
                iTmpY = MAX(k, 0);
                iTmpY = MIN(iTmpY, height - 1);

                img_index = iTmpY * width + iTmpX;

                fSum += tempBuf[img_index] * ky[k_index];
            }
            dst[i * width + j] = fSum;
        }
    }
}