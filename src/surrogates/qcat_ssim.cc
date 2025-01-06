

// #define K1 0.01
// #define K2 0.03

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <qcat_ssim.h>

/////////////////// 1D

double SSIM_1d_calcWindow_float(float* data, float* other, int offset0, int windowSize0)
{
    int i0;
    int np=0; //Number of points

    float xMin=data[offset0];
    float xMax=data[offset0];
    float yMin=other[offset0];
    float yMax=other[offset0];
    double xSum=0;
    double ySum=0;

    for(i0=offset0; i0<offset0+windowSize0; i0++) {
        np++;
        if(xMin>data[i0])
            xMin=data[i0];
        if(xMax<data[i0])
            xMax=data[i0];
        if(yMin>other[i0])
            yMin=other[i0];
        if(yMax<other[i0])
            yMax=other[i0];
        xSum+=data[i0];
        ySum+=other[i0];
    }


    double xMean=xSum/np;
    double yMean=ySum/np;

    double var_x = 0, var_y = 0, var_xy = 0;

    for(i0=offset0; i0<offset0+windowSize0; i0++) {
        var_x += (data[i0] - xMean)*(data[i0] - xMean);
        var_y += (other[i0] - yMean)*(other[i0] - yMean);
        var_xy += (data[i0] - xMean)*(other[i0] - yMean);
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;

    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    return ssim;
}

double SSIM_1d_calcWindow_double(double* oriData, double* decData, int offset0, int windowSize0)
{
    int i0;
    int np=0; //Number of points

    double* data = oriData;
    double* other = decData;

    double xMin=data[offset0];
    double xMax=data[offset0];
    double yMin=other[offset0];
    double yMax=other[offset0];
    double xSum=0;
    double ySum=0;

    for(i0=offset0; i0<offset0+windowSize0; i0++) {
        np++;
        if(xMin>data[i0])
            xMin=data[i0];
        if(xMax<data[i0])
            xMax=data[i0];
        if(yMin>other[i0])
            yMin=other[i0];
        if(yMax<other[i0])
            yMax=other[i0];
        xSum+=data[i0];
        ySum+=other[i0];
    }


    double xMean=xSum/np;
    double yMean=ySum/np;

    double var_x = 0, var_y = 0, var_xy = 0;

    for(i0=offset0; i0<offset0+windowSize0; i0++) {
        var_x += (data[i0] - xMean)*(data[i0] - xMean);
        var_y += (other[i0] - yMean)*(other[i0] - yMean);
        var_xy += (data[i0] - xMean)*(other[i0] - yMean);
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;

    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    return ssim;
}

double SSIM_1d_windowed_float(float* oriData, float* decData, size_t size0, int windowSize0, int windowShift0) {
    int offset0;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }

    //offsetInc0=windowSize0/2;
    offsetInc0=windowShift0;


    for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
        nw++;
        ssimSum+=SSIM_1d_calcWindow_float(oriData, decData, offset0, windowSize0);
    }

    return ssimSum/nw;
}

double SSIM_1d_windowed_double(double* oriData, double* decData, size_t size0, int windowSize0, int windowShift0)
{
    int offset0;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }

    //offsetInc0=windowSize0/2;
    offsetInc0=windowShift0;


    for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
        nw++;
        ssimSum+=SSIM_1d_calcWindow_double(oriData, decData, offset0, windowSize0);
    }

    return ssimSum/nw;
}

//////////////////// 2D

double SSIM_2d_windowed_float(float* oriData, float* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1)
{
    int offset0,offset1;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1;

    float* data = oriData;
    float* other = decData;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;

    for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

        for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
            nw++;
            double ssim = SSIM_2d_calcWindow_float(data, other, size0, offset0, offset1, windowSize0, windowSize1);
            ssimSum+=ssim;
        }
    }

    return ssimSum/nw;
}

double SSIM_2d_calcWindow_float(float* data, float *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1)
{
    int i0,i1,index;
    int np=0; //Number of points
    float xMin=data[offset0+size0*offset1];
    float xMax=data[offset0+size0*offset1];
    float yMin=other[offset0+size0*offset1];
    float yMax=other[offset0+size0*offset1];
    double xSum=0;
    double ySum=0;

    for(i1=offset1; i1<offset1+windowSize1; i1++) {
        for(i0=offset0; i0<offset0+windowSize0; i0++) {
            np++;
            index=i0+size0*i1;
            if(xMin>data[index])
                xMin=data[index];
            if(xMax<data[index])
                xMax=data[index];
            if(yMin>other[index])
                yMin=other[index];
            if(yMax<other[index])
                yMax=other[index];
            xSum+=data[index];
            ySum+=other[index];
        }
    }

    double xMean=xSum/np;
    double yMean=ySum/np;

    double var_x = 0, var_y = 0, var_xy = 0;

    for(i1=offset1; i1<offset1+windowSize1; i1++) {
        for(i0=offset0; i0<offset0+windowSize0; i0++) {
            index=i0+size0*i1;
            var_x += (data[index] - xMean)*(data[index] - xMean);
            var_y += (other[index] - yMean)*(other[index] - yMean);
            var_xy += (data[index] - xMean)*(other[index] - yMean);
        }
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;


    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    return ssim;
}

double SSIM_2d_windowed_double(double* oriData, double* decData, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowShift0, int windowShift1)
{
    int offset0,offset1;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1;

    double* data = oriData;
    double* other = decData;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;

    for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

        for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
            nw++;
            ssimSum+=SSIM_2d_calcWindow_double(data, other, size0, offset0, offset1, windowSize0, windowSize1);
        }
    }

    return ssimSum/nw;
}

double SSIM_2d_calcWindow_double(double* data, double *other, size_t size0, int offset0, int offset1, int windowSize0, int windowSize1) {
    int i0,i1,index;
    int np=0; //Number of points
    double xMin=data[offset0+size0*offset1];
    double xMax=data[offset0+size0*offset1];
    double yMin=other[offset0+size0*offset1];
    double yMax=other[offset0+size0*offset1];
    double xSum=0;
    double ySum=0;

    for(i1=offset1; i1<offset1+windowSize1; i1++) {
        for(i0=offset0; i0<offset0+windowSize0; i0++) {
            np++;
            index=i0+size0*i1;
            if(xMin>data[index])
                xMin=data[index];
            if(xMax<data[index])
                xMax=data[index];
            if(yMin>other[index])
                yMin=other[index];
            if(yMax<other[index])
                yMax=other[index];
            xSum+=data[index];
            ySum+=other[index];
        }
    }

    double xMean=xSum/np;
    double yMean=ySum/np;

    double var_x = 0, var_y = 0, var_xy = 0;

    for(i1=offset1; i1<offset1+windowSize1; i1++) {
        for(i0=offset0; i0<offset0+windowSize0; i0++) {
            index=i0+size0*i1;
            var_x += (data[index] - xMean)*(data[index] - xMean);
            var_y += (other[index] - yMean)*(other[index] - yMean);
            var_xy += (data[index] - xMean)*(other[index] - yMean);
        }
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;

    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    
    if(ssim<0)
		printf("ssim=%f", ssim);

    return ssim;
}


//////////////////////// 3D

double SSIM_3d_windowed_float(float* oriData, float* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2)
{
    int offset0,offset1,offset2;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1,offsetInc2;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }
    if(windowSize2>size2) {
        printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    //offsetInc2=windowSize2/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;
    offsetInc2=windowShift2;

    for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2) { //MOVING WINDOW

        for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

            for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
                nw++;
                ssimSum+=SSIM_3d_calcWindow_float(oriData, decData, size1, size0, offset0, offset1, offset2, windowSize0, windowSize1, windowSize2);

            }
        }
    }

    return ssimSum/nw;
}

double SSIM_3d_calcWindow_float(float* data, float* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2) {
    int i0,i1,i2,index;
    int np=0; //Number of points
    float xMin=data[offset0+size0*(offset1+size1*offset2)];
    float xMax=data[offset0+size0*(offset1+size1*offset2)];
    float yMin=other[offset0+size0*(offset1+size1*offset2)];
    float yMax=other[offset0+size0*(offset1+size1*offset2)];
    double xSum=0;
    double ySum=0;

    for(i2=offset2; i2<offset2+windowSize2; i2++) {
        for(i1=offset1; i1<offset1+windowSize1; i1++) {
            for(i0=offset0; i0<offset0+windowSize0; i0++) {
                np++;
                index=i0+size0*(i1+size1*i2);
                if(xMin>data[index])
                    xMin=data[index];
                if(xMax<data[index])
                    xMax=data[index];
                if(yMin>other[index])
                    yMin=other[index];
                if(yMax<other[index])
                    yMax=other[index];
                xSum+=data[index];
                ySum+=other[index];
            }
        }
    }


    double xMean=xSum/np;
    double yMean=ySum/np;
    double var_x = 0, var_y = 0, var_xy = 0;

    for(i2=offset2; i2<offset2+windowSize2; i2++) {
        for(i1=offset1; i1<offset1+windowSize1; i1++) {
            for(i0=offset0; i0<offset0+windowSize0; i0++) {
                index=i0+size0*(i1+size1*i2);
                var_x += (data[index] - xMean)*(data[index] - xMean);
                var_y += (other[index] - yMean)*(other[index] - yMean);
                var_xy += (data[index] - xMean)*(other[index] - yMean);
            }
        }
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;


    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;

    return ssim;
}


double SSIM_3d_windowed_double(double* oriData, double* decData, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowShift0, int windowShift1, int windowShift2)
{
    int offset0,offset1,offset2;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1,offsetInc2;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }
    if(windowSize2>size2) {
        printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    //offsetInc2=windowSize2/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;
    offsetInc2=windowShift2;

    for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2) { //MOVING WINDOW

        for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

            for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
                nw++;
                ssimSum+=SSIM_3d_calcWindow_double(oriData, decData, size1, size0, offset0, offset1, offset2, windowSize0, windowSize1, windowSize2);

            }
        }
    }

    return ssimSum/nw;
}

double SSIM_3d_calcWindow_double(double* data, double* other, size_t size1, size_t size0, int offset0, int offset1, int offset2, int windowSize0, int windowSize1, int windowSize2)
{
    int i0,i1,i2,index;
    int np=0; //Number of points
    double xMin=data[offset0+size0*(offset1+size1*offset2)];
    double xMax=data[offset0+size0*(offset1+size1*offset2)];
    double yMin=other[offset0+size0*(offset1+size1*offset2)];
    double yMax=other[offset0+size0*(offset1+size1*offset2)];
    double xSum=0;
    double ySum=0;

    for(i2=offset2; i2<offset2+windowSize2; i2++) {
        for(i1=offset1; i1<offset1+windowSize1; i1++) {
            for(i0=offset0; i0<offset0+windowSize0; i0++) {
                np++;
                index=i0+size0*(i1+size1*i2);
                if(xMin>data[index])
                    xMin=data[index];
                if(xMax<data[index])
                    xMax=data[index];
                if(yMin>other[index])
                    yMin=other[index];
                if(yMax<other[index])
                    yMax=other[index];
                xSum+=data[index];
                ySum+=other[index];
            }
        }
    }

    double xMean=xSum/np;
    double yMean=ySum/np;
    double var_x = 0, var_y = 0, var_xy = 0;

    for(i2=offset2; i2<offset2+windowSize2; i2++) {
        for(i1=offset1; i1<offset1+windowSize1; i1++) {
            for(i0=offset0; i0<offset0+windowSize0; i0++) {
                index=i0+size0*(i1+size1*i2);
                var_x += (data[index] - xMean)*(data[index] - xMean);
                var_y += (other[index] - yMean)*(other[index] - yMean);
                var_xy += (data[index] - xMean)*(other[index] - yMean);
            }
        }
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;


    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;

    return ssim;
}

double SSIM_4d_windowed_float(float* oriData, float* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3)
{
    int offset0,offset1,offset2,offset3;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1,offsetInc2,offsetInc3;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }
    if(windowSize2>size2) {
        printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);
    }
    if(windowSize3>size3) {
        printf("ERROR: windowSize3 = %d > %zu\n", windowSize3, size3);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    //offsetInc2=windowSize2/2;
    //offsetInc3=windowSize3/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;
    offsetInc2=windowShift2;
    offsetInc3=windowShift3;

    for(offset3=0; offset3+windowSize3<=size3; offset3+=offsetInc3) { //MOVING WINDOW

        for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2) { //MOVING WINDOW

            for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

                for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
                    nw++;
                    ssimSum+=SSIM_4d_calcWindow_float(oriData, decData, size2, size1, size0, offset0, offset1, offset2, offset3, windowSize0, windowSize1, windowSize2, windowSize3);
                }
            }
        }
    }

    return ssimSum/nw;
    return 0;
}

double SSIM_4d_calcWindow_float(float* data, float* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3)
{
    int i0,i1,i2,i3,index;
    int np=0; //Number of points
    float xMin=data[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    float xMax=data[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    float yMin=other[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    float yMax=other[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    double xSum=0;
    double ySum=0;
    for(i3=offset3; i3<offset3+windowSize3; i3++) {
        for(i2=offset2; i2<offset2+windowSize2; i2++) {
            for(i1=offset1; i1<offset1+windowSize1; i1++) {
                for(i0=offset0; i0<offset0+windowSize0; i0++) {
                    np++;
                    index=i0+size0*(i1+size1*(i2+size2*i3));
                    if(xMin>data[index])
                        xMin=data[index];
                    if(xMax<data[index])
                        xMax=data[index];
                    if(yMin>other[index])
                        yMin=other[index];
                    if(yMax<other[index])
                        yMax=other[index];
                    xSum+=data[index];
                    ySum+=other[index];
                }
            }
        }
    }

    double xMean=xSum/np;
    double yMean=ySum/np;
    double var_x = 0, var_y = 0, var_xy = 0;

    for(i3=offset3; i3<offset3+windowSize3; i3++) {
        for(i2=offset2; i2<offset2+windowSize2; i2++) {
            for(i1=offset1; i1<offset1+windowSize1; i1++) {
                for(i0=offset0; i0<offset0+windowSize0; i0++) {
                    index=i0+size0*(i1+size1*(i2+size2*i3));
                    var_x += (data[index] - xMean)*(data[index] - xMean);
                    var_y += (other[index] - yMean)*(other[index] - yMean);
                    var_xy += (data[index] - xMean)*(other[index] - yMean);
                }
            }
        }
    }

    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;

    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    return ssim;
}

double SSIM_4d_windowed_double(double* oriData, double* decData, size_t size3, size_t size2, size_t size1, size_t size0, int windowSize0, int windowSize1, int windowSize2, int windowSize3, int windowShift0, int windowShift1, int windowShift2, int windowShift3)
{
    int offset0,offset1,offset2,offset3;
    int nw=0; //Number of windows
    double ssimSum=0;
    int offsetInc0,offsetInc1,offsetInc2,offsetInc3;

    if(windowSize0>size0) {
        printf("ERROR: windowSize0 = %d > %zu\n", windowSize0, size0);
    }
    if(windowSize1>size1) {
        printf("ERROR: windowSize1 = %d > %zu\n", windowSize1, size1);
    }
    if(windowSize2>size2) {
        printf("ERROR: windowSize2 = %d > %zu\n", windowSize2, size2);
    }
    if(windowSize3>size3) {
        printf("ERROR: windowSize3 = %d > %zu\n", windowSize3, size3);
    }

    //offsetInc0=windowSize0/2;
    //offsetInc1=windowSize1/2;
    //offsetInc2=windowSize2/2;
    //offsetInc3=windowSize3/2;
    offsetInc0=windowShift0;
    offsetInc1=windowShift1;
    offsetInc2=windowShift2;
    offsetInc3=windowShift3;

    for(offset3=0; offset3+windowSize3<=size3; offset3+=offsetInc3) { //MOVING WINDOW

        for(offset2=0; offset2+windowSize2<=size2; offset2+=offsetInc2) { //MOVING WINDOW

            for(offset1=0; offset1+windowSize1<=size1; offset1+=offsetInc1) { //MOVING WINDOW

                for(offset0=0; offset0+windowSize0<=size0; offset0+=offsetInc0) { //MOVING WINDOW
                    nw++;
                    ssimSum+=SSIM_4d_calcWindow_double(oriData, decData, size2, size1, size0, offset0, offset1, offset2, offset3, windowSize0, windowSize1, windowSize2, windowSize3);
                }
            }
        }
    }

    return ssimSum/nw;
    return 0;
}

double SSIM_4d_calcWindow_double(double* data, double* other, size_t size2, size_t size1, size_t size0, int offset0, int offset1, int offset2, int offset3,int windowSize0, int windowSize1, int windowSize2, int windowSize3)
{
    int i0,i1,i2,i3,index;
    int np=0; //Number of points
    double xMin=data[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    double xMax=data[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    double yMin=other[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    double yMax=other[offset0+size0*(offset1+size1*(offset2+size2*offset3))];
    double xSum=0;
    double ySum=0;
    for(i3=offset3; i3<offset3+windowSize3; i3++) {
        for(i2=offset2; i2<offset2+windowSize2; i2++) {
            for(i1=offset1; i1<offset1+windowSize1; i1++) {
                for(i0=offset0; i0<offset0+windowSize0; i0++) {
                    np++;
                    index=i0+size0*(i1+size1*(i2+size2*i3));
                    if(xMin>data[index])
                        xMin=data[index];
                    if(xMax<data[index])
                        xMax=data[index];
                    if(yMin>other[index])
                        yMin=other[index];
                    if(yMax<other[index])
                        yMax=other[index];
                    xSum+=data[index];
                    ySum+=other[index];
                }
            }
        }
    }

    double xMean=xSum/np;
    double yMean=ySum/np;
    double var_x = 0, var_y = 0, var_xy = 0;

    for(i3=offset3; i3<offset3+windowSize3; i3++) {
        for(i2=offset2; i2<offset2+windowSize2; i2++) {
            for(i1=offset1; i1<offset1+windowSize1; i1++) {
                for(i0=offset0; i0<offset0+windowSize0; i0++) {
                    index=i0+size0*(i1+size1*(i2+size2*i3));
                    var_x += (data[index] - xMean)*(data[index] - xMean);
                    var_y += (other[index] - yMean)*(other[index] - yMean);
                    var_xy += (data[index] - xMean)*(other[index] - yMean);
                }
            }
        }
    }
    var_x /= np;
    var_y /= np;
    var_xy /= np;

    double xSigma=sqrt(var_x);
    double ySigma=sqrt(var_y);
    double xyCov = var_xy;

    double c1,c2;
    if(xMax-xMin==0) {
        c1=K1*K1;
        c2=K2*K2;
    } else {
        c1=K1*K1*(xMax-xMin)*(xMax-xMin);
        c2=K2*K2*(xMax-xMin)*(xMax-xMin);
    }
    double c3=c2/2;

    double luminance=(2*xMean*yMean+c1)/(xMean*xMean+yMean*yMean+c1);
    double contrast=(2*xSigma*ySigma+c2)/(xSigma*xSigma+ySigma*ySigma+c2);
    double structure=(xyCov+c3)/(xSigma*ySigma+c3);
    double ssim=luminance*contrast*structure;
    return ssim;
}