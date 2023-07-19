#include <stdio.h>
#include <iostream>
#include "device_launch_parameters.h"
#include <cuda_runtime.h>
#define pi 3.141592653589793
#define CUDA_CHECK(err)                                                                            \
    do {                                                                                           \
        cudaError_t err_ = (err);                                                                  \
        if (err_ != cudaSuccess) {                                                                 \
            std::printf("CUDA error %d at %s:%d\n", err_, __FILE__, __LINE__);                     \
        }                                                                                          \
    } while (0)
using namespace std;

__device__ void full_model(double x, double y, double z, double u, double R, double* C, double* S, double* fac, double* f1, double* f2, double* f3, double* J1, double* J2, double* J3, double* J4, double* J5, double* J6) {
    double r = pow(x * x + y * y + z * z, 0.5);
    double miu = z / r;
    double lambda = 2 * pi * (y < 0) + pow(-1, (y < 0)) * acos(x / pow(x * x + y * y, 0.5));
    double P[41][41] = { 0 };//{1,0,miu,pow(1-miu*miu,0.5)};
    double nor, U_, dUdr_, dUdmiu_, dUdlambda_, d2Udr2_, d2Udmiu2_, d2Udlambda2_;
    P[0][0] = 1; P[1][0] = miu; P[1][1] = pow(1 - miu * miu, 0.5);
    double dP[40][41] = { 0 }, d2P[39][41] = { 0 };
    double dUdr = 0, dUdmiu = 0, dUdlambda = 0, U = 0, d2Udr2 = 0, d2Udmiu2 = 0, d2Udlambda2 = 0;
    for (int n = 2; n <= 40; n++) {
        for (int m = 0; m <= n; m++) {
            if (m == 0) {
                P[n][m] = ((2 * n - 1) * miu * P[n - 1][0] - (n - 1) * P[n - 2][0]) / n;
            }
            if (m == n) {
                P[n][m] = (2 * n - 1) * pow(1 - miu * miu, 0.5) * P[n - 1][n - 1];
            }
            if (m == (n - 1)) {
                P[n][m] = (2 * n - 1) * miu * P[n - 1][n - 1];
            }
            if (m != 0 && m != n && m != (n - 1)) {
                P[n][m] = (2 * n - 1) * miu * P[n - 1][m] / (n - m) - (n + m - 1) * P[n - 2][m] / (n - m);
            }
            if (n > 2) {
                dP[n - 1][m] = ((m - n) * P[n][m] + n * miu * P[n - 1][m]) / (1 - miu * miu);
            }
            if (n > 3) {
                d2P[n - 2][m] = ((-1 + n + 2 * miu * miu + 3 * (n - 2) * miu * miu + (n - 2) * (n - 2) * miu * miu) * P[n - 2][m] + (1 + m - n) * ((5 + 2 * (n - 2)) * miu * P[n - 1][m] + (m - n) * P[n][m])) / pow(miu * miu - 1, 2);
                if (m <= n - 2) {
                    nor = pow(((m == 0) * 1 + (m != 0) * 2) * (2 * n - 3) * fac[n - m - 2] / fac[n + m - 2], 0.5);
                    U_ = pow(R / r, n - 2) * nor * P[n - 2][m] * (C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) + S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda));
                    U = U + U_;

                    dUdr_ = pow(R / r, n - 2) * (n - 1) * nor * P[n - 2][m] * (C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) + S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda));
                    dUdmiu_ = pow(R / r, n - 2) * nor * dP[n - 2][m] * (C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) + S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda));
                    dUdlambda_ = pow(R / r, n - 2) * nor * P[n - 2][m] * (-C[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda) * m + S[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) * m);

                    dUdr = dUdr + dUdr_; dUdmiu = dUdmiu + dUdmiu_; dUdlambda = dUdlambda + dUdlambda_;
                    
                    d2Udr2_ = pow(R / r, n - 2) * (n - 1) * n * nor * P[n - 2][m] * (C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) + S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda));
                    d2Udmiu2_ = pow(R / r, n - 2) * nor * d2P[n - 2][m] * (C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) + S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda));
                    d2Udlambda2_ = pow(R / r, n - 2) * nor * P[n - 2][m] * (-C[(n + 3) * (n - 2) / 2 + m] * cos(m * lambda) * m * m - S[(n + 3) * (n - 2) / 2 + m] * sin(m * lambda) * m * m);

                    d2Udr2 = d2Udr2 + d2Udr2_; d2Udmiu2 = d2Udmiu2 + d2Udmiu2_; d2Udlambda2 = d2Udlambda2 + d2Udlambda2_;
                    

                }
            }
        }
    }
    U = (u / r) * (1 + U);

    dUdr = (-u / pow(r, 2)) * (dUdr + 1); dUdmiu = (u / r) * (dUdmiu); dUdlambda = (u / r) * dUdlambda;
    *f1 = dUdr * (x / r) + dUdmiu * (-x * z / pow(r, 3)) + dUdlambda * (-y / (x * x + y * y));
    *f2 = dUdr * (y / r) + dUdmiu * (-y * z / pow(r, 3)) + dUdlambda * (x / (x * x + y * y));
    *f3 = dUdr * (z / r) + dUdmiu * ((x * x + y * y) / pow(r, 3)) + 0;
    
    d2Udr2 = (u / pow(r, 3)) * (2 + d2Udr2); d2Udmiu2 = (u / r) * d2Udmiu2;
    d2Udlambda2 = (u / r) * d2Udlambda2;

    *J1 = d2Udr2 * (pow(x, 2) / pow(r, 2)) + dUdr * (-pow(x, 2) / pow(r, 3) + 1 / r) + d2Udmiu2 * pow(-x * z / pow(r, 3), 2) + dUdmiu * (3 * pow(x, 2) * z / pow(r, 5) - z / pow(r, 3)) + d2Udlambda2 * (pow(-y / (pow(x, 2) + pow(y, 2)), 2)) + dUdlambda * (2 * x * y / pow(pow(x, 2) + pow(y, 2), 2));
    *J2 = d2Udr2 * (x * y / pow(r, 2)) + dUdr * (-x * y / pow(r, 3)) + d2Udmiu2 * (-x * z / pow(r, 3) * (-y * z / pow(r, 3))) + dUdmiu * (3 * x * y * z / pow(r, 5)) + d2Udlambda2 * (-y / (pow(x, 2) + pow(y, 2)) * x / (pow(x, 2) + pow(y, 2))) + dUdlambda * ((pow(y, 2) - pow(x, 2)) / pow(pow(x, 2) + pow(y, 2), 2));
    *J3 = d2Udr2 * (x * z / pow(r, 2)) + dUdr * (-x * z / pow(r, 3)) + d2Udmiu2 * (-x * z / pow(r, 3) * (-pow(z, 2) / pow(r, 3) + 1 / r)) + dUdmiu * (3 * x * pow(z, 2) / pow(r, 5) - x / pow(r, 3)) + 0 + 0;
    *J4 = d2Udr2 * (pow(y, 2) / pow(r, 2)) + dUdr * (-pow(y, 2) / pow(r, 3) + 1 / r) + d2Udmiu2 * (-y * z / pow(r, 3) * (-y * z / pow(r, 3))) + dUdmiu * (3 * pow(y, 2) * z / pow(r, 5) - z / pow(r, 3)) + d2Udlambda2 * (pow(x / (pow(x, 2) + pow(y, 2)), 2)) + dUdlambda * (-2 * y * x / pow(pow(x, 2) + pow(y, 2), 2));
    *J5 = d2Udr2 * (y * z / pow(r, 2)) + dUdr * (-y * z / pow(r, 3)) + d2Udmiu2 * (-y * z / pow(r, 3) * (-pow(z, 2) / pow(r, 3) + 1 / r)) + dUdmiu * (3 * y * pow(z, 2) / pow(r, 5) - y / pow(r, 3)) + 0 + 0;
    *J6 = d2Udr2 * (pow(z, 2) / pow(r, 2)) + dUdr * (-pow(z, 2) / pow(r, 3) + 1 / r) + d2Udmiu2 * (pow(-pow(z, 2) / pow(r, 3) + 1 / r, 2)) + dUdmiu * (3 * pow(z, 3) / pow(r, 5) - 3 * z / pow(r, 3)) + 0 + 0;
  
}

__global__ void dCreatetm(double *d_tm, int m) {
    int i = threadIdx.x;
	d_tm[i]=-cos(double(i)/(m-1)*pi);
}

__global__ void dCreatebasfunq(double*d_q, double* d_tm) {

    int i = threadIdx.x;
    int j = threadIdx.y;
    int jnum = blockDim.y;
    d_q[i * jnum + j] = cos(j * acos(d_tm[i]));
    
}

__global__ void dCreatebasfunqex(double* d_q, double* d_tm) {
 
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y,inum=blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * (inum * bjnum * jnum) + bj * jnum + i *bjnum* jnum + j;
    int indexi = bi * inum + i;
    int indexj = bj * jnum + j;
    d_q[index] = cos(indexj * acos(d_tm[indexi]));

}

__global__ void dCreatebasfunfqv0(double* d_fq, double* d_tm,int m) {//生成行优先存储的fq
    int i = threadIdx.x;
    int j = threadIdx.y;
    int jnum = blockDim.y;
    if (j == 0) {
        d_fq[i * jnum + j] = (d_tm[i] + 1)/2;
    }
    else if (j==1) {
        d_fq[i * jnum + j] = (pow(d_tm[i], 2) - 1) / 2;
    }
    else if (j == m-1) {
        d_fq[i * jnum + j] = (cos((j + 1) * acos(d_tm[i])) - cos((j + 1) * acos(-1.0))) / (4 * j + 4) - (cos((j - 1) * acos(d_tm[i])) - cos((j - 1) * acos(-1.0))) / (4 * j - 4);
    }
    else {
        d_fq[i * jnum + j] = (cos((j+1) * acos(d_tm[i])) - cos((j+1) * acos(-1.0))) / (2 * j+2) - (cos((j - 1) * acos(d_tm[i])) - cos((j - 1) * acos(-1.0))) / (2 * j - 2);
    }
    
}

__global__ void dCreatebasfunfq(double* d_fq, double* d_tm, int m) {//生成列优先存储的fq
    int i = threadIdx.x;
    int j = threadIdx.y;
    int inum = blockDim.x;
    if (j == 0) {
        d_fq[j * inum + i] = (d_tm[i] + 1) / 2;
    }
    else if (j == 1) {
        d_fq[j * inum + i] = (pow(d_tm[i], 2) - 1) / 2;
    }
    else if (j == m - 1) {
        d_fq[j * inum + i] = (cos((j + 1) * acos(d_tm[i])) - cos((j + 1) * acos(-1.0))) / (4 * j + 4) - (cos((j - 1) * acos(d_tm[i])) - cos((j - 1) * acos(-1.0))) / (4 * j - 4);
    }
    else {
        d_fq[j * inum + i] = (cos((j + 1) * acos(d_tm[i])) - cos((j + 1) * acos(-1.0))) / (2 * j + 2) - (cos((j - 1) * acos(d_tm[i])) - cos((j - 1) * acos(-1.0))) / (2 * j - 2);
    }

}

__global__ void dCreatebasfunfqexv0(double* d_fq, double* d_tm, int m) {
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * (inum * bjnum * jnum) + bj * jnum + i * bjnum * jnum + j;
    int indexi = bi * inum + i;
    int indexj = bj * jnum + j;
    if (indexj == 0) {
        d_fq[index] = (d_tm[indexi] + 1) / 2;
    }
    else if (indexj == 1) {
        d_fq[index] = (pow(d_tm[indexi], 2) - 1) / 2;
    }
    else if (indexj == m - 1) {
        d_fq[index] = (cos((indexj + 1) * acos(d_tm[indexi])) - cos((indexj + 1) * acos(-1.0))) / (4 * indexj + 4) - (cos((indexj - 1) * acos(d_tm[indexi])) - cos((indexj - 1) * acos(-1.0))) / (4 * indexj - 4);
    }
    else {
        d_fq[index] = (cos((indexj + 1) * acos(d_tm[indexi])) - cos((indexj + 1) * acos(-1.0))) / (2 * indexj + 2) - (cos((indexj - 1) * acos(d_tm[indexi])) - cos((indexj - 1) * acos(-1.0))) / (2 * indexj - 2);
    }
}

__global__ void dCreatebasfunfqex(double* d_fq, double* d_tm, int m) {
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int indexi = bi * inum + i;
    int indexj = bj * jnum + j;
    int index = indexi + indexj * binum * inum;
    if (indexj == 0) {
        d_fq[index] = (d_tm[indexi] + 1) / 2;
    }
    else if (indexj == 1) {
        d_fq[index] = (pow(d_tm[indexi], 2) - 1) / 2;
    }
    else if (indexj == m - 1) {
        d_fq[index] = (cos((indexj + 1) * acos(d_tm[indexi])) - cos((indexj + 1) * acos(-1.0))) / (4 * indexj + 4) - (cos((indexj - 1) * acos(d_tm[indexi])) - cos((indexj - 1) * acos(-1.0))) / (4 * indexj - 4);
    }
    else {
        d_fq[index] = (cos((indexj + 1) * acos(d_tm[indexi])) - cos((indexj + 1) * acos(-1.0))) / (2 * indexj + 2) - (cos((indexj - 1) * acos(d_tm[indexi])) - cos((indexj - 1) * acos(-1.0))) / (2 * indexj - 2);
    }

}

__global__ void dCreatebasfundqv0(double* d_dq, double* d_tm, int m) {//生成行优先存储的dq
    int i = threadIdx.x;
    int j = threadIdx.y;
    int jnum = blockDim.y;
    double a = (j%(m-1)==0) ? 0.5 : 1;
    if (i == 0) {
        d_dq[i * jnum + j] = pow(-1.0,j+1)*j*j*a;
    }
    else if (i == m-1) {
        d_dq[i * jnum + j] = j*j*a;
    }
    else {
        d_dq[i * jnum + j] = (j * sin(j * acos(d_tm[i]))) / pow((1 - pow(d_tm[i],2)),0.5)*a;
    }
}

__global__ void dCreatebasfundq(double* d_dq, double* d_tm, int m) {//生成列优先存储的dq
    int i = threadIdx.x;
    int j = threadIdx.y;
    int inum = blockDim.x;
    double a = (j % (m - 1) == 0) ? 0.5 : 1;
    if (i == 0) {
        d_dq[j * inum + i] = pow(-1.0, j + 1) * j * j * a;
    }
    else if (i == m - 1) {
        d_dq[j * inum + i] = j * j * a;
    }
    else {
        d_dq[j * inum + i] = (j * sin(j * acos(d_tm[i]))) / pow((1 - pow(d_tm[i], 2)), 0.5) * a;
    }
}

__global__ void dCreatebasfundqexv0(double* d_dq, double* d_tm, int m) {
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * (inum * bjnum * jnum) + bj * jnum + i * bjnum * jnum + j;
    int indexi = bi * inum + i;
    int indexj = bj * jnum + j;
    double a = (indexj % (m - 1) == 0) ? 0.5 : 1;
    if (indexi == 0) {
        d_dq[index] = pow(-1.0, indexj + 1) * indexj * indexj * a;
    }
    else if (indexi == m - 1) {
        d_dq[index] = indexj * indexj * a;
    }
    else {
        d_dq[index] = (indexj * sin(indexj * acos(d_tm[indexi]))) / pow((1 - pow(d_tm[indexi], 2)), 0.5) * a;
    }
}

__global__ void dCreatebasfundqex(double* d_dq, double* d_tm, int m) {
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int indexi = bi * inum + i;
    int indexj = bj * jnum + j;
    int index = indexi + indexj * binum * inum;
    double a = (indexj % (m - 1) == 0) ? 0.5 : 1;
    if (indexi == 0) {
        d_dq[index] = pow(-1.0, indexj + 1) * indexj * indexj * a;
    }
    else if (indexi == m - 1) {
        d_dq[index] = indexj * indexj * a;
    }
    else {
        d_dq[index] = (indexj * sin(indexj * acos(d_tm[indexi]))) / pow((1 - pow(d_tm[indexi], 2)), 0.5) * a;
    }
}

__global__ void dCreateqnr(double* d_qnr, int m) {
    int i = threadIdx.x;
    if (i == 0) {
        d_qnr[i] = 1.0/(m-1);
    }
    else if (i == m - 1) {
        d_qnr[i] = 1.0 / (m - 1);
    }
    else {
        d_qnr[i] = 2.0 / (m - 1);
    }

}

__global__ void dCreatet(double *d_t,double *d_tm,double ti,double dt) {
    int i = threadIdx.x;
	d_t[i] = (d_tm[i]*dt+2*ti+dt)/2;
}

__global__ void dCreateGx(double *d_Gx,double *d_y_oldv,double dt,int m) {
    int i = threadIdx.x;
    d_Gx[i] = d_y_oldv[i]*dt/2;
    d_Gx[m+i] = d_y_oldv[m+i] * dt/2;
    d_Gx[2*m + i] = d_y_oldv[2*m + i] * dt/2;
}

__global__ void dCreateGvA(double* d_Gv, double* d_A, double* d_y_oldx, double dt ,int m,double u, double R, double* d_C, double* d_S, double* d_fac) {
    int i = threadIdx.x;
    double f1 = 0, f2 = 0, f3 = 0;
    double A1 = 0, A2 = 0, A3 = 0, A4 = 0, A5 = 0, A6 = 0;
    full_model(d_y_oldx[i], d_y_oldx[m + i], d_y_oldx[2 * m + i], u, R, d_C, d_S, d_fac, &f1, &f2, &f3, &A1, &A2, &A3, &A4, &A5, &A6);
    d_Gv[i] = f1 * dt / 2;
    d_Gv[m + i] = f2 * dt / 2;
    d_Gv[2 * m + i] = f3 * dt / 2;
    d_A[i] = A1 * dt / 2;
    d_A[m + i] = A2 * dt / 2;
    d_A[2 * m + i] = A3 * dt / 2;
    d_A[3 * m + i] = A4 * dt / 2;
    d_A[4 * m + i] = A5 * dt / 2;
    d_A[5 * m + i] = A6 * dt / 2;
}

__global__ void dCreateA(double* d_A, double* d_y_old,  int m, double dt, double u) {
    int i = threadIdx.x;
    double r = pow((pow(d_y_old[i], 2) + pow(d_y_old[m + i], 2) + pow(d_y_old[2 * m + i], 2)), 2.5)*2/dt;
    d_A[i] = -1 * u * (pow(d_y_old[m + i], 2) + pow(d_y_old[2 * m + i], 2) - 2 * pow(d_y_old[i], 2)) / r;
    d_A[m + i] = 3 * u *d_y_old[i]  * d_y_old[m + i]/ r;
    d_A[2 * m + i] = 3 * u * d_y_old[i] * d_y_old[2*m + i] / r;
    d_A[3 * m + i] = -1 * u * (pow(d_y_old[i], 2) + pow(d_y_old[2 * m + i], 2) - 2 * pow(d_y_old[m+i], 2)) / r;
    d_A[4 * m + i] = 3 * u * d_y_old[m+i] * d_y_old[2 * m + i] / r;
    d_A[5 * m + i] = -1 * u * (pow(d_y_old[i], 2) + pow(d_y_old[m + i], 2) - 2 * pow(d_y_old[2*m + i], 2)) / r;
}

__global__ void dCreateJ(double* d_J, double* d_AA, double* d_Te, double* d_P, int m, double dt) {
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * inum + bj * binum * inum * jnum + i + j * binum * inum;
    int pa = (bi == bj) ? -1 : 0;
    double ta = ((bj - 3) == bi) ? dt / 2 : 0;
    int Aa = (bi > 2) ? 1 : 0;
    int Ab = (bj < 3) ? 1 : 0;
    int Ac = (bi - 3) + bj + ((bj - 3) && bj);
    int Ad = (Ac > 0 && Ac < 6) ? Ac : 0;
    d_J[index] = pa * d_P[i + j * inum] + ta * d_Te[i + j * inum] + Aa * Ab * d_AA[Ad * m * m + i + j * inum];
}

__global__ void dCreateJC(double* d_J, double* d_Te, double* d_P, int m, double dt) {
    int mul = m / 32;
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * inum + bj * binum * inum * jnum + i + j * binum * inum;
    int myi = i + inum * (bi % mul);
    int myj = j + jnum * (bj % mul);
    int mybi = bi / mul;
    int mybj = bj / mul;
    int pa = (mybi == mybj) ? -1 : 0;
    double ta = ((mybj - 3) == mybi) ? dt / 2 : 0;
    d_J[index] = pa * d_P[myi + myj * inum * mul] + ta * d_Te[myi + myj * inum * mul];
}

__global__ void dCreateJex(double* d_J, double* d_AA, double* d_Te, double* d_P, int m, double dt) {//待续
    int mul = m / 32;
    int i = threadIdx.x; int bi = blockIdx.x;
    int j = threadIdx.y; int bj = blockIdx.y;
    int jnum = blockDim.y, inum = blockDim.x; int bjnum = gridDim.y, binum = gridDim.x;
    int index = bi * inum + bj * binum * inum * jnum + i + j * binum * inum;
    int myi = i + inum * (bi % mul);
    int myj = j + jnum * (bj % mul);
    int mybi = bi / mul;
    int mybj = bj / mul;
    int pa = (mybi == mybj) ? -1 : 0;
    double ta = ((mybj - 3) == mybi) ? dt / 2 : 0;
    int Aa = (mybi > 2) ? 1 : 0;
    int Ab = (mybj < 3) ? 1 : 0;
    int Ac = (mybi - 3) + mybj + ((mybj - 3) && mybj);
    int Ad = (Ac > 0 && Ac < 6) ? Ac : 0;
    d_J[index] = pa * d_P[myi + myj * inum*mul] + ta * d_Te[myi + myj * inum*mul] + Aa * Ab * d_AA[Ad * m * m + myi + myj * inum*mul];
}

__global__ void dCreatey0(double* d_y0, double x0,double y0, double z0, double vx0,double vy0,double vz0,int m) {
    int i = threadIdx.x;
    d_y0[i] = x0;
    d_y0[m+i] = y0;
    d_y0[2*m+i] = z0;
    d_y0[3*m+i] = vx0;
    d_y0[4*m+i] = vy0;
    d_y0[5*m+i] = vz0;
}


/////////////////////////////////// 上面的代码是CUDA核函数，下面的是暴露给主机的接口/////////////////////////////////////



extern void hCreatetm(double* d_tm, int m)
{
    cudaError_t cudaStatus;
    dCreatetm << <1, m >> > (d_tm, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatetm launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfunq(double* d_q, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    dim3 block(m, m);
    dCreatebasfunq << <1, block >> > (d_q, d_tm);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfunq launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfunqex(double* d_q, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    int n = (m < 32) ? m : 32;
    dim3 grid(m/32,m/32),block(n, n);
    dCreatebasfunqex << <grid, block >> > (d_q, d_tm);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfunq launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfunfq(double* d_fq, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    dim3 block(m, m);
    dCreatebasfunfq << <1, block >> > (d_fq, d_tm, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfunfq launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfunfqex(double* d_fq, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    int n = (m < 32) ? m : 32;
    dim3 grid(m/32, m/32), block(n, n);
    dCreatebasfunfqex << <grid, block >> > (d_fq, d_tm, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfunfqex launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfundq(double* d_dq, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    dim3 block(m, m);
    dCreatebasfundq << <1, block >> > (d_dq, d_tm, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfundq launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatebasfundqex(double* d_dq, double* d_tm, int m)
{
    cudaError_t cudaStatus;
    int n = (m < 32) ? m : 32;
    dim3 grid(m/32, m/32), block(n, n);
    dCreatebasfundqex << <grid, block >> > (d_dq, d_tm, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatebasfundqex launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreateqnr(double* d_qnr,int m)
{
    cudaError_t cudaStatus;
    dCreateqnr << <1, m >> > (d_qnr, m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateqnr launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatet(double* d_t, double* d_tm, int m, double ti, double dt)
{
    cudaError_t cudaStatus;
    dCreatet << <1, m >> > (d_t, d_tm, ti, dt);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatet launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreateGA(double* d_G, double* d_A, double* d_y_old,int m,double dt,double u, double R, double* d_C, double* d_S, double* d_fac)
{
    cudaError_t cudaStatus;
    dCreateGx << <1, m >> > (d_G, &d_y_old[3*m],dt,m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateGx launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
    dCreateGvA << <1, m >> > (&d_G[3*m], d_A, d_y_old, dt,m,u,R,d_C, d_S, d_fac);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateGvA launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreateA(double* d_A, double* d_y_old, int m, double dt, double u)
{
    cudaError_t cudaStatus;
    dCreateA << <1, m >> > (d_A, d_y_old, m,dt,u);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateA launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}
extern void hCreateJ(double* d_J, double* d_AA, double* d_Te, double* d_P, int m,int n, double dt)
{
    cudaError_t cudaStatus;
    dim3 grid(n, n), block(m,m);
    dCreateJ << <grid, block >> > (d_J,d_AA,d_Te,d_P,m,dt);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateJ launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreateJC(double* d_J, double* d_Te, double* d_P, int m, int n, double dt)
{
    cudaError_t cudaStatus;
    int bm = (m < 32) ? m : 32;
    dim3 grid(n * m / 32, n * m / 32), block(bm, bm);
    dCreateJC << <grid, block >> > (d_J, d_Te, d_P, m, dt);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateJC launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreateJex(double* d_J, double* d_AA, double* d_Te, double* d_P, int m, int n, double dt)
{
    cudaError_t cudaStatus;
    int bm = (m < 32) ? m : 32;
    dim3 grid(n*m/32, n*m/32), block(bm, bm);
    dCreateJex << <grid, block >> > (d_J, d_AA, d_Te, d_P, m, dt);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreateJex launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

extern void hCreatey0(double* d_y0,double x0, double y0, double z0, double vx0, double vy0, double vz0, int m)
{
    cudaError_t cudaStatus;
    dCreatey0 << <1, m >> > (d_y0, x0,y0,z0,vx0,vy0,vz0,m);
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "hCreatey0 launch failed: %s\n", cudaGetErrorString(cudaStatus));
    }
}

