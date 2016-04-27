#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <cstdio>
#include <string>
#include <string.h>
#include <iostream>
#include "OSUFlow.h"
#include "system/cmd_arg_reader.h"
#include "system/path.h"
#include "macros.h"
#include "file/nrrd.h"
#include "statistics.h"
#include "cp_time.h"
#include "thread/ThreadClass.h"
#include "omp.h"
#include "online_quad_bezier_fit.h"
using namespace std;
using namespace JCLib; // https://github.com/chunmingchen/JCLib

vector<string> fileAry;
int W, H, D, Ts;
#define POS_ID(x,y,z) ((x)+W*((y)+H*(z)));
bool saveErrHist=false;


inline VECTOR3 &operator +=(VECTOR3 & v0, const VECTOR3 & v1)
{
    v0[0]+=v1[0]; v0[1]+=v1[1]; v0[2]+=v1[2];
    return v0;
}

inline VECTOR3 operator*(const VECTOR3 & v0, const VECTOR3 & v1)
{
    VECTOR3 y(v0[0]*v1[0], v0[1]*v1[1], v0[2]*v1[2]); return y;
}
inline VECTOR3 sqrt(const VECTOR3 &x)
{
    VECTOR3 y(sqrt(x[0]), sqrt(x[1]), sqrt(x[2]));
}
template<typename T> inline T &getElem(T &x, int d) {return x;}
inline float & getElem(VECTOR3 &v, int d) {return v[d];}

////////////////////////////////////////////////////////


template<class T, int DIMS>
class ErrorModeling{

public:
    void test_online() {
        //float seq[] = {1, 2, 3, 100, 5, 6, 7, 8, 9, 10};
        float seq[] = {1, 8, 27, 100, 125, 196, 343, 512, 729, 1111};
        //float seq[] = {1, 4, 9, 16, 25, 36, 49, 64, 81, 100};
        int n = 10;
        int sampling = n-1;
        int i;
        vector<OnlineQuadBezierFit<T> > onlineAry(1);
        W=H=D=1;

        //////////////////////////////////////
        /// test online version
        //////////////////////////////////////
        for (i=0; i<n; i++)
        {
            onlineAry[0].addData(seq[i], (double)i/(n-1));
        }

        fitOnlineQuadBezier(quadBezierAry, rmsErrAry, onlineAry, n);
        cout << "online:  Ctrl=" <<  quadBezierAry[0] << ", stderr=" << rmsErrAry[0] << endl;

        //////////////////////////////////////
        /// test offline version
        //////////////////////////////////////
        vector<vector<T> > fieldAry(n, vector<T>(1));
        for (i=0; i<n; i++)
        {
            fieldAry[i][0] = seq[i];
        }
        fitQuadBezier(fieldAry, 0, n);
        cout <<"offline: Ctrl=" << quadBezierAry[0] << ", stderr=" << rmsErrAry[0] << endl;

        // print out
        for (i=0; i<n; i++)
        {
            T y0 = onlineAry[0].y0;
            T yn = onlineAry[0].y1;
            T ctrl = quadBezierAry[0];
            float t = (float)i/(n-1);
            float u = (float)(n-1-i)/(n-1);
            T fitted = y0*(u*u) + ctrl*(2*t*u) + yn*(t*t);
            cout << fitted << endl;
        }
    }


    void reportTime() {
        printf("sampling=%d, threads=%d, Fit times: %d\n", GET_ARG_INT("sampling"), GET_ARG_INT("threads"), fitTimes);
        printf("Average online add time, fit time (ms): %.5lf, %.5lf\n", totalAddTime.getValue()*1e-3/addTimes, totalFitTime.getValue()*1e-3/fitTimes );
    }

private:  // only for offline:

    // online timing
    AtomicLong totalAddTime, totalFitTime;
    int addTimes, fitTimes;

    vector<vector<T> > quadFuncAry; // xyz, abc, dim
    vector<T> rmsErrAry;  // mean square error;  xyz, dim
    vector<T> stdErrAry;  // std error;  xyz, dim
    vector<T> meanErrAry;  // std error;  xyz, dim
    vector<vector<T> > errAry;  // std error;  xyz, dim

    // fitting with quad bezier
    vector<T> quadBezierAry; // control point
    vector<T> y0Ary, ynAry;

    // flowfieldAry: one layer of the original vector data
    // z: the z of the layer
    // start: file id
    void fitQuadBezier(vector<vector<T> > &flowfieldAry, int z, int n)
    {
        if (z==0) {
            quadBezierAry = vector<T> (W*H*D);
            stdErrAry =  vector<T> (W*H*D);
            rmsErrAry =  vector<T> (W*H*D);
            //errAry = vector<vector<T> >(w*h*d, vector<T>(n));
            meanErrAry = vector<T> (W*H*D);
            y0Ary = vector<T> (W*H*D);
            ynAry = vector<T> (W*H*D);
        }

        println ("Fitting quadratic function");
        int x, y, i;
        //for (z = 0; z < d; z++)
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                {
                    int id_in = POS_ID(x,y,0);
                    int id_out = POS_ID(x,y,z);

                    vector<T> y_ary(n);
                    for (i=0; i<n; i++)
                        y_ary[i] = flowfieldAry[i][id_in];

                    for (int d=0; d<DIMS; d++) //dim
                    {
                        float y1=getElem(y_ary[0], d);
                        float yn=getElem(y_ary[n-1], d);
                        getElem(y0Ary[id_out], d) = y1;
                        getElem(ynAry[id_out], d) = yn;
                        // t: x,  u: 1-x
                        double sum_t1u1y=0, sum_t1u3=0, sum_t3u1=0, sum_t2u2=0;
                        //println("y=");
                        for (i=0; i<n; i++)
                        {
                            double t = (double)i/(n-1);
                            double u = 1-t;
                            double u2 = u*u;
                            double t2 = t*t;
                            sum_t1u1y 	+= getElem(y_ary[i], d)*t*u;
                            sum_t1u3 	+= t*u*u2;
                            sum_t2u2	+= t2*u2;
                            sum_t3u1	+= t*t2*u;
                        }
                        float ctrl = (sum_t1u1y - sum_t1u3 * y1 - sum_t3u1 * yn) / 2. / sum_t2u2 ;
                        getElem(quadBezierAry[id_out], d) = ctrl;

                        //println("ctrl=%f,   t1u1y=%lf, t1u3=%lf, t2u2=%lf, t3u1=%lf", ctrl, sum_t1u1y, sum_t1u3, sum_t2u2, sum_t3u1);

                        //printf("d=%d, asserting...%f\n", d, abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-yn));
                        //assert(abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-y_ary[n-1][d])<1e-5);

                        // gen err std
                        //println("Val, Err:");
                        vector<float> err_ary(n), err_ary2(n);
                        for (i=0; i<n; i++)
                        {
                            float truth = getElem(y_ary[i], d);
                            float t = (float)i/(n-1);
                            float u = (float)(n-1-i)/(n-1);
                            float fitted = u*u*y1 + 2*t*u*ctrl + t*t*yn;
                            err_ary[i] = truth-fitted;
                            err_ary2[i] = err_ary[i]*err_ary[i];
                            //println("%f %f", fitted, err_ary[i]);
                        }
                        double mean = JCLib::getMean(err_ary.begin(), err_ary.end());
                        double std = JCLib::getDeviation(err_ary.begin(), err_ary.end());
                        double rms = sqrt(JCLib::getSum(err_ary2.begin(), err_ary2.end())/(n-1));
                        getElem(stdErrAry[id_out],d) = std;
                        getElem(rmsErrAry[id_out],d) = rms;
                        getElem(meanErrAry[id_out],d) = mean;
                        //println("std=%lf mean=%lf rms=%f", std, mean, rms);

                        // errAry
                        //for (i=0; i<n; i++)
                        //	errAry[id_out][i][d]=err_ary[i];

                    }


                    //VECTOR3 err_std;
                    //VECTOR3 sum_err(0,0,0);
                    //for (i=0; i<n; i++)
                    //	sum_err = sum_err+(par[0]*i*i + par[1]*i + par[2] - y_ary[i]);

                    //getchar();

                }
    }



};

int main(int argc, const char **argv) {
    ErrorModeling<float, 1>em;
    em.test_online();
    return 0;

}
