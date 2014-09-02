#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <cstdio>
#include <string>
#include <string.h>
#include "OSUFlow.h"
#include "system/cmd_arg_reader.h"
#include "system/path.h"
#include "macros.h"
#include "file/nrrd.h"
#include "statistics.h"
#include "cp_time.h"
#include "thread/ThreadClass.h"
#include "omp.h"
using namespace std;
using namespace JCLib;

vector<string> fileAry;
int W, H, D, Ts;
#define POS_ID(x,y,z) ((x)+W*((y)+H*(z)));
bool saveErrHist=false;

///////////////////////
struct SINGLETON{
    float v;
    SINGLETON() {v=0;}
    inline float &operator[](int i) { assert(i==0); return v; }
    inline SINGLETON operator+(const SINGLETON &x)  { SINGLETON y; y.v = v+x.v; return y; }
    inline SINGLETON operator-(const SINGLETON &x)  { SINGLETON y; y.v = v-x.v; return y; }
    inline SINGLETON operator*(const SINGLETON &x)  { SINGLETON y; y.v = v*x.v; return y; }
    inline SINGLETON operator*(float f)  { SINGLETON y; y.v = v*f; return y; }
    inline SINGLETON operator*(double f)  { SINGLETON y; y.v = v*f; return y; }
    inline SINGLETON &operator+=(const SINGLETON &x)  { v += x.v; return *this; }
};

inline VECTOR3 &operator +=(VECTOR3 & v0, const VECTOR3 & v1)
{
    v0[0]+=v1[0]; v0[1]+=v1[1]; v0[2]+=v1[2];
    return v0;
}
inline VECTOR3 mult(const VECTOR3 & v0, const VECTOR3 & v1)
{
    VECTOR3 y(v0[0]*v1[0], v0[1]*v1[1], v0[2]*v1[2]); return y;
}
inline SINGLETON mult(const SINGLETON & v0, const SINGLETON & v1)
{
    SINGLETON y; y.v = v0.v * v1.v; return y;
}
inline VECTOR3 sqrt(const VECTOR3 &x)
{
    VECTOR3 y(sqrt(x[0]), sqrt(x[1]), sqrt(x[2]));
}
inline SINGLETON sqrt(const SINGLETON x)
{
    SINGLETON y; y.v = sqrt(x.v); return y;
}


////////////////////////////////////////////////////////

// global: fileAry, w,h,d,t
void load_list(string list_filename) {
	int i, j;
	string path = getPath(list_filename);

	FILE *fp = fopen(list_filename.c_str(), "rt");
	char str[1024];
	fgets(str, 1024, fp);
    i = sscanf(str, "%d %d %d %d", &W, &H, &D, &Ts);
	assert(i==4);
	fgets(str, 1024, fp); // dummy line


	fileAry.clear();
    for (i = 0; i < Ts; i++) {
		fgets(str, 1024, fp);
		char *p = strtok(str, " \r\n");
		//printf("%s\n", p);

		if (isFilenameOnly(p))
			fileAry.push_back(path + p);
		else
			fileAry.push_back(p);
	}
}



template<class T, int DIMS>
class ErrorModeling{

    vector<vector<T> > quadFuncAry; // xyz, abc, dim
    vector<T> rmsErrAry;  // mean square error;  xyz, dim
    vector<T> stdErrAry;  // std error;  xyz, dim
    vector<T> meanErrAry;  // std error;  xyz, dim
    vector<vector<T> > errAry;  // std error;  xyz, dim

    // online timing
    AtomicLong totalAddTime, totalFitTime;
    int addTimes, fitTimes;

    // fitting with quad with connected end points
    void fitQuadFuncEndPoints(vector<vector<T> > &flowfieldAry, int start, int sampling)
    {
        int n=sampling+1;
        quadFuncAry = vector<vector<T> > (W*H*D, vector<T>(3));
        stdErrAry =  vector<T> (W*H*D);
        rmsErrAry =  vector<T> (W*H*D);
        errAry = vector<vector<T> >(W*H*D, vector<T>(n));
        meanErrAry = vector<T> (W*H*D);

        println ("Fitting quadratic function");
        int x, y, z, i;
        for (z = 0; z < D; z++)
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                {
                    int id = POS_ID(x,y,z);
                    vector<T> y_ary(n);
                    for (i=0; i<n; i++)
                        y_ary[i] = flowfieldAry[i][id];

                    for (int d=0; d<DIMS; d++) //dim
                    {
                        float y1=y_ary[0][d];
                        float yn=y_ary[n-1][d];
                        float sum_y=0, sum_y2=0, sum_x2y=0, sum_xy=0;
                        //println("y=");
                        for (i=0; i<n; i++)
                        {
                            sum_y += y_ary[i][d];
                            sum_y2+= pow(y_ary[i][d], 2);
                            sum_xy += i* y_ary[i][d];
                            sum_x2y += i*i * y_ary[i][d];
                            //println("%f", y_ary[i][d]);
                        }
                        float sum_x = (n-1)*n/2;
                        float sum_x2 = (n-1)*(n)*(2*n-1)/6.f;
                        float sum_x3 = sum_x*sum_x;
                        float sum_x4 = (n-1)*n*(2*n-1)*(3*(n-1)*(n-1)+3*(n-1)-1)/30.f;

                        float par[3];
                        par[0] = (sum_x3*(yn-y1)/(n-1) - sum_x2y - yn*sum_x2 + 2*y1*sum_x2 + (n-1)*sum_xy - (n-1)*y1*sum_x)
                                / (2*sum_x3*(n-1)-sum_x4-sum_x2*(n-1)*(n-1));
                        par[1] = (	(yn-y1)-par[0]*(n-1)*(n-1))/(n-1);
                        par[2] = y_ary[0][d];
                        quadFuncAry[id][0][d] = par[0];
                        quadFuncAry[id][1][d] = par[1];
                        quadFuncAry[id][2][d] = par[2];

                        //println("a=%f  b=%f  c=%f", par[0], par[1], par[2]);

                        //printf("d=%d, asserting...%f\n", d, abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-yn));
                        //assert(abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-y_ary[n-1][d])<1e-5);

                        // gen err std
                        //println("Err:");
                        vector<float> err_ary(n), err_ary2(n);
                        for (i=0; i<n; i++)
                        {
                            float truth = y_ary[i][d];
                            float fitted = par[0]*(i*i) + par[1]*i + par[2];
                            err_ary[i] = truth-fitted;
                            err_ary2[i] = err_ary[i]*err_ary[i];
                            //println("%f", err_ary[i]);
                        }
                        double mean = JCLib::getMean(err_ary.begin(), err_ary.end());
                        double std = JCLib::getDeviation(err_ary.begin(), err_ary.end());
                        double rms = sqrt(JCLib::getSum(err_ary2.begin(), err_ary2.end())/(n-1));
                        stdErrAry[id][d] = std;
                        rmsErrAry[id][d] = rms;
                        meanErrAry[id][d] = mean;
                        //println("std=%lf", std);

                        // errAry
                        for (i=0; i<n; i++)
                            errAry[id][i][d]=err_ary[i];

                    }


                    //VECTOR3 err_std;
                    //VECTOR3 sum_err(0,0,0);
                    //for (i=0; i<n; i++)
                    //	sum_err = sum_err+(par[0]*i*i + par[1]*i + par[2] - y_ary[i]);

                    //getchar();

                }
    }

    void fitQuadFunc(vector<vector<T> > &flowfieldAry, int start, int sampling)
    {
        int n=sampling+1;
        quadFuncAry = vector<vector<T> > (W*H*D, vector<T>(3));
        stdErrAry =  vector<T>  (W*H*D);
        rmsErrAry = vector<T> (W*H*D);

        println ("Fitting quadratic function");
        int x, y, z, i;
        for (z = 0; z < D; z++)
            for (y = 0; y < H; y++)
                for (x = 0; x < W; x++)
                {
                    int id = POS_ID(x,y,z);
                    vector<T> y_ary(n);
                    for (i=0; i<n; i++)
                        y_ary[i] = flowfieldAry[i][id];

                    vector<VECTOR3> &par = quadFuncAry[id];
                    for (int d=0; d<DIMS; d++) //dim
                    {
                        float y1=y_ary[0][d];
                        float yn=y_ary[n-1][d];
                        float sum_y=0, sum_x2y=0, sum_xy=0;
                        //println("y=");
                        for (i=0; i<n; i++)
                        {
                            sum_y += y_ary[i][d];
                            sum_xy += i* y_ary[i][d];
                            sum_x2y += i*i * y_ary[i][d];
                            //println("%f", y_ary[i][d]);
                        }
                        float sum_x = (n-1)*n/2;
                        float sum_x2 = (n-1)*(n)*(2*n-1)/6.f;
                        float sum_x3 = sum_x*sum_x;
                        float sum_x4 = (n-1)*n*(2*n-1)*(3*(n-1)*(n-1)+3*(n-1)-1)/30.f;

                        MATRIX3 A(VECTOR3(sum_x4, sum_x3, sum_x2), VECTOR3(sum_x3, sum_x2, sum_x), VECTOR3(sum_x2, sum_x, n));
                        VECTOR3 b(sum_x2y, sum_xy, sum_y);
                        MATRIX3 invA;
                        A.inverse(invA); // invA = inv(A)
                        VECTOR3 par_ = invA*b;
                        par[0][d] = par_[0];
                        par[1][d] = par_[1];
                        par[2][d] = par_[2];

                        //println("a=%f  b=%f  c=%f", par[0][d], par[1][d], par[2][d]);

                    }


                    //getchar();

                }
    }


    // fitting with quad bezier
    vector<T> quadBezierAry; // control point
    vector<T> y0Ary, ynAry;


    // flowfieldAry: one layer of the original vector data
    // z: the z of the layer
    // start: file id
    void fitQuadBezier(vector<vector<T> > &flowfieldAry, int z, int sampling)
    {
        int n=sampling+1;
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
                        float y1=y_ary[0][d];
                        float yn=y_ary[n-1][d];
                        y0Ary[id_out][d] = y1;
                        ynAry[id_out][d] = yn;
                        // t: x,  u: 1-x
                        double sum_t1u1y=0, sum_t1u3=0, sum_t3u1=0, sum_t2u2=0;
                        //println("y=");
                        for (int t=0; t<n; t++)
                        {
                            int u = n-1-t;
                            int u2 = u*u;
                            int t2 = t*t;
                            sum_t1u1y 	+= y_ary[t][d]*t*u;
                            sum_t1u3 	+= t*u*u2;
                            sum_t2u2	+= t2*u2;
                            sum_t3u1	+= t*t2*u;
                            //println("%f", y_ary[t][d]);
                        }
                        float ctrl = (sum_t1u1y  * (n-1) * (n-1) - sum_t1u3 * y1 - sum_t3u1 * yn) / 2.f / sum_t2u2 ;
                        quadBezierAry[id_out][d] = ctrl;

                        //println("ctrl=%f,   t1u1y=%lf, t1u3=%lf, t2u2=%lf, t3u1=%lf", ctrl, sum_t1u1y, sum_t1u3, sum_t2u2, sum_t3u1);

                        //printf("d=%d, asserting...%f\n", d, abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-yn));
                        //assert(abs(par[0][d]*(n-1)*(n-1)+par[1][d]*(n-1)+par[2][d]-y_ary[n-1][d])<1e-5);

                        // gen err std
                        //println("Val, Err:");
                        vector<float> err_ary(n), err_ary2(n);
                        for (i=0; i<n; i++)
                        {
                            float truth = y_ary[i][d];
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
                        stdErrAry[id_out][d] = std;
                        rmsErrAry[id_out][d] = rms;
                        meanErrAry[id_out][d] = mean;
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


    void saveFittedFlowfields(vector<vector<T> >& funcAry, int start, int sampling, const char *folder)
    {
        println("Saving fitted flow fields...");
        int i;
        FILE *fp;
#if 0 // store each timestep
        int skip = 1;
        vector<T> flowfield(W*H*D);
        for (i=0; i<=sampling; i+=skip)
        {
            // gen field
#if 0 // normal polynomial
            for (int id=0; id<W*H*D; id++)
            {
                flowfield[id] = funcAry[id][0]*(i*i) + funcAry[id][1]*i + funcAry[id][2];
            }
#else // bezier
            // gen err std
            //println("Val, Err:");
            int n = sampling+1;
            for (int id=0; id<W*H*D; id++)
            {
                T y0 = y0Ary[id];
                T yn = ynAry[id];
                T ctrl = quadBezierAry[id];
                float t = (float)i/(n-1);
                float u = (float)(n-1-i)/(n-1);
                T fitted = y0*(u*u) + ctrl*(2*t*u) + yn*(t*t);
                flowfield[id] = fitted;
            }
#endif

            // save file
            if (DIMS==3) {
                fp = fopen(strprintf("%s/sampling%d_%02d.vec", folder, sampling, start+i).c_str(), "wb");
                {
                    int dim[3];
                    dim[0] = W; dim[1] = H; dim[2] = D;
                    fwrite(dim, 4, 3, fp);
                }
            } else
                fp = fopen(strprintf("%s/sampling%d_%02d.raw", folder, sampling, start+i).c_str(), "wb");
            fwrite(&flowfield[0], sizeof(T), W*H*D, fp);
            fclose(fp);
        }
#endif

        // save quad bezier control point
        fp = fopen(strprintf("%s/sampling%d_%02d_bctrl.raw", folder, sampling, start).c_str(), "wb");
        fwrite(&quadBezierAry[0], sizeof(T), W*H*D, fp);
        fclose(fp);
    #if 1
        // save rms err
        fp = fopen(strprintf("%s/sampling%d_%02d_rmserr.raw", folder, sampling, start).c_str(), "wb");
        fwrite(&rmsErrAry[0], sizeof(T), W*H*D, fp);
        fclose(fp);
        // save errstd
        //fp = fopen(strprintf("%s/sampling%d_%02d_stderr.raw", folder, sampling, start).c_str(), "wb");
        //fwrite(&stdErrAry[0], 4, W*H*D*3, fp);
        //fclose(fp);
        // save err mean
        //fp = fopen(strprintf("%s/sampling%d_%02d_meanerr.raw", folder, sampling, start).c_str(), "wb");
        //fwrite(&meanErrAry[0], 4, W*H*D*3, fp);
        //fclose(fp);
    #endif
    #if 0
        // save a
        fp = fopen(strprintf("%s/sampling%d_%02d_a.raw", folder, sampling, start).c_str(), "wb");
        for (int id=0; id<w*h*d; id++)
            fwrite(&funcAry[id][0], sizeof(T), 1, fp);
        fclose(fp);
    #endif
    #if 0
        // errAry
        fp = fopen(strprintf("%s/sampling%d_%02d_err.raw", folder, sampling, start).c_str(), "wb");
        for (int id=0; id<w*h*d; id++)
            fwrite(&errAry[id][0], sizeof(T), (sampling+1), fp);
        fclose(fp);
    #endif
    }

public:
    void run_offline(int sampling) {
        int i, j, z;
        println("allocating size=%d", sampling);
        vector<vector<T> > flowfieldAry(sampling+1, vector<T>(W * H * 1));
        string out_path = GET_ARG_STRING("out_path").c_str();

        int count = 0;
    #if 1
        for (i = 0; i < Ts; i += sampling)
        {
            vector<FILE *> fp_ary;
            if (i+sampling >= Ts) {
                break;
            }
            for (j=0; j <= sampling; j++)
            {
                int id = i + j;
                if (id >= Ts) {
                    break;
                }

                println("Opening %s", fileAry[id].c_str());
                FILE *fp = fopen(fileAry[id].c_str(), "rb");
                if (!fp) {
                    perror("load flow field");
                    exit(1);
                }
                if (DIMS==3) { // for .vec format
                    fseek(fp, 12, SEEK_SET);
                }
                fp_ary.push_back(fp);
            }

            for (z = 0; z < D; z++ )
            {
                println("Reading files for z=%d", z);
                for (j = 0; j <= sampling; j++) {
                    //println("Reading %s for z=%d", fileAry[id].c_str(), z);
                    fread(&flowfieldAry[j][0], sizeof(T), W * H * 1, fp_ary[j]);
                    //print_array(&flowfieldAry[j][0][0], 3);
                }

                //computeMaxError(flowfieldAry, i, sampling);
                //computeMaxAngleMagError(flowfieldAry, i, sampling);
                //genConvex4Error(flowfieldAry, i, sampling);

                //fitQuadFunc(flowfieldAry, i, sampling);
                //saveFittedFlowfields(quadFuncAry, i, sampling, "fitted_quad");

                //fitQuadFuncEndPoints(flowfieldAry, i, sampling);
                //saveFittedFlowfields(quadFuncAry, i, sampling, out_path.c_str());

                fitQuadBezier(flowfieldAry, z, sampling);


            }
            for (j = 0; j <= sampling ; j++) {
                fclose(fp_ary[j]);
            }

            if (i==0) {
                std::string cmd ;
                if (DIMS==3)
                    cmd = strprintf("cp %s %s/sampling%d_%02d.vec", fileAry[i].c_str(), out_path.c_str(), sampling, i);
                else
                    cmd = strprintf("cp %s %s/sampling%d_%02d.raw", fileAry[i].c_str(), out_path.c_str(), sampling, i);
                println("%s", cmd.c_str());
                system(cmd.c_str());

            }
            {
                std::string cmd ;
                if (DIMS==3)
                    cmd = strprintf("cp %s %s/sampling%d_%02d.vec", fileAry[i+sampling].c_str(), out_path.c_str(), sampling, i+sampling);
                else
                    cmd = strprintf("cp %s %s/sampling%d_%02d.raw", fileAry[i+sampling].c_str(), out_path.c_str(), sampling, i+sampling);
                println("%s", cmd.c_str());
                system(cmd.c_str());
            }
            saveFittedFlowfields(quadFuncAry, i, sampling, out_path.c_str());

            count ++;
        }
    #else
        for (i = 0; i < Ts; i += sampling)
        {
            if (i+sampling >= Ts) {
                break;
            }
            count ++;
        }

    #endif

        FILE *fp = fopen(strprintf("%s/all_bezier_rms.list", out_path.c_str()).c_str(), "w");
        fprintf(fp, "%d %d %d %d\n", W, H, D, count);
        fprintf(fp, "%lg\n", 1.f/(double)sampling);
        for (i=0; i< count ; i++)
        {
            fprintf(fp, "sampling%d_%02d.vec sampling%d_%02d_bctrl.raw sampling%d_%02d_rmserr.raw\n",
                sampling, sampling*i, sampling, sampling*i, sampling, sampling*i);

            int i1 = i+1;
            if (i==count-1)
                fprintf(fp, "sampling%d_%02d.vec sampling%d_%02d_bctrl.raw sampling%d_%02d_rmserr.raw\n",
                    sampling, sampling*i1, sampling, sampling*i, sampling, sampling*i);

        }
        fclose(fp);
    }


    struct OnlineCache{
        vector<T> ysum;
        vector<T> ytsum;
        vector<T> yt2sum;
        vector<T> y2sum;
        vector<T> y0;
        vector<T> y1;
    };

    void initOnlineQuadBezier(OnlineCache &online, int size) {
        online.ysum = vector<T>(size);    memset(&online.ysum[0], 0, size*sizeof(T));
        online.ytsum = vector<T>(size);   memset(&online.ytsum[0], 0, size*sizeof(T));
        online.yt2sum = vector<T>(size);  memset(&online.yt2sum[0], 0, size*sizeof(T));
        online.y2sum = vector<T>(size);   memset(&online.y2sum[0], 0, size*sizeof(T));
        online.y0 = vector<T>(size);   memset(&online.y0[0], 0, size*sizeof(T));
        online.y1 = vector<T>(size);   memset(&online.y1[0], 0, size*sizeof(T));
    }

    void addOnlineQuadBezier(OnlineCache &online, vector<T> flowfield, float t)
    {
        assert(flowfield.size() == online.ysum.size());
        Timer timer;
        timer.start();

#pragma omp parallel for
        for (int i=0; i<online.ysum.size(); i++)
        {
            T y = flowfield[i];
            T yt = y*t;
            online.ysum[i] += y;
            online.ytsum[i] += yt;
            online.yt2sum[i] += yt*t;
            online.y2sum[i] += mult(y,y); // element-wise multiplication for VECTOR3
        }
        timer.end();

        if (t==0) {
            memcpy(&online.y0[0], &flowfield[0], flowfield.size()*sizeof(T));
        }
        if (t==1.f) {
            memcpy(&online.y1[0], &flowfield[0], flowfield.size()*sizeof(T));
        }

        totalAddTime+=timer.getElapsedUS();
        addTimes++;
    }

    // output : ctrlAry, stderrAry
    void fitOnlineQuadBezier(vector<T> &ctrlAry, vector<T> &stderrAry, OnlineCache &online, int sampling)
    {
        int n=sampling+1;
        float sum_t1u1y=0, sum_t1u3=0, sum_t3u1=0, sum_t2u2=0;
        float sum_u4=0, sum_t4=0;

        int i;
        // t: x,  u: 1-x
        for (i=0; i<=sampling; i++) {
            float t = (float)i/sampling;
            float u = 1-t;
            float u2 = u*u;
            float t2 = t*t;
            sum_t1u3 	+= t*u*u2;
            sum_t2u2	+= t2*u2;
            sum_t3u1	+= t*t2*u;
            sum_u4      += u2*u2;
            sum_t4      += t2*t2;

        }

        Timer timer;
        timer.start();

#pragma omp parallel for
        for (i=0; i<online.ysum.size(); i++)
        {
            T &p0 = online.y0[i];
            T &p2 = online.y1[i];
            // control
            T ctrl = (online.ytsum[i] - online.yt2sum[i] - p0*sum_t1u3 - p2*sum_t3u1 ) * (.5/sum_t2u2) ;
            ctrlAry[i] = ctrl;

            // sum est y
            T sum_est_y2; // assume initialized
#if 0
            for (int j=0; j<=sampling; j++)
            {
                float t = (float)j/sampling;
                float u = 1-t;
                T est_y = p0*(u*u) + ctrl*( u*t*2. ) + p2*(t*t);
                sum_est_y2 += mult(est_y, est_y);
            }
#else
            sum_est_y2 = mult(p0,p0)*sum_u4 + mult(ctrl, ctrl)*(4.*sum_t2u2) + mult(p2, p2)*sum_t4 +
                    mult(p0, ctrl)*(4.*sum_t1u3) + mult(p0, p2)*(2.*sum_t2u2) + mult(ctrl, p2)*(4.*sum_t3u1);
#endif

            T t1 = mult(p0, online.ysum[i]);
            T t2 =  mult((ctrl-p0), online.ytsum[i])*2.f;
            T t3 = mult((p0-ctrl*2.f+p2), online.yt2sum[i]);

            T sum_yiyest = mult(p0, online.ysum[i]) + mult((ctrl-p0), online.ytsum[i])*2.f + mult((p0-ctrl*2.f+p2), online.yt2sum[i]);
            stderrAry[i] = sqrt(( online.y2sum[i] - sum_yiyest*2.f + sum_est_y2 ) * (1.f/(n-1)));
        }

        timer.end();
        totalFitTime+=timer.getElapsedUS();
        fitTimes++;
    }


    void run_online(int sampling) {
        fitTimes = 0;  addTimes = 0;

        int i, j, z;
        OnlineCache online;

        int size = W*H*D;
        println("sampling=%d", sampling);
        vector<T> flowfield(size);
        initOnlineQuadBezier(online, size);
        // output
        this->quadBezierAry = vector<T>(size);
        this->rmsErrAry = vector<T>(size);


        string out_path = GET_ARG_STRING("out_path").c_str();
        int count = 0;
        int sample_base = 0;
    #if 1
        // the simulation loop
        initOnlineQuadBezier(online, size);
        for (i = 0; i < Ts; i ++)
        {

            // start file reading
            {
                println("Opening %s", fileAry[i].c_str());
                FILE *fp = fopen(fileAry[i].c_str(), "rb");
                if (!fp) {
                    perror("load flow field");
                    exit(1);
                }
                if (getFileExtension( fileAry[i] )=="vec") {
                    // skip first 3 integers
                    fseek(fp, 12, SEEK_SET);
                }

                printf("Reading data...");
                fread(&flowfield[0], sizeof(T), W * H * D, fp);
                printf("Done\n");

                fclose(fp);
            }


            addOnlineQuadBezier(online, flowfield, (i-sample_base)/(float)sampling);

            if (i%sampling==0) {
                if (i>0) {
                    fitOnlineQuadBezier(this->quadBezierAry, this->rmsErrAry , online, sampling);
                    vector<vector<T> > fittedAry; // dum
                    if (!out_path.empty())
                        saveFittedFlowfields(fittedAry, sample_base, sampling, out_path.c_str());
                    count ++;
                }

                initOnlineQuadBezier(online, size);
                sample_base = i;

                // start saving
                if (!out_path.empty()) {
                    std::string cmd ;
                    if (DIMS==3)
                        cmd = strprintf("cp %s %s/sampling%d_%02d.vec", fileAry[i].c_str(), out_path.c_str(), sampling, i);
                    else
                        cmd = strprintf("cp %s %s/sampling%d_%02d.raw", fileAry[i].c_str(), out_path.c_str(), sampling, i);
                    println("%s", cmd.c_str());
                    system(cmd.c_str());
                }

            }
            reportTime();
        }
    #else
        for (i = 0; i < Ts; i += sampling)
        {
            if (i+sampling >= Ts) {
                break;
            }
            count ++;
        }

    #endif

        FILE *fp = fopen(strprintf("%s/all_bezier_rms.list", out_path.c_str()).c_str(), "w");
        fprintf(fp, "%d %d %d %d\n", W, H, D, count);
        fprintf(fp, "%lg\n", 1.f/(double)sampling);
        for (i=0; i< count ; i++)
        {
            fprintf(fp, "sampling%d_%02d.vec sampling%d_%02d_bctrl.raw sampling%d_%02d_rmserr.raw\n",
                sampling, sampling*i, sampling, sampling*i, sampling, sampling*i);

            int i1 = i+1;
            if (i==count-1)
                fprintf(fp, "sampling%d_%02d.vec sampling%d_%02d_bctrl.raw sampling%d_%02d_rmserr.raw\n",
                    sampling, sampling*i1, sampling, sampling*i, sampling, sampling*i);

        }
        fclose(fp);
    }

    void test_online() {
        //float seq[] = {1, 2, 3, 100, 5, 6, 7, 8, 9, 10};
        float seq[] = {1, 4, 9, 100, 25, 36, 49, 64, 81, 100};
        int sampling = 9, n = 10;
        int i;
        vector<vector<T> > fieldAry(n, vector<T>(1));
        OnlineCache online;
        vector<T> ctrlAry(1);
        stdErrAry =  vector<T> (1);
        quadBezierAry = vector<T> (1);
        rmsErrAry =  vector<T> (1);
        meanErrAry = vector<T> (1);
        W=H=D=1;

        initOnlineQuadBezier(online, 1);
        for (i=0; i<=sampling; i++)
        {
            fieldAry[i][0][0] = seq[i];
            addOnlineQuadBezier(online, fieldAry[i], (float)i/sampling);
        }

        fitOnlineQuadBezier(quadBezierAry, rmsErrAry, online, sampling);
        printf("online:  Ctrl=%f, stderr=%f\n",  quadBezierAry[0][0], rmsErrAry[0][0]);

        fitQuadBezier(fieldAry, 0, sampling);
        printf("offline: Ctrl=%f, stderr=%f\n", quadBezierAry[0][0], rmsErrAry[0][0]);

    }

    void reportTime() {
        printf("sampling=%d, threads=%d, Fit times: %d\n", GET_ARG_INT("sampling"), GET_ARG_INT("threads"), fitTimes);
        printf("Average online add time, fit time (ms): %.5lf, %.5lf\n", totalAddTime.getValue()*1e-3/addTimes, totalFitTime.getValue()*1e-3/fitTimes );
    }
};

int main(int argc, const char **argv) {
    CmdArgReader::init(argc, argv, "-sampling=4 -list=all.list -2d=0 -errhist=1 -scalar=0 -threads=4");

    omp_set_num_threads(GET_ARG_INT("threads"));

#if 0 // debug
    {
        ErrorModeling<SINGLETON, 1>em;
        em.test_online();
        return 0;
    }
#endif

	saveErrHist = GET_ARG_INT("errhist");

	load_list(GET_ARG_STRING("list"));

    int sampling = GET_ARG_INT("sampling");

	//OSUFlow osuflow;
	//osuflow.LoadData(GET_ARG_STRING("list").c_str(), false); // time-varying
    if (GET_ARG_INT("scalar")) {
        SINGLETON x; //test
        x[0]=1;
        ErrorModeling<SINGLETON, 1> em;
        em.run_online(sampling);
        em.reportTime();
    } else {
        ErrorModeling<VECTOR3, 3> em;
        em.run_online(sampling);
        em.reportTime();
    }
}
