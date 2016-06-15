#include <vector>
#include <cstdio>
#include <cassert>

using namespace std;

template <typename T>
class OnlineQuadBezierFit {

public:
  T ysum, ytsum, yt2sum, y2sum, y0, y1;

  OnlineQuadBezierFit()
    : ysum(0), ytsum(0), yt2sum(0), y2sum(0), y0(0), y1(0)
  {}

  void reset() {
      ysum = ytsum = yt2sum = y2sum = y0 = y1 = 0;
  }

  // t: in range [0..1].  0: first value, 1: last value
  // y: input data value at t
  void addData(const T& y, double t)
  {
    T yt = y*t;
    ysum += y;
    ytsum += yt;
    yt2sum += yt*t;
    y2sum += y*y; // element-wise multiplication for VECTOR3

    if (t==0) {
      y0 = y;
    }
    if (t==1.) {
      y1 = y;
    }

  }

};

// input: online data info , sequence length (n)
// output : ctrlAry, stderrAry
template <typename T>
void fitOnlineQuadBezier(vector<T> &ctrlAry, vector<T> &stderrAry, vector<OnlineQuadBezierFit<T> > &onlineAry, int n)
{
    double sum_t1u1y=0, sum_t1u3=0, sum_t3u1=0, sum_t2u2=0;
    double sum_u4=0, sum_t4=0;

    // allocate output arrays
    int pts = onlineAry.size();
    ctrlAry.clear(); ctrlAry.resize(pts);
    stderrAry.clear(); stderrAry.resize(pts);

    int i;
    // t: x,  u: 1-x
    for (i=0; i<n; i++) {
        double t = (double)i/(n-1);
        double u = 1-t;
        double u2 = u*u;
        double t2 = t*t;
        sum_t1u3 	+= t*u*u2;
        sum_t2u2	+= t2*u2;
        sum_t3u1	+= t*t2*u;
        sum_u4      += u2*u2;
        sum_t4      += t2*t2;

    }

#pragma omp parallel for
    for (i=0; i<pts; i++)
    {
        T &p0 = onlineAry[i].y0;
        T &p2 = onlineAry[i].y1;
        // control
        T ctrl = (onlineAry[i].ytsum - onlineAry[i].yt2sum - p0*sum_t1u3 - p2*sum_t3u1 ) * (.5/sum_t2u2) ;
        ctrlAry[i] = ctrl;

        // sum est y
        T sum_est_y2; // assume initialized
#if 0
        for (int j=0; j<=sampling; j++)
        {
            double t = (double)j/sampling;
            double u = 1-t;
            T est_y = p0*(u*u) + ctrl*( u*t*2. ) + p2*(t*t);
            sum_est_y2 += mult(est_y, est_y);
        }
#else
        sum_est_y2 = (p0*p0)*sum_u4 + (ctrl*ctrl)*(4.*sum_t2u2) + (p2*p2)*sum_t4 +
                p0*ctrl*(4.*sum_t1u3) + (p0*p2)*(2.*sum_t2u2) + (ctrl*p2)*(4.*sum_t3u1);
#endif

        T t1 = p0 * onlineAry[i].ysum;
        T t2 =  (ctrl-p0) * onlineAry[i].ytsum*2.;
        T t3 = (p0-ctrl*2.f+p2) * onlineAry[i].yt2sum;

        T sum_yiyest = (p0*onlineAry[i].ysum) + ((ctrl-p0)*onlineAry[i].ytsum)*2.f + ((p0-ctrl*2.f+p2) * onlineAry[i].yt2sum);
        stderrAry[i] = sqrt(max((T)0, (onlineAry[i].y2sum - sum_yiyest*2.f + sum_est_y2 ) * (1.f/(n-1)))); // max(): neg val happens when result is very close to 0
    }

}
