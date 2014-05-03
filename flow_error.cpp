#include <algorithm>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include "OSUFlow.h"
#include "system/cmd_arg_reader.h"
#include "system/path.h"
#include "macros.h"
#include "file/nrrd.h"
#include "statistics.h"
using namespace std;
using namespace JCLib;

vector<string> fileAry;
int w, h, d, t;
#define POS_ID(x,y,z) ((x)+w*((y)+h*(z)));
bool saveErrHist=false;

//#define DEBUG_CONVEX

// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
float det( VECTOR2 &O,  VECTOR2 &A,  VECTOR2 &B)
{
	return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
}

struct VECTOR2_comp{
	bool operator() ( const VECTOR2 &p1_,  const VECTOR2 &p2_) {
		VECTOR2 p1 = p1_;
		VECTOR2 p2 = p2_;
		return p1[0] < p2[0] || (p1[0]==p2[0] && p1[1]<p2[1]);
	}
};

// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
vector<VECTOR2> convex_hull(vector<VECTOR2> P, int *lower_hull_size)
{
	int n = P.size(), k = 0;
	vector<VECTOR2> H(2*n);

	// Sort points lexicographically
	sort(P.begin(), P.end(), VECTOR2_comp());

	// Build lower hull
	for (int i = 0; i < n; i++) {
		while (k >= 2 && det(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}
	*lower_hull_size = k;

	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
		while (k >= t && det(H[k-2], H[k-1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k);
	return H;
}


VECTOR2 intersect(VECTOR2 &p1, VECTOR2 &p2, VECTOR2 &p3, VECTOR2 &p4) {
// Store the values for fast access and easy
// equations-to-code conversion
float x1 = p1[0], x2 = p2[0], x3 = p3[0], x4 = p4[0];
float y1 = p1[1], y2 = p2[1], y3 = p3[1], y4 = p4[1];

float d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
// If d is zero, there is no intersection
if (d == 0) return p2;

// Get the x and y
float pre = (x1*y2 - y1*x2), post = (x3*y4 - y3*x4);
float x = ( pre * (x3 - x4) - (x1 - x2) * post ) / d;
float y = ( pre * (y3 - y4) - (y1 - y2) * post ) / d;

// Return the point of intersection
return VECTOR2(x,y);
}

////////////////////////////////////////////////////////

// global: fileAry, w,h,d,t
void load_list(string list_filename) {
	int i, j;
	string path = getPath(list_filename);

	FILE *fp = fopen(list_filename.c_str(), "rt");
	char str[1024];
	fgets(str, 1024, fp);
	i = sscanf(str, "%d %d %d %d", &w, &h, &d, &t);
	assert(i==4);
	fgets(str, 1024, fp); // dummy line


	fileAry.clear();
	for (i = 0; i < t; i++) {
		fgets(str, 1024, fp);
		char *p = strtok(str, " \r\n");
		//printf("%s\n", p);

		if (isFilenameOnly(p))
			fileAry.push_back(path + p);
		else
			fileAry.push_back(p);
	}
}

void computeMaxError(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	println ("Compute Max Error");
	vector<float> maxErrField(w * h * d);
	int x, y, z, i;
	for (z = 0; z < d; z++)
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++) {
				int id = POS_ID(x,y,z);
				VECTOR3 dist = flowfieldAry[sampling][id]
						- flowfieldAry[0][id];
				vector<float> errAry(sampling);
				for (i = 0; i < sampling; i++) {
					VECTOR3 off = dist * (float(i)/sampling);
					VECTOR3 &x0 = flowfieldAry[0][id];
					VECTOR3 interp(x0[0] + off[0], x0[1] + off[1],
							x0[2] + off[2]);
					errAry[i] = (interp - flowfieldAry[i][id]).GetMag();
					if (isinf(errAry[i]))
						errAry[i] = -1;
				}
				maxErrField[id] = *max_element(errAry.begin(), errAry.end());
			}

	string fname;
	fname = strprintf("maxErrField%d-%d", sampling, start).c_str();
	FILE *fp = fopen( ("err/"+fname+".raw").c_str(), "wb");
	assert(fp);
	fwrite(&maxErrField[0], sizeof(float), w * h * d, fp);
	fclose(fp);

	string nrrd_fname;
	write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), w, h, d, "float");
	println("Done");
}


void computeMaxAngleMagError(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	println ("Compute Max angle/mag Error");
	vector<float> maxMagErrField(w * h * d);
	vector<float> maxAngErrField(w * h * d);
	int x, y, z, i;
	for (z = 0; z < d; z++)
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++) {
				int id = POS_ID(x,y,z);
				VECTOR3 dist = flowfieldAry[sampling][id]
						- flowfieldAry[0][id];
				vector<float> magErrAry(sampling);
				vector<float> angErrAry(sampling);
				for (i = 0; i < sampling; i++) {
					VECTOR3 off = dist * (float(i)/sampling);
					VECTOR3 &x0 = flowfieldAry[0][id];
					VECTOR3 interp(x0[0] + off[0], x0[1] + off[1],
							x0[2] + off[2]);
					magErrAry[i] = fabs(interp.GetMag() - flowfieldAry[i][id].GetMag());
					if (isinf(magErrAry[i]))
						magErrAry[i] = -1;

					float costh = dot(interp, flowfieldAry[i][id]) / interp.GetMag() / flowfieldAry[i][id].GetMag();
					angErrAry[i] = acos(costh)*180.0/PI;
					if (isinf(angErrAry[i]) || isnan(angErrAry[i]))
						angErrAry[i] = -1;

				}
				maxMagErrField[id] = *max_element(magErrAry.begin(), magErrAry.end());
				maxAngErrField[id] = *max_element(angErrAry.begin(), angErrAry.end());
				//printf("%f\n", maxAngErrField[id]);
			}

	string fname;
	fname = strprintf("maxMagErrField%d-%d", sampling, start).c_str();
	FILE *fp = fopen( ("err/"+fname+".raw").c_str(), "wb");
	assert(fp);
	fwrite(&maxMagErrField[0], sizeof(float), w * h * d, fp);
	fclose(fp);

	write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), w, h, d, "float");


	fname = strprintf("maxAngErrField%d-%d", sampling, start).c_str();
	fp = fopen( ("err/"+fname+".raw").c_str(), "wb");
	assert(fp);
	fwrite(&maxAngErrField[0], sizeof(float), w * h * d, fp);
	fclose(fp);

	string nrrd_fname;
	write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), w, h, d, "float");
	println("Done");

}

void genConvex4Error(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	println ("Compute convex4 Error");
	vector<float> errProximateRatioField(w*h*d);
	int x, y, z, i, dim;
	//d=2; //!!!!!!!!!!
	for (z = 0; z < d; z++)
	{
		if (y%10==0)
			println("z=%d", z);
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++)
			{
				int id = POS_ID(x,y,z);
				// interp
				VECTOR3 dist = flowfieldAry[sampling][id]
						- flowfieldAry[0][id];
				vector<VECTOR3> interpAry(sampling+1);
				for (i = 0; i <= sampling; i++) {
					VECTOR3 off = dist * (float(i)/sampling);
					VECTOR3 &x0 = flowfieldAry[0][id];
					VECTOR3 interp(x0[0] + off[0], x0[1] + off[1],
							x0[2] + off[2]);
					interpAry[i] = interp;
				}

				vector<float> errProximateRatioAry;
				for (dim=0; dim<3; dim++)
				{
					vector<VECTOR2> convex;
					convex.push_back(VECTOR2(0, flowfieldAry[0][id][dim]));
					convex.push_back(VECTOR2(1, flowfieldAry[1][id][dim]));
					int lower_hull_size, upper_hull_size;
					for (i=2; i<=sampling; i++)
					{
						convex.push_back(VECTOR2(i, flowfieldAry[i][id][dim]));

						convex = convex_hull(convex, &lower_hull_size);
						upper_hull_size = convex.size()-lower_hull_size+1;
#ifdef DEBUG_CONVEX
						printf("lower=%d upper=%d\n", lower_hull_size, upper_hull_size);
#endif

						if (upper_hull_size >= 5)
						{
							//println("upper cut");
							int s = convex.size();
							//   merge last 4 points -> 3 points
							//println("%f %f ", convex[s-4][0], convex[s-4][1]);
							//println("%f %f ", convex[s-3][0], convex[s-3][1]);
							//println("%f %f ", convex[s-2][0], convex[s-2][1]);
							//println("%f %f ", convex[s-1][0], convex[s-1][1]);
							convex[s-3] = intersect(convex[s-4], convex[s-3], convex[s-2], convex[s-1]);;
							//println("->%f %f ", convex[s-3][0], convex[s-3][1]);
							convex.erase(convex.begin()+(s-2));
							upper_hull_size--;
						}

						if (lower_hull_size >= 5)
						{
							//println("lower cut");
							// merge first 4 points -> 3 points
							convex[1] = intersect(convex[0], convex[1], convex[2], convex[3]);;
							convex.erase(convex.begin()+2);
							lower_hull_size --;
						}

					}
#ifdef DEBUG_CONVEX
					printf("convex:\n");
					for (i=0; i<convex.size(); i++)
						println("%f %f", convex[i][0], convex[i][1]);
					for (i=0; i<=sampling; i++)
						println("%f", flowfieldAry[i][id][dim]);
#endif
					// convex: <lower> <upper>
					vector<float> valueLowAry(sampling+1), valueHighAry(sampling+1);
					for (i=1; i<convex.size(); i++)
					{
						for (int k=convex[i-1][0]; k<=convex[i][0]; k++)
						{
							VECTOR2 &left = convex[i-1];
							VECTOR2 &right = convex[i];
							float y = left[1]+(k-left[0])/(right[0]-left[0])*(right[1]-left[1]);
							valueLowAry[k] =  y ;
						}
						if (convex[i][0]==sampling)
							break;
					}
					for (i--; i<convex.size(); i++)
					{
						for (int k=convex[i-1][0]; k>=convex[i][0]; k--)
						{
							VECTOR2 &left = convex[i-1];
							VECTOR2 &right = convex[i];
							float y = left[1]+(k-left[0])/(right[0]-left[0])*(right[1]-left[1]);
							valueHighAry[k] = y;
						}
					}

#ifdef DEBUG_CONVEX
					printf("reconstruct low:\n");
					for (i=0; i<=sampling; i++)
						println("%f", valueLowAry[i]);
					printf("reconstruct high:\n");
					for (i=0; i<=sampling; i++)
						println("%f", valueHighAry[i]);
#endif
					for (i=1; i<sampling; i++) {
						float real_err = fabs(flowfieldAry[i][id][dim]-interpAry[i][dim]);
						float approximate = fabs(valueHighAry[i]-valueLowAry[i]);
						float ratio = real_err/approximate;
						if (!isinf(ratio))
							errProximateRatioAry.push_back(real_err/approximate);
					}

				}
				errProximateRatioField[id] = getMean(errProximateRatioAry.begin(), errProximateRatioAry.end());
			}
	}


	string fname;
	fname = strprintf("meanProximateRatio%d-%d", sampling, start).c_str();
	FILE *fp = fopen( ("err/"+fname+".raw").c_str(), "wb");
	assert(fp);
	fwrite(&errProximateRatioField[0], sizeof(float), w * h * d, fp);
	fclose(fp);

	write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), w, h, d, "float");

}


vector<vector<VECTOR3 > > quadFuncAry; // xyz, abc, dim
vector<VECTOR3> rmsErrAry;  // mean square error;  xyz, dim
vector<VECTOR3> stdErrAry;  // std error;  xyz, dim
vector<VECTOR3> meanErrAry;  // std error;  xyz, dim
vector<vector<VECTOR3> > errAry;  // std error;  xyz, dim
// fitting with quad with connected end points
void fitQuadFuncEndPoints(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	int n=sampling+1;
	quadFuncAry = vector<vector<VECTOR3> > (w*h*d, vector<VECTOR3>(3));
	stdErrAry =  vector<VECTOR3> (w*h*d);
	rmsErrAry =  vector<VECTOR3> (w*h*d);
	errAry = vector<vector<VECTOR3> >(w*h*d, vector<VECTOR3>(n));
	meanErrAry = vector<VECTOR3> (w*h*d);

	println ("Fitting quadratic function");
	int x, y, z, i;
	for (z = 0; z < d; z++)
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++)
			{
				int id = POS_ID(x,y,z);
				vector<VECTOR3> y_ary(n);
				for (i=0; i<n; i++)
					y_ary[i] = flowfieldAry[i][id];

				for (int d=0; d<3; d++) //dim
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
					double rms = sqrt(JCLib::getMean(err_ary2.begin(), err_ary2.end()));
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

void fitQuadFunc(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	int n=sampling+1;
	quadFuncAry = vector<vector<VECTOR3> > (w*h*d, vector<VECTOR3>(3));
	stdErrAry =  vector<VECTOR3>  (w*h*d);
	rmsErrAry = vector<VECTOR3> (w*h*d);

	println ("Fitting quadratic function");
	int x, y, z, i;
	for (z = 0; z < d; z++)
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++)
			{
				int id = POS_ID(x,y,z);
				vector<VECTOR3> y_ary(n);
				for (i=0; i<n; i++)
					y_ary[i] = flowfieldAry[i][id];

				vector<VECTOR3> &par = quadFuncAry[id];
				for (int d=0; d<3; d++) //dim
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
vector<VECTOR3> quadBezierAry;
void fitQuadBezier(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
	int n=sampling+1;
	quadBezierAry = vector<VECTOR3> (w*h*d);
	stdErrAry =  vector<VECTOR3> (w*h*d);
	rmsErrAry =  vector<VECTOR3> (w*h*d);
	errAry = vector<vector<VECTOR3> >(w*h*d, vector<VECTOR3>(n));
	meanErrAry = vector<VECTOR3> (w*h*d);

	println ("Fitting quadratic function");
	int x, y, z, i;
	for (z = 0; z < d; z++)
		for (y = 0; y < h; y++)
			for (x = 0; x < w; x++)
			{
				int id = POS_ID(x,y,z);
				vector<VECTOR3> y_ary(n);
				for (i=0; i<n; i++)
					y_ary[i] = flowfieldAry[i][id];

				for (int d=0; d<3; d++) //dim
				{
					float y1=y_ary[0][d];
					float yn=y_ary[n-1][d];
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
					quadBezierAry[id][d] = ctrl;

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
					double rms = sqrt(JCLib::getMean(err_ary2.begin(), err_ary2.end()));
					stdErrAry[id][d] = std;
					rmsErrAry[id][d] = rms;
					meanErrAry[id][d] = mean;
					//println("std=%lf mean=%lf rms=%f", std, mean, rms);

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


void saveFittedFlowfields(vector<vector<VECTOR3> >& funcAry, int start, int sampling, const char *folder)
{
	println("Saving fitted flow fields...");
	int i;
	FILE *fp;
#if 0
#if 1
	int skip = sampling; //1
#endif
	vector<VECTOR3> flowfield(w*h*d);
	for (i=0; i<=sampling; i+=skip)
	{
		// gen field
		int id;
		for (id=0; id<w*h*d; id++)
		{
			flowfield[id] = funcAry[id][0]*(i*i) + funcAry[id][1]*i + funcAry[id][2];
		}

		// save file
		fp = fopen(strprintf("%s/sampling%d_%02d.vec", folder, sampling, start+i).c_str(), "wb");
		{
			int dim[3];
			dim[0] = w; dim[1] = h; dim[2] = d;
			fwrite(dim, 4, 3, fp);
		}
		fwrite(&flowfield[0], 4, w*h*d*3, fp);
		fclose(fp);
	}
#endif

	// save quad bezier control point
	fp = fopen(strprintf("%s/sampling%d_%02d_bctrl.raw", folder, sampling, start).c_str(), "wb");
	fwrite(&quadBezierAry[0], 4, w*h*d*3, fp);
	fclose(fp);
#if 1
	// save rms err
	fp = fopen(strprintf("%s/sampling%d_%02d_rmserr.raw", folder, sampling, start).c_str(), "wb");
	fwrite(&rmsErrAry[0], 4, w*h*d*3, fp);
	fclose(fp);
	// save errstd
	fp = fopen(strprintf("%s/sampling%d_%02d_stderr.raw", folder, sampling, start).c_str(), "wb");
	fwrite(&stdErrAry[0], 4, w*h*d*3, fp);
	fclose(fp);
	// save err mean
	fp = fopen(strprintf("%s/sampling%d_%02d_meanerr.raw", folder, sampling, start).c_str(), "wb");
	fwrite(&meanErrAry[0], 4, w*h*d*3, fp);
	fclose(fp);
#endif
#if 0
	// save a
	fp = fopen(strprintf("%s/sampling%d_%02d_a.raw", folder, sampling, start).c_str(), "wb");
	for (int id=0; id<w*h*d; id++)
		fwrite(&funcAry[id][0], 4, 3, fp);
	fclose(fp);
#endif
#if 0
	// errAry
	fp = fopen(strprintf("%s/sampling%d_%02d_err.raw", folder, sampling, start).c_str(), "wb");
	for (int id=0; id<w*h*d; id++)
		fwrite(&errAry[id][0], 4, 3*(sampling+1), fp);
	fclose(fp);
#endif
}

void run(int sampling) {
	int i, j;
	println("allocating size=%d", sampling);
	vector<vector<VECTOR3> > flowfieldAry(sampling+1, vector<VECTOR3>(w * h * d));
	string out_path = GET_ARG_STRING("out_path").c_str();

	for (i = 0; i < t; i += sampling) {
		for (j = 0; j <= sampling; j++) {
			int id = i + j;
			if (id >= t)
				return;

			println("Loading %s", fileAry[id].c_str());
			FILE *fp = fopen(fileAry[id].c_str(), "rb");
			if (!fp) {
				perror("load flow field");
				exit(1);
			}
			int dum[3];
			fread(dum, 4, 3, fp);
			fread(&flowfieldAry[j][0], sizeof(VECTOR3), w * h * d, fp);
			fclose(fp);

			//print_array(&flowfieldAry[j][0][0], 3);
		}

		//computeMaxError(flowfieldAry, i, sampling);
		//computeMaxAngleMagError(flowfieldAry, i, sampling);
		//genConvex4Error(flowfieldAry, i, sampling);

		//fitQuadFunc(flowfieldAry, i, sampling);
		//saveFittedFlowfields(quadFuncAry, i, sampling, "fitted_quad");

		//fitQuadFuncEndPoints(flowfieldAry, i, sampling);
		//saveFittedFlowfields(quadFuncAry, i, sampling, out_path.c_str());

		fitQuadBezier(flowfieldAry, i, sampling);
		if (i==0) {
			std::string cmd = strprintf("cp %s %s/sampling%d_%02d.vec", fileAry[i].c_str(), out_path.c_str(), sampling, i);
			println("%s", cmd.c_str());
			system(cmd.c_str());
		}
		{
			std::string cmd = strprintf("cp %s %s/sampling%d_%02d.vec", fileAry[i+sampling].c_str(), out_path.c_str(), sampling, i+sampling);
			println("%s", cmd.c_str());
			system(cmd.c_str());
		}
		saveFittedFlowfields(quadFuncAry, i, sampling, out_path.c_str());


	}

}

int main(int argc, const char **argv) {
	CmdArgReader::init(argc, argv, "-sampling=4 -list=all.list -2d=0 -errhist=1");

	saveErrHist = GET_ARG_INT("errhist");

	load_list(GET_ARG_STRING("list"));
	//OSUFlow osuflow;
	//osuflow.LoadData(GET_ARG_STRING("list").c_str(), false); // time-varying
	run(GET_ARG_INT("sampling"));
}
