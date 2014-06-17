
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

void computeMaxError(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
    println ("Compute Max Error");
    vector<float> maxErrField(W * H * D);
    int x, y, z, i;
    for (z = 0; z < D; z++)
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++) {
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
    fwrite(&maxErrField[0], sizeof(float), W*H*D,  fp);
    fclose(fp);

    string nrrd_fname;
    write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), W, H, D, "float");
    println("Done");
}


void computeMaxAngleMagError(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
    println ("Compute Max angle/mag Error");
    vector<float> maxMagErrField(W * H * D);
    vector<float> maxAngErrField(W * H*  D);
    int x, y, z, i;
    for (z = 0; z < D; z++)
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++) {
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
    fwrite(&maxMagErrField[0], sizeof(float), W*H*D, fp);
    fclose(fp);

    write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), W, H, D, "float");


    fname = strprintf("maxAngErrField%d-%d", sampling, start).c_str();
    fp = fopen( ("err/"+fname+".raw").c_str(), "wb");
    assert(fp);
    fwrite(&maxAngErrField[0], sizeof(float), W*H*D, fp);
    fclose(fp);

    string nrrd_fname;
    write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), W, H, D, "float");
    println("Done");

}

void genConvex4Error(vector<vector<VECTOR3> > &flowfieldAry, int start, int sampling)
{
    println ("Compute convex4 Error");
    vector<float> errProximateRatioField(W*H*D);
    int x, y, z, i, dim;
    //d=2; //!!!!!!!!!!
    for (z = 0; z < D; z++)
    {
        if (y%10==0)
            println("z=%d", z);
        for (y = 0; y < H; y++)
            for (x = 0; x < W; x++)
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
    fwrite(&errProximateRatioField[0], sizeof(float), W*H*D, fp);
    fclose(fp);

    write_nrrd_3d( ("err/"+fname+".nhdr").c_str(), (fname+".raw").c_str(), W, H, D, "float");

}

