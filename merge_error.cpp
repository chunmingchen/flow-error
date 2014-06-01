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
#include "statistics.h"

#include "PathlineLoader.h"

using namespace std;
using namespace JCLib;

int main(int argc, char **argv)
{
    char *file_fw = argv[1];
    char *file_bw = argv[2];
    PathlineLoader pathline_fw(file_fw);
    PathlineLoader pathline_bw(file_bw);

    pathline_fw.forwardOrder();
    pathline_bw.forwardOrder();

    auto itr_trace_bw = pathline_bw.traceList.begin();

    assert(pathline_fw.traceList.size() == pathline_bw.traceList.size());

    // out
    PathlineLoader out_pathlines;
    SeedTraceList4D &outTraceList = out_pathlines.traceList;

    for ( auto ptrace_fw : pathline_fw.traceList )
    {
        // out
        jcListTimeSeedTrace *outTrace = new jcListTimeSeedTrace;
        outTrace->varList = new list<VECTOR3 *> ;

        auto ptrace_bw = *itr_trace_bw;

        int n = ptrace_fw->size();
        assert( n == ptrace_fw->varList->size() );
        if ( n != ptrace_bw->size() ) {
            // fw and bw finishes differently
            outTraceList.push_back(outTrace);
            ++itr_trace_bw;
            printf("Trace length not match\n");
            continue;
        }
        assert( n == ptrace_bw->varList->size() );

        auto itr_point_fw = ptrace_fw->begin();
        auto itr_point_bw = ptrace_bw->begin();
        auto itr_var_fw = ptrace_fw->varList->begin();
        auto itr_var_bw = ptrace_bw->varList->begin();



        while( itr_point_fw != ptrace_fw->end() && itr_point_bw != ptrace_bw->end() )
        {
            VECTOR4 &p1 = **itr_point_fw;
            VECTOR3 &v1 = **itr_var_fw;
            VECTOR4 &p2 = **itr_point_bw;
            VECTOR3 &v2 = **itr_var_bw;

            if ( p1[3]+1e-3 < p2[3] ) {
                ++ itr_point_fw; ++itr_var_fw;
                printf("Trace time not match %f %f\n", p1[3], p2[3]);
                continue;
            }
            if ( p2[3]+1e-3 < p1[3] ) {
                ++ itr_point_bw; ++itr_var_bw;
                printf("(Trace time not match: %f %f)", p1[3], p2[3]);
                continue;
            }
            if (fabs(p1[3]-p2[3] ) > 1e-5) {
                println("(Allow time diff %f)", fabs(p1[3]-p2[3]));
            }



//            if( fabs(p1[3] - p2[3]) > 1e-5 ) {
//                println("time not match! %f, %f", p1[3], p2[3]);
//                exit(1);
//            }

#ifdef _DEBUG
            printf(" p1:"); print_array(&p1[0], 4);
            printf(" v1:"); print_array(&v1[0], 3);
            printf(" p2:"); print_array(&p2[0], 4);
            printf(" v2:"); print_array(&v2[0], 3);
#endif
            VECTOR4 outp(0,0,0,p1[3]);
            VECTOR3 outv;
            for (int d=0; d<3; d++)
            {
                if (v1[d]+v2[d]==0) {
                    if (p2[d] != p1[d]) {
                        JCLib::print_array(&p1[0], 4);
                        JCLib::print_array(&p2[0], 4);
                        printf("p1!=p2!! d=%d\n", d);
                        exit(1);
                    }

                    outv[d] = v1[d];
                    outp[d] = p1[d];
                } else {
                    outv[d] = v1[d]*v2[d]/ (v1[d]+v2[d]);
                    outp[d] = (p1[d]*v2[d] + p2[d]*v1[d]) / (v1[d]+v2[d]);  //outp[d] = outv[d] * (p1[d]/(double)v1[d] + p2[d]/(double)v2[d]);
                }
            }
#ifdef _DEBUG
            printf(" outp':"); print_array(&outp[0], 3);
            printf(" outv':"); print_array(&outv[0], 3);
#endif

            //VECTOR3 vsum = v1+v2;
            //double p1p2_vsum = v1[0]*v2[0]/vsum[0] + v1[1]*v2[1]/vsum[1] + v1[2]*v2[2]/vsum[2]  ;
            //double c = exp(-0.5 * p1p2_vsum)  / sqrt( pow(2*M_PI,3) * dot(vsum, vsum)) ;
            //outp = outp * c;
            //outv = outv * c;
            //printf(" outp:"); print_array(&outp[0], 3);
            //printf(" outv:"); print_array(&outv[0], 3);
            //printf("\n");

            outTrace->push_back(new VECTOR4(outp));
            outTrace->varList->push_back(new VECTOR3(outv));

            ++itr_point_fw;
            ++itr_var_fw;
            ++itr_point_bw;
            ++itr_var_bw;
        }

        outTraceList.push_back(outTrace);
        ++itr_trace_bw;
    }

    if (argc>3)
        out_pathlines.saveDataNormal(argv[3]);
    else
        out_pathlines.saveDataNormal("merged.out");
    out_pathlines.saveDataVTK("merged.vtp");
}
