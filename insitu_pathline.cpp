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


int main(int argc, const char **argv)
{
    CmdArgReader::init(argc, argv, "-sampling=4 -list=all.list -2d=0 -errhist=1 -scalar=0 -threads=4");

    //omp_set_num_threads(GET_ARG_INT("threads"));

    OSUFlow osuflow;
    osuflow.LoadData(GET_ARG_STRING("list").c_str(), true);

    // set seeds
    {
        VECTOR3 minb, maxb;
        osuflow.Boundary(minb, maxb);

        int seeds = maxb[0]*maxb[1]*maxb[2];
        vector<VECTOR3> seedAry(seeds);
        int x,y,z;
        for (z=0; z<maxb[2]; z++)
            for (y=0; y<maxb[1]; y++)
                for (x=0; x<maxb[0]; x++)
                {
                    seedAry.push_back(VECTOR3(x+.5f,y+.5f,z+.5f));
                }

        osuflow.SetSeedPoints(&seedAry[0], num_seeds);
    }

    osuflow.GenPathLines(listSeedTraces, FORWARD, seeds);
}
