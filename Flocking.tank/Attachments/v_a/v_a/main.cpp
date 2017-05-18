// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

double Computation(const DTVectorCollection2D &v);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTVectorCollection2D v;
    Read(inputFile,"v",v);

    // The computation.
    double computed;
    clock_t t_before = clock();
    computed = Computation(v);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("Real Number","Seq_Var");

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    return 0;
}

double Computation(const DTVectorCollection2D &v)
{
    double toReturn;

    // v.Points() - DTPointCollection2D - underlying points
    // v.Points().Data() - DTDoubleArray - A 2xN array
    // v.Points().NumberOfPoints() - int - N
    // v.Vectors() - DTDoubleArray - A 2xN array

    DTDoubleArray vectors = v.Vectors();
    int N = v.Points().NumberOfPoints();
    double vmag = sqrt(vectors(0,0)*vectors(0,0)+vectors(1,0)*vectors(1,0));
    
    double v_a_frontpart = (1/(N*vmag));
    double v_a_x = 0;
    double v_a_y = 0;
    for (int i = 0; i < N; i++) {
        v_a_x += vectors(0,i);
        v_a_y += vectors(1,i);
    }
    toReturn = v_a_frontpart*sqrt(v_a_x*v_a_x + v_a_y*v_a_y);
    
    return toReturn;
}
