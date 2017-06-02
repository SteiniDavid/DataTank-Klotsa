// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTRegion2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTVectorCollection2D Computation(const DTRegion2D &box,int N,double v,int seed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTRegion2D box;
    Read(inputFile,"box",box);
    int N = int(inputFile.ReadNumber("N"));
    double v = inputFile.ReadNumber("v");
    int seed = int(inputFile.ReadNumber("seed"));
    
    // The computation.
    DTVectorCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(box,N,v,seed);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("VectorCollection2D","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

#include "DTRandom.h"

DTVectorCollection2D Computation(const DTRegion2D &box,int N,double v,int seed)
{
    DTVectorCollection2D toReturn;
    
    // box.isSet - false means invalid.
    // box.xmin, xmax, ymin, ymax
    
    DTMutableDoubleArray points(2,N);
    DTMutableDoubleArray vectors(2,N);
    
    double L = box.xmax;
    double angle;
    
    DTRandom R(seed);
    for (int i=0; i<N; i++) {
        points(0,i) = R.UniformHalf()*L;
        points(1,i) = R.UniformHalf()*L;
        
        angle = R.UniformHalf()*2*M_PI;
        vectors(0,i) = cos(angle)*v;
        vectors(1,i) = sin(angle)*v;
    }
    
    toReturn = DTVectorCollection2D(DTPointCollection2D(points), vectors);
    
    return toReturn;
}

