// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPointCollection2D.h"
#include "DTPointValueCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTPointValueCollection2D Computation(double radius,const DTPointCollection2D &points);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double radius = inputFile.ReadNumber("radius");
    DTPointCollection2D points;
    Read(inputFile,"points",points);
    
    // The computation.
    DTPointValueCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(radius,points);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("PointValueCollection2D","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

DTPointValueCollection2D Computation(double radius,const DTPointCollection2D &points)
{
    DTPointValueCollection2D toReturn;
    
    
    
    return toReturn;
}
