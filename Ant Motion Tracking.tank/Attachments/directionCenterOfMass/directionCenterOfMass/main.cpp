// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPointValueCollection2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTVectorCollection2D Computation(const DTPointValueCollection2D &centerMass);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTPointValueCollection2D centerMass;
    Read(inputFile,"centerMass",centerMass);
    
    // The computation.
    DTVectorCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(centerMass);
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

DTVectorCollection2D Computation(const DTPointValueCollection2D &centerMass)
{
    
    // To extract a single entry, you can use
    // DTPointCollection Pts = ...
    // for (int i = 0;i<Pts.NumberOfPoints();i++) {
    //     DTPoint2D p = Pts(i);
    //     ....
    // }
    //
    // More efficient
    // DTDoubleArray data = Pts.Data();
    // DTPoint2D p;
    // for (int i = 0;i<data.n();i++) {
    //     p.x = data(0,i);
    //     p.y = data(1,i);
    //     ....
    // }
    
    DTPointCollection2D points = centerMass;
    for (int i = 0; i <points.NumberOfPoints(); i++) {
        DTPoint2DValue point = points(i);
        double xOrig = point(0,i);
        double yOrig;
        double xNext;
        double yNext;
        
    }
    
    
    
    
    DTVectorCollection2D toReturn;
    
    // centerMass.Points() - DTPointCollection2D - where values are defined.
    // centerMass.Points().Data() - DTDoubleArray - A 2xN array.
    // centerMass.Values() - DTDoubleArray - A list with N entries.
    
    return toReturn;
}
