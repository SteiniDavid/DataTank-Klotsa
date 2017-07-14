// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"
#include "DTPath2DValues.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTDoubleArray Computation(double dt,const DTPath2DValues &path);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    double dt = inputFile.ReadNumber("dt");
    DTPath2DValues path;
    Read(inputFile,"path",path);
    
    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(dt,path);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("Array","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

DTDoubleArray Computation(double dt,const DTPath2DValues &path)
{
 
    DTPath2D thePath = path.Path();
    DTDoubleArray values = path.Values();
    int numPoints = thePath.NumberOfPoints();
    
    DTMutableDoubleArray toReturn(4,numPoints);
    int posInReturn = 0;
    
    DTDoubleArray loops = thePath.Data();
    int loc = 0;
    int len = loops.n();
    int ptN,loopLength,loopStarts,loopEnds;
    double x,y;
    int loopN = 1;
    
    while (loc<len) {
        loopLength = int(loops(1,loc));
        loopStarts = loc+1;
        loopEnds = loc+1+loopLength;
        // loopClosed = (loops(0,loopStarts)==loops(0,loopEnds-1) && loops(1,loopStarts)==loops(1,loopEnds-1));
        for (ptN=loopStarts;ptN<loopEnds;ptN++) {
            x = loops(0,ptN);
            y = loops(1,ptN);
            toReturn(0,posInReturn) = x;
            toReturn(1,posInReturn) = y;
            toReturn(2,posInReturn) = loopN;
            toReturn(3,posInReturn) = values(ptN); //how do i show in finder most easily. 
            posInReturn++;
        }
        loc = loopEnds; // Prepare for the next loop.
        loopN++;
    }

    
    
    // path.Path() - DTPath2D.
    // path.Values() - DTDoubleArray - same packing as the path.
    
    return toReturn; //well might as well include time since it cant hurt right? yes
}
