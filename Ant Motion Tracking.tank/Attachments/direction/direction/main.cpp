// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPath2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTVectorCollection2D Computation(const DTPath2D &outline);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTPath2D outline;
    Read(inputFile,"outline",outline);

    // The computation.
    DTVectorCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(outline);
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


/*
 A standard traversal code will look like:
 
 DTDoubleArray loops = path.Data();
 int loc = 0;
 int len = loops.n();
 int ptN,loopLength,loopStarts,loopEnds;
 while (loc<len) {
 loopLength = int(loops(1,loc));
 loopStarts = loc+1;
 loopEnds = loc+1+loopLength;
 // loopClosed = (loops(0,loopStarts)==loops(0,loopEnds-1) && loops(1,loopStarts)==loops(1,loopEnds-1));
 for (ptN=loopStarts;ptN<loopEnds;ptN++) {
 // x = loops(0,ptN);
 // y = loops(1,ptN);
 }
 loc = loopEnds; // Prepare for the next loop.
 }
 */


DTVectorCollection2D Computation(const DTPath2D &outline)
{
    DTVectorCollection2D toReturn;

    int howManyLoops = outline.NumberOfLoops();
    DTMutableDoubleArray centerOfMass(2,howManyLoops);
    for (int i=0;i<howManyLoops;i++) {
        DTPath2D loop = ExtractLoop(outline,i);
        DTDoubleArray points = loop.Data();
        //Finding center of mass
        double xVal = 0, yVal = 0;
        for (int j=1; j < points.n(); j++) {
            xVal += points(0,j);
            yVal += points(1,j);
        }
        double xCenter = xVal/(points.n()-1);
        double yCenter = yVal/(points.n()-1);
        centerOfMass(0,i) = xCenter;
        centerOfMass(1,i) = yCenter;
    }
    
    DTMutableDoubleArray vectorHeads(2,howManyLoops);
    //Next we need to find the longest distance between a center of mass and a particle.
    //We will use that to find which direction the loop is likely pointing.
    //This will be used to fingure out a velocity vector later on.
    for (int i=0;i<howManyLoops;i++) {
        DTPath2D loop = ExtractLoop(outline,i);
        DTDoubleArray points = loop.Data();
        double xVal = 0, yVal = 0;
        double xMax = 0, yMax = 0;
        double xCenter = centerOfMass(0,i);
        double yCenter = centerOfMass(1,i);
        double maxDist = 0;
        for (int j=1; j < points.n(); j++) {
            xVal = points(0,j);
            yVal = points(1,j);
            double dx = xVal-xCenter;
            double dy = yVal-yCenter;
            double dist = dx*dx+dy*dy;
            //I don't know how to extract just the x or j component for a specific index
            if(dist> maxDist){
                xMax = xVal;
                yMax = yVal;
                maxDist = dist;
            }
        }
        
        //Should or when should I make the vector collection
        //Should it be mutable? do I pass in the center of masses in the start
        //
        vectorHeads(0,i) = xMax-xCenter;
        vectorHeads(1,i) = yMax-yCenter;
    }

    DTVectorCollection2D antVectors(DTPointCollection2D(centerOfMass),vectorHeads);
    
    //look through and allign vectors with direction they are going in.
    //If it was going in the wrong direction it should be that the vector is reversed for
    //that loop.
    
    
    
    // DTMutableVectorCollection2D(const DTPointCollection2D &,const DTMutableDoubleArray &);
    
    // outline.Data() - DTDoubleArray - Packed representation (see header file).
    // outline.NumberOfLoops() - how many components there are.
    // outline.NumberOfPoints() - Number of xy values.
    // outline.ExtractLoop(int) - Get a single component.

    return antVectors;
}
