// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPointCollection2D.h"
#include "DTSeriesPointCollection2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTVectorCollection2D Computation(const DTSeriesPointCollection2D &allPoints);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection2D allPoints;
    Read(inputFile,"allPoints",allPoints);
    
    // The computation.
    DTVectorCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(allPoints);
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

struct SingleList
{
    SingleList() {howManyPoints = 0; points = DTMutableDoubleArray(2,10);} // Initialize
    
    DTMutableDoubleArray points; // 2xN array of xy values
    int howManyPoints;
    
    DTPoint2D lastPoint;
    int timeIndex;
    void AddPoint(const DTPoint2D &addThis,int timeI);
};

void SingleList::AddPoint(const DTPoint2D &addThis,int timeI)
{
    if (howManyPoints==points.n()) {
        // at the end of the array
        points = IncreaseSize(points,points.Length());
    }
    points(0,howManyPoints) = addThis.x;
    points(1,howManyPoints) = addThis.y;
    howManyPoints++;
    lastPoint = addThis;
    timeIndex = timeI;
}


DTVectorCollection2D Computation(const DTSeriesPointCollection2D &allPoints)
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
    
    DTMutableList<SingleList> listOfSegments;
    
    DTDoubleArray timeValues = allPoints.TimeValues();
    int timeN, howManyTimes = timeValues.Length();
    DTPointCollection2D current;
    int frameNumber = 0;
    for (timeN=0;timeN<howManyTimes;timeN++) {
        current = allPoints.Get(timeValues(timeN));
        
        // For each point in current, find the closest point in the listOfSegments
        for (int i = 0; i < current.NumberOfPoints(); i++) {
            
            // DTPoint2D point = next(i);
            // double xOrig = point.x;
            double yOrig;
            double xNext;
            double yNext;
            
        }
    
        frameNumber++;
    }
    
    
    
    
    
    DTVectorCollection2D toReturn;
    
    // centerMass.Points() - DTPointCollection2D - where values are defined.
    // centerMass.Points().Data() - DTDoubleArray - A 2xN array.
    // centerMass.Values() - DTDoubleArray - A list with N entries.
    
    return toReturn;
}
