// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPath2D.h"
#include "DTPointCollection2D.h"
#include "DTSeriesPointCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTPath2D Computation(const DTSeriesPointCollection2D &allPoints,double maxDist);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection2D allPoints;
    Read(inputFile,"allPoints",allPoints);
    double maxDist = inputFile.ReadNumber("maxDist");
    
    // The computation.
    DTPath2D computed;
    clock_t t_before = clock();
    computed = Computation(allPoints,maxDist);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("Path2D","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

class SingleList
{
public:
    SingleList() {howManyPoints = 0; points = DTMutableDoubleArray(2,10);} // Initialize
    
    DTMutableDoubleArray points; // 2xN array of xy values
    int howManyPoints;
    
    DTPoint2D lastPoint;
    int timeIndex;
    void AddPoint(const DTPoint2D &addThis,int timeI)
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
    
    DTPath2D ConvertToPath();
};

DTPath2D SingleList::ConvertToPath()
{
    // The first column needs to be the length
    DTMutableDoubleArray loop(2,howManyPoints+1);
    // [ 0 x1 .... xN]
    // [ N y1 .... yN]
    loop(0,0) = 0;
    loop(1,0) = howManyPoints;
    
    for (int i=0;i<howManyPoints;i++) {
        loop(0,i+1) = points(0,i);
        loop(1,i+1) = points(1,i);
    }
    
    return DTPath2D(loop);
}

DTPath2D Computation(const DTSeriesPointCollection2D &allPoints,double maxDist)
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
    
    DTMutableList<SingleList> listOfSegments(10);
    int howManySegments = 0;
    
    DTDoubleArray timeValues = allPoints.TimeValues();
    int timeN, howManyTimes = timeValues.Length();
    DTPointCollection2D current;
    int frameNumber = 0;
    DTMutableDoubleArray dist(100,100);
    
    for (timeN=0;timeN<howManyTimes;timeN++) {
        current = allPoints.Get(timeValues(timeN));
        
        // For each point in current, find the closest point in the listOfSegments
        if (current.IsEmpty()) continue;
        // For each point in the segments, find a close match in the current
        for (int j = 0; j < howManySegments; j++) {
            for (int k = 0; k < current.NumberOfPoints(); k++) {
                DTPoint2D currentPoint = current(k);
                dist(j,k) = (currentPoint.x-listOfSegments(j).points(0,j))+(currentPoint.y-listOfSegments(j).points(1,j));
            }
        }
        // Compute a distance matrix
        // dist(i,j) = distance from point i in the current point collection to point j in the segment list
        
        // Find I,J such that dist(I,J) is the minimum distance
        // Connect points I and J, i.e. point I in the current should be added to point J in the segment list
        // and then set dist(I,:) and dist(:,J) to the largest value (maxDist) so that it won't be picked again.
        // Repeat the last 3 lines until the minimum distance is >= maxDist
        
        // Any leftover points need to be added to the segments
        
        for (int i = 0; i < current.NumberOfPoints(); i++) {
            
            // DTPoint2D currentPoint = current.DTPoint2D();
            //DTPoint2D lastPoint = listOfSegments.
            // double xOrig = point.x;
            double yOrig;
            double xNext;
            double yNext;
            
        }
        
        frameNumber++;
    }
    
    // At this point go over the segments, tie them together into a path.
    DTMutableList<DTPath2D> segmentPaths(howManySegments);
    int i;
    for (i=0;i<howManySegments;i++) {
        segmentPaths(i) = listOfSegments(i).ConvertToPath();
    }
    return Combine(segmentPaths);
}
