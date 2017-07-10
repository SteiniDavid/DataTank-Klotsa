// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPath2DValues.h"
#include "DTPointCollection2D.h"
#include "DTSeriesPointCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTPath2DValues Computation(const DTSeriesPointCollection2D &allPoints,double maxDist);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection2D allPoints;
    Read(inputFile,"allPoints",allPoints);
    double maxDist = inputFile.ReadNumber("maxDist");

    // The computation.
    DTPath2DValues computed;
    clock_t t_before = clock();
    computed = Computation(allPoints,maxDist);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("Path2DValues","Seq_Var");

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
    SingleList() {
        howManyPoints = 0; points = DTMutableDoubleArray(2,10);
    } // Initialize
    
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
    
    DTPath2DValues ConvertToPath(DTDoubleArray timeValues){ //i made a new version.
        // The first column needs to be the length
        DTMutableDoubleArray loop(2,howManyPoints+1);
        // [ 0 x1 .... xN]
        // [ N y1 .... yN]
        
        //[ 0 x1 .... xN 0 x1 ... xM ...]
        //[ N y1 .... yN M y1 ... yM ...]
        loop(0,0) = 0;
        loop(1,0) = howManyPoints;
        
        for (int i=0;i<howManyPoints;i++) {
            loop(0,i+1) = points(0,i);
            loop(1,i+1) = points(1,i);
        }
        
        DTMutableDoubleArray timePlusOne(howManyPoints+1);
        timePlusOne(0) = howManyPoints;
        for (int i=0;i<howManyPoints;i++) {
            timePlusOne(i+1) = timeValues(i);
        }
        
        DTPath2D path = DTPath2D(loop);
        return DTPath2DValues(path, timePlusOne);
    }
    
    //a way of getting a more compact explanation in the debugger
    void pall(void) const {
        cerr << "time = " << timeIndex << endl;
        ExtractColumns(points,DTRange(0,howManyPoints)).pall();
    }
};

DTPath2DValues Computation(const DTSeriesPointCollection2D &allPoints,double maxDist)
{
    DTMutableList<SingleList> listOfSegments(31);
    int howManySegments = 0;
    DTDoubleArray timeValues = allPoints.TimeValues();
    int timeN, howManyTimes = timeValues.Length();
    DTPointCollection2D current;
    DTMutableDoubleArray dist(1000,1000);
    dist = maxDist; //sets everything to max dist initally
    
    //Starting initiation
    int numPoints = allPoints.Get(timeValues(111)).NumberOfPoints();
    for (int i = 0; i < numPoints-1; i++) {
        DTPointCollection2D currentPoints = allPoints.Get(timeValues(111));
        DTPoint2D currentPoint = currentPoints(i);
        listOfSegments(i).AddPoint(currentPoint, 2);
    }
    //The start time is hardcoded for the video and will need to be intelligently modified later
    
    DTProgress progress;
    for (timeN=3;timeN<howManyTimes;timeN++) {
        current = allPoints.Get(timeValues(timeN));
        progress.UpdatePercentage(timeN/double(howManyTimes));
        
        if (current.IsEmpty()) continue; //If its empty continue
        howManySegments = listOfSegments.Length();
        
        // Computes distance matrix
        // dist(i,j) = distance from point i in the current point collection to point j in the segment list
        for (int j = 0; j < howManySegments; j++) {
            DTPoint2D endPoint = listOfSegments(j).lastPoint;
            for (int k = 0; k < current.NumberOfPoints(); k++) {
                DTPoint2D currentPoint = current(k);
                dist(j,k) = Norm(currentPoint-endPoint);
            }
        }
        
        cerr << timeN << endl;
        
        // Finds j,k such that dist(j,k) is the minimum distance
        while (1) {
            ssize_t minAt;
            double minDist = Minimum(dist,minAt); // minAt = I + J*dist.m()
            if (minDist>=maxDist) { //When true everything is matched up
                break;
            }
            int j = minAt%dist.m();
            int k = minAt/dist.m();
            
            listOfSegments(j).AddPoint(current(k), timeN);
            for (int i=0;i < dist.n(); i++) {
                dist(j,i) = maxDist;
            }
            for (int i=0;i < dist.m(); i++) {
                dist(i,k) = maxDist;
            }
        }
    }
    
    //Goes over the segment paths and stiches them togther using a function of the structure
    DTMutableList<DTPath2DValues> segmentPaths(howManySegments);
    for (int i=0;i<howManySegments;i++) {
        segmentPaths(i) = listOfSegments(i).ConvertToPath(timeValues);
    }
    return Combine(segmentPaths);
}

