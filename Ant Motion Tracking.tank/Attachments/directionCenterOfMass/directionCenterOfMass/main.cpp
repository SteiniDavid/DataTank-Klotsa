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
    
    DTPath2D ConvertToPath(){
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

    void pall(void) const;
};

void SingleList::pall(void) const
{
    cerr << "time = " << timeIndex << endl;
    ExtractColumns(points,DTRange(0,howManyPoints)).pall();
}


DTPath2D Computation(const DTSeriesPointCollection2D &allPoints,double maxDist)
{
    
    DTMutableList<SingleList> listOfSegments(10);
    int howManySegments = 0;
    
    DTDoubleArray timeValues = allPoints.TimeValues();
    int timeN, howManyTimes = timeValues.Length();
    DTPointCollection2D current;
    int frameNumber = 0;
    DTMutableDoubleArray dist(100,100);
    dist = maxDist;
    
    //Starting
    int numPoints = allPoints.Get(timeValues(2)).NumberOfPoints();
    for (int i = 0; i < numPoints; i++) {
        DTPointCollection2D currentPoints = allPoints.Get(timeValues(2));
        DTPoint2D currentPoint = currentPoints(i);
        
        listOfSegments(i).AddPoint(currentPoint, 2);
    }
    
    for (timeN=3;timeN<howManyTimes;timeN++) {
        current = allPoints.Get(timeValues(timeN));
        
        // For each point in current, find the closest point in the listOfSegments
        if (current.IsEmpty()) continue;
        howManySegments = listOfSegments.Length();
        // For each point in the segments, find a close match in the current

        // Compute a distance matrix
        // dist(i,j) = distance from point i in the current point collection to point j in the segment list
        for (int j = 0; j < howManySegments; j++) {
            DTPoint2D endPoint = listOfSegments(j).lastPoint; //   DTPoint2D(listOfSegments(j).points(0,j),listOfSegments(j).points(1,j));
            for (int k = 0; k < current.NumberOfPoints(); k++) {
                DTPoint2D currentPoint = current(k);
                // dist(j,k) = Norm(DTPoint2D(currentPoint.x-listOfSegments(j).points(0,j)),(currentPoint.y-listOfSegments(j).points(1,j)));
                dist(j,k) = Norm(currentPoint-endPoint);
                // cerr << j << "," << k << " : " << dist(j,k) << endl;
            }
        }
        
        // Find I,J such that dist(I,J) is the minimum distance
        
        while (1) {
            ssize_t minAt;
            double minDist = Minimum(dist,minAt); // minAt = I + J*dist.m()
            if (minDist>=maxDist) {
                // Done
                //DTErrorMessage("not finished");
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
       
        
        
        // Connect points I and J, i.e. point I in the current should be added to point J in the segment list
        // and then set dist(I,:) and dist(:,J) to the largest value (maxDist) so that it won't be picked again.
        // Repeat the last 3 lines until the minimum distance is >= maxDist
        
        
        // Any leftover points need to be added to the segments
        
       
    }
    
    // At this point go over the segments, tie them together into a path.
    DTMutableList<DTPath2D> segmentPaths(howManySegments);
    int i;
    for (i=0;i<howManySegments;i++) {
        segmentPaths(i) = listOfSegments(i).ConvertToPath();
    }
    return Combine(segmentPaths);
}
