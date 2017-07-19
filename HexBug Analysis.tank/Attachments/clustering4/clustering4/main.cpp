// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"
#include "DTPointCollection2D.h"
#include "DTPointValueCollection2D.h"
#include "DTSeriesPointCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

//////////////////////////////////////////////////////////////////////////////
//    DT_RetGroup
//////////////////////////////////////////////////////////////////////////////

struct DT_RetGroup {
    DTPointValueCollection2D cluster;
    double avgClusterSize;
    double percentInCluster;

    void pinfo(void) const;
    void pinfoIndent(string) const;

    static void WriteStructure(DTDataStorage &,string);
};

void DT_RetGroup::pinfo(void) const
{
    pinfoIndent("");
}

void DT_RetGroup::pinfoIndent(string pad) const
{
    cerr << pad << "cluster = "; cluster.pinfo();
    cerr << pad << "avgClusterSize = " << avgClusterSize << endl;
    cerr << pad << "percentInCluster = " << percentInCluster << endl;
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save("cluster",name+"_1N");
    output.Save("PointValueCollection2D",name+"_1T");

    output.Save("avgClusterSize",name+"_2N");
    output.Save("Real Number",name+"_2T");

    output.Save("percentInCluster",name+"_3N");
    output.Save("Real Number",name+"_3T");

    output.Save(3,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
    Write(output,name+"_cluster",var.cluster);
    output.Save(var.avgClusterSize,name+"_avgClusterSize");
    output.Save(var.percentInCluster,name+"_percentInCluster");
    Write(output,name,DTDoubleArray()); // So that DataTank can see the variable.
}

//////////////////////////////////////////////////////////////////////////////
//    End of structure definitions.
//////////////////////////////////////////////////////////////////////////////

DT_RetGroup Computation(const DTSeriesPointCollection2D &points,double radius);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTSeriesPointCollection2D points;
    Read(inputFile,"points",points);
    double radius = inputFile.ReadNumber("radius");

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    DT_RetGroup computed;

    // The computation.
    clock_t t_before = clock();
    computed = Computation(points,radius);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    outputFile.Save("Group","Seq_Var");
    DT_RetGroup::WriteStructure(outputFile,"SeqInfo_Var");
    Write(outputFile,"Var",computed);

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    outputFile.SaveIndex();

    return 0;
}

DT_RetGroup Computation(const DTSeriesPointCollection2D &points,double radius)
{
    DTDoubleArray timeValues = points.TimeValues();
    int howManyTimes = timeValues.Length();
    int numPoints = points.Get(1).NumberOfPoints();
    
    DT_RetGroup state;
    
    for (int time = 0; time < howManyTimes; time++) {
        DTPointCollection2D currentPts = points.Get(time);
        
        DTMutableDoubleArray values(numPoints);
        
        for (int i = 0; i < numPoints; i++) {
            values(i) = i;
        }
        
        for (int i = 1; i < numPoints-1; i++) {
            for (int j = 0; j < numPoints; j++) {
                double dist = Norm(currentPts(i) - currentPts(j));
                if (dist <= radius*radius && i != j) {
                    //cout << "match" << endl;
                    if (values(i) < values(j)) {
                        values(j)=values(i);
                    } else  {
                        values(i)=values(j);
                    }
                }
                //cout << points(i).x << " , " << points(j).x << " , " << dist << endl;
            }
        }
        
        for (int i = 0; i < values.Length(); i++) {
            int sumI = 0;
            for (int j = 0; j < values.Length(); j++) {
                if (values(i) == values(j)) {
                    sumI++;
                }
            }
            if (sumI < 3) {
                values(i) = -1;
                //cout << "set to -1" << endl;
            }
            
        }
        //calculating average cluster size
        double numInCluster = 0;
        double numClusters = 0;
        for (int i = 0; i < values.Length(); i++) {
            int sumI = 0;
            for (int j = 0; j < values.Length(); j++) {
                if (i == values(j)) {
                    numInCluster++;
                    sumI++;
                }
            }
            if (sumI > 2) {
                numClusters++;
            }
        }
        
        state.cluster = DTPointValueCollection2D(currentPts, values);
        state.avgClusterSize = numInCluster/numClusters;;
        state.percentInCluster = numInCluster/numPoints;
        
    }
    
    // toReturn.cluster - 2D Point Value Collection
    // toReturn.avgClusterSize - Real Number
    // toReturn.percentInCluster - Real Number

    return state;
}
