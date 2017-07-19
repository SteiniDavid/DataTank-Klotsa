// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointCollection2D.h"
#include "DTPointValueCollection2D.h"
#include "DTProgress.h"
#include "DTSeriesGroup.h"
#include "DTSeriesPointCollection2D.h"

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
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTSeriesPointCollection2D &points,double radius,
                 DTSeriesGroup<DT_RetGroup> &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTSeriesPointCollection2D points;
    Read(inputFile,"points",points);
    double radius = inputFile.ReadNumber("radius");

    // Output series.
    DTSeriesGroup<DT_RetGroup> computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"radius",radius);
    }


    // The computation.
    clock_t t_before = clock();
    Computation(points,radius,computed);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    outputFile.SaveIndex();

    return 0;
}

//////////////////////////////////////////////////////////////////////////////
//    Computational routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTSeriesPointCollection2D &points,double radius,
                 DTSeriesGroup<DT_RetGroup> &computed)
{
    DTDoubleArray timeValues = points.TimeValues();
    int howManyTimes = timeValues.Length();
    int numPoints = points.Get(1).NumberOfPoints();
    
    for (int time = 0; time < howManyTimes; time++) {
        DT_RetGroup state;
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
        
        computed.Add(state, time);
    }

    
    
    
    
    

    DTProgress progress;

    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.

}
