// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDoubleArray.h"
#include "DTPointCollection2D.h"
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
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save(0,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
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
    
    int numPoints = points.Get(timeValues(1)).NumberOfPoints();
    DTMutableDoubleArray values(numPoints);
    for (int i = 0; i < numPoints; i++) {
        values(i) = 0;
    }
    for (int timeStep = 0; timeStep < howManyTimes; timeStep++) {
        
        DTPointCollection2D currentPoints = points.Get(timeValues(2));
        for (int i = 1; i < numPoints; i++) {
            for (int j = 0; j < numPoints-1; j++) {
                double dist = Norm(currentPoints(i) - currentPoints(j));
                if (dist <= radius*radius && i != j) {
                    cout << "match" << endl;
                    values(i)=1;
                    values(j)=1;
                }
                cout << currentPoints(i).x << " , " << currentPoints(j).x << " , " << dist << endl;
            }
        }
        
    }
    
    
    DT_RetGroup toReturn;


    return toReturn;
}
