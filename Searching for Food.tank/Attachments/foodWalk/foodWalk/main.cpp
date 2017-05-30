// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTMesh2D.h"
#include "DTPoint2D.h"
#include "DTPointCollection2D.h"
#include "DTProgress.h"
#include "DTRegion2D.h"
#include "DTSeriesGroup.h"

//////////////////////////////////////////////////////////////////////////////
//    DT_RetGroup
//////////////////////////////////////////////////////////////////////////////

struct DT_RetGroup {
    DTPointCollection2D position;
    DTMesh2D foodBeacon;
    DTMesh2D homeBeacon;
    
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
    cerr << pad << "position = "; position.pinfo();
    cerr << pad << "foodBeacon = "; foodBeacon.pinfo();
    cerr << pad << "homeBeacon = "; homeBeacon.pinfo();
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save("position",name+"_1N");
    output.Save("PointCollection2D",name+"_1T");
    
    output.Save("foodBeacon",name+"_2N");
    output.Save("Mesh2D",name+"_2T");
    
    output.Save("homeBeacon",name+"_3N");
    output.Save("Mesh2D",name+"_3T");
    
    output.Save(3,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
    Write(output,name+"_position",var.position);
    Write(output,name+"_foodBeacon",var.foodBeacon);
    Write(output,name+"_homeBeacon",var.homeBeacon);
    Write(output,name,DTDoubleArray()); // So that DataTank can see the variable.
}

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTPointCollection2D &initial,const DTPoint2D &food,
                 double dt,double endTime,int seed,
                 const DTDictionary &parameters,const DTRegion2D &region,
                 DTSeriesGroup<DT_RetGroup> &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTPointCollection2D initial;
    Read(inputFile,"initial",initial);
    DTPoint2D food;
    Read(inputFile,"food",food);
    double dt = inputFile.ReadNumber("dt");
    double endTime = inputFile.ReadNumber("endTime");
    int seed = int(inputFile.ReadNumber("seed"));
    DTDictionary parameters;
    Read(inputFile,"parameters",parameters);
    DTRegion2D region;
    Read(inputFile,"region",region);
    
    // Output series.
    DTSeriesGroup<DT_RetGroup> computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"initial",initial);
        WriteOne(outputFile,"food",food);
        WriteOne(outputFile,"dt",dt);
        WriteOne(outputFile,"endTime",endTime);
        WriteOne(outputFile,"seed",seed);
        WriteOne(outputFile,"parameters",parameters);
        WriteOne(outputFile,"region",region);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(initial,food,dt,endTime,seed,parameters,region,computed);
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

#include "DTRandom.h"
#include "DTPointValueCollection2D.h"

void createHomeBeacon(const DTPointCollection2D &val, const DTDictionary &parameters);

void Computation(const DTPointCollection2D &initial,const DTPoint2D &food,
                 double dt,double endTime,int seed,
                 const DTDictionary &parameters,const DTRegion2D &region,
                 DTSeriesGroup<DT_RetGroup> &computed)
{
    DTMutablePointCollection2D ants = initial.Copy();
    int nunAnts = ants.NumberOfPoints();
    
    DTMutableDoubleArray antArray = ants.Data();
    double x,y;
    
    DT_RetGroup state;
    DTRandom rndNumber(seed);
    double xy[2];
    
    double sqrtDT = sqrt(dt);
    
    double noise = parameters("noise");
    double radius = parameters("radius");
    
    //Reading in the parameters of region
    double xMax = region.xmax;
    double yMax = region.ymax;
    double xMin = region.xmin;
    double yMin = region.ymin;
    
    double t = 0;
    while (t<endTime) {
        
        for (int i = 0; i<nunAnts; i++) {
            x = antArray(0,i);
            y = antArray(1,i);
            
            
            //Random Walk
            rndNumber.Normal(xy,2);
            x += sqrtDT*noise*xy[0];
            y += sqrtDT*noise*xy[1];
            
            //Hitting the boundary has to be handled
            //It will be treated as a ballistic reflection
            if (x > xMax) {
                x = 2*xMax-x;
            }
            if (y > yMax) {
                y = 2*yMax-y;
            }
            if (x < xMin) {
                x = 2*xMin-x;
            }
            if (y < yMin) {
                y = 2*yMin-y;
            }
            
            //once values have been updated
            antArray(0,i) = x;
            antArray(1,i) = y;
        }
        
        t += dt;
        
        // Save
        state.position = ants;
        state.foodBeacon;
        state.homeBeacon;
        computed.Add(state,t);
    }
    
    DTPointCollection2D temp = initial;
    // createHomeBeacon(temp);
    // tmep is now empty
    
    DTProgress progress;
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
}

void createHomeBeacon(const DTPointCollection2D &initial, const DTDictionary &parameters){
    DTMutableIntArray swarmID(initial.NumberOfPoints());
    DTMutableIntArray swarmBiggestID(initial.NumberOfPoints());
    double radius = parameters("radius");
    DTDoubleArray antArray = initial.Data();
    int numPoints = initial.NumberOfPoints();

    for (int i = 0; i < numPoints; i++) {
        
    }
    
}





