// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTProgress.h"
#include "DTRegion2D.h"
#include "DTSeriesGroup.h"
#include "DTVectorCollection2D.h"

//////////////////////////////////////////////////////////////////////////////
//    DT_RetGroup
//////////////////////////////////////////////////////////////////////////////

struct DT_RetGroup {
    DTVectorCollection2D vectors;
    
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
    cerr << pad << "vectors = "; vectors.pinfo();
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save("vectors",name+"_1N");
    output.Save("VectorCollection2D",name+"_1T");
    
    output.Save(1,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
    Write(output,name+"_vectors",var.vectors);
    Write(output,name,DTDoubleArray()); // So that DataTank can see the variable.
}

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTRegion2D &box,const DTVectorCollection2D &v_0,
                 double nu,double r,int maxtime,
                 DTSeriesGroup<DT_RetGroup> &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTRegion2D box;
    Read(inputFile,"box",box);
    DTVectorCollection2D v_0;
    Read(inputFile,"v_0",v_0);
    double nu = inputFile.ReadNumber("nu");
    double r = inputFile.ReadNumber("r");
    int maxtime = int(inputFile.ReadNumber("maxtime"));
    
    // Output series.
    DTSeriesGroup<DT_RetGroup> computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"box",box);
        WriteOne(outputFile,"v_0",v_0);
        WriteOne(outputFile,"nu",nu);
        WriteOne(outputFile,"r",r);
        WriteOne(outputFile,"maxtime",maxtime);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(box,v_0,nu,r,maxtime,computed);
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

void Computation(const DTRegion2D &box,const DTVectorCollection2D &v_0,
                 double nu,double r,int maxtime,
                 DTSeriesGroup<DT_RetGroup> &computed)
{
    // Insert your code here.
    
    DTProgress progress;
    DT_RetGroup returnStructure;
    
    returnStructure.vectors = v_0;
    computed.Add(returnStructure, 0.0);
    
    DTMutableDoubleArray x = v_0.Points().Data().Copy();
    DTMutableDoubleArray vectors = v_0.Vectors().Copy();
    DTMutableDoubleArray tempvectors = vectors;
    
    double v = sqrt(vectors(0,0)*vectors(0,0)+vectors(1,0)*vectors(1,0));
    int N = x.n(); // v_0.Points().NumberOfPoints();
    double dx,dy;
    double L = box.xmax;
    DTRandom R(676757);

    for (int timeIteration = 0; timeIteration<maxtime; timeIteration++) {
        progress.UpdatePercentage(timeIteration/double(maxtime));
        for (int i = 0; i< N; i++) {
            int counter = 0;
            double sumx = 0;
            double sumy = 0;
            for (int j = 0; j< N; j++) {
                dx = x(0,i)-x(0,j);
                if(dx > L/2) {
                    dx -= L;
                }
                if (dx < -L/2) {
                    dx += L;
                }
                
                dy = x(1,i)-x(1,j);
                if(dy > L/2) {
                    dy -= L;
                }
                if (dy < -L/2) {
                    dy += L;
                }
                double d = sqrt(dx*dx+dy*dy);
                if (d <= r) {
                    counter++;
                    sumx += vectors(0,j);
                    sumy += vectors(1,j);
                }
            }
            
            double vx = sumx/counter;
            double vy = sumy/counter;
            
            double theta = atan2(vy, vx);
            theta = theta + (R.UniformHalf()-.5)*nu;
            tempvectors(0,i) = cos(theta)*v;
            tempvectors(1,i) = sin(theta)*v;
        }
        for (int k= 0; k< N; k++) {
            x(0,k) = x(0,k) + tempvectors(0,k);
            if (x(0,k) > L){
                x(0,k) -= L;
            }
            if (x(0,k) < 0) {
                x(0,k) += L;
            }

            x(1,k) = x(1,k) + tempvectors(1,k);
            if (x(1,k) > L){
                x(1,k) -= L;
            }
            if (x(1,k) < 0) {
                x(1,k) += L;
            }
        }
        std::swap(tempvectors,vectors);
        
        returnStructure.vectors = DTVectorCollection2D(DTPointCollection2D(x),vectors);
        computed.Add(returnStructure, timeIteration+1);
    }
    
}