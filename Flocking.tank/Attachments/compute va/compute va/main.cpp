// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTPlot1D.h"
#include "DTRegion2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>

DTPlot1D Computation(const DTRegion2D &box,const DTVectorCollection2D &v_0,double nu,
                     double r,int maxtime);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTRegion2D box;
    Read(inputFile,"box",box);
    DTVectorCollection2D v_0;
    Read(inputFile,"v_0",v_0);
    double nu = inputFile.ReadNumber("nu");
    double r = inputFile.ReadNumber("r");
    int maxtime = int(inputFile.ReadNumber("maxtime"));

    // The computation.
    DTPlot1D computed;
    clock_t t_before = clock();
    computed = Computation(box,v_0,nu,r,maxtime);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("Plot1D","Seq_Var");

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    return 0;
}

#include "DTRandom.h"

double ComputeVA(const DTVectorCollection2D &v);

DTPlot1D Computation(const DTRegion2D &box,const DTVectorCollection2D &v_0,double nu,
                     double r,int maxtime)
{
    DTProgress progress;
    
    DTMutableDoubleArray x = v_0.Points().Data().Copy();
    DTMutableDoubleArray vectors = v_0.Vectors().Copy();
    DTMutableDoubleArray tempvectors = vectors;
    
    double v = sqrt(vectors(0,0)*vectors(0,0)+vectors(1,0)*vectors(1,0));
    int N = x.n(); // v_0.Points().NumberOfPoints();
    double dx,dy;
    double L = box.xmax;
    double L2 = L/2;
    DTRandom R(676757);
    double Lr = L-r;
    
    double *xD = x.Pointer();
    
    DTMutableDoubleArray xValues(maxtime+1);
    DTMutableDoubleArray yValues(maxtime+1);
    
    xValues(0) = 0;
    yValues(0) = ComputeVA(v_0);

    for (int timeIteration = 0; timeIteration<maxtime; timeIteration++) {
        progress.UpdatePercentage(timeIteration/double(maxtime));
        for (int i = 0; i< N; i++) {
            int counter = 0;
            double sumx = 0;
            double sumy = 0;
            double xi = x(0,i);
            double yi = x(1,i);
            double xj;
            double xipr = xi+r;
            double ximr = xi-r;
            
            if (xi<r || xi>Lr) {
                for (int j = 0; j< N; j++) {
                    dx = xi - xD[2*j]; // dx = xi-x(0,j);
                    if(dx > L2) {
                        dx -= L;
                    }
                    if (dx < -L2) {
                        dx += L;
                    }
                    if (dx > r) continue;
                    if (dx < -r) continue;
                    
                    dy = yi-xD[2*j+1]; // x(1,i)-x(1,j);
                    if(dy > L2) {
                        dy -= L;
                    }
                    if (dy < -L2) {
                        dy += L;
                    }
                    if (dy > r) continue;
                    if (dy < -r) continue;
                    
                    double d = sqrt(dx*dx+dy*dy);
                    if (d <= r) {
                        counter++;
                        sumx += vectors(0,j);
                        sumy += vectors(1,j);
                    }
                }
            }
            else {
                for (int j = 0; j< N; j++) {
                    xj = xD[2*j];
                    if (xj>xipr || xj<ximr) continue;
                    
                    dx = xi - xj; // xi-x(0,j);
                    // if (dx<-r || dx>r) continue;
                    
                    dy = yi-xD[2*j+1]; // x(1,i)-x(1,j);
                    if(dy > L2) {
                        dy -= L;
                    }
                    if (dy < -L2) {
                        dy += L;
                    }
                    if (dy > r) continue;
                    if (dy < -r) continue;
                    
                    double d = sqrt(dx*dx+dy*dy);
                    if (d <= r) {
                        counter++;
                        sumx += vectors(0,j);
                        sumy += vectors(1,j);
                    }
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
        
        DTVectorCollection2D result = DTVectorCollection2D(DTPointCollection2D(x),vectors);
        
        xValues(timeIteration+1) = timeIteration+1;
        yValues(timeIteration+1) = ComputeVA(result);
    }
    
    return DTPlot1D(xValues,yValues);
}

double ComputeVA(const DTVectorCollection2D &v)
{
    double toReturn;
    
    // v.Points() - DTPointCollection2D - underlying points
    // v.Points().Data() - DTDoubleArray - A 2xN array
    // v.Points().NumberOfPoints() - int - N
    // v.Vectors() - DTDoubleArray - A 2xN array
    
    DTDoubleArray vectors = v.Vectors();
    int N = v.Points().NumberOfPoints();
    double vmag = sqrt(vectors(0,0)*vectors(0,0)+vectors(1,0)*vectors(1,0));
    
    double v_a_frontpart = (1/(N*vmag));
    double v_a_x = 0;
    double v_a_y = 0;
    for (int i = 0; i < N; i++) {
        v_a_x += vectors(0,i);
        v_a_y += vectors(1,i);
    }
    toReturn = v_a_frontpart*sqrt(v_a_x*v_a_x + v_a_y*v_a_y);
    
    return toReturn;
}

