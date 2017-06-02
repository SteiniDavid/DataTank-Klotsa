// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTProgress.h"
#include "DTRegion2D.h"
#include "DTSeriesVectorCollection2D.h"
#include "DTVectorCollection2D.h"

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTRegion2D &box,const DTVectorCollection2D &v_0,
                 double nu,double r,double maxtime,
                 DTSeriesVectorCollection2D &computed);

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
    double maxtime = inputFile.ReadNumber("maxtime");
    
    // Output series.
    DTSeriesVectorCollection2D computed(outputFile,"Var");
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
                 double nu,double r,double maxtime,
                 DTSeriesVectorCollection2D &computed)
{
    // Insert your code here.
    
    DTProgress progress;

    computed.Add(v_0, 0); //add the inital velocity vectors to start of the array of vector collections
    
    DTMutableDoubleArray x = v_0.Points().Data().Copy(); //extract out the points, data, and then copy it so that it does not reference the same object, whihc is declared as constant and cant be changed
    DTMutableDoubleArray vectors = v_0.Vectors().Copy(); //the vectors can also be extracted out of it and copied because you have to for it to work.
    DTMutableDoubleArray tempvectors = vectors; //copy of vectors created, notice that these are DTMutableDoubleArray objects, it means what it looks like
    
    double v = sqrt(vectors(0,0)*vectors(0,0)+vectors(1,0)*vectors(1,0)); //finds velocity of first vector
    int N = x.n(); // v_0.Points().NumberOfPoints(); //finds the number of vectors
    double dx,dy; //initializing dx and dy, used later
    double L = box.xmax; //gives you the max x value for the region passed in in the beginning, used for bounding
    DTRandom R(676757); //creates random object based on set seed, arbitrary
    
    for (int timeIteration = 0; timeIteration<maxtime; timeIteration++) { //loops for time until time is up
        progress.UpdatePercentage(timeIteration/double(maxtime)); //updates what the progress is, not really necessary but nice
        for (int i = 0; i< N; i++) { //for number of vectors in a given frame
            int counter = 0; //set counter to zero
            double sumx = 0; //set the sum of x components to zero
            double sumy = 0; //set the sum of y components to zero
            for (int j = 0; j< N; j++) { //again goes for ever particle in array, to do comparison with each to each
                dx = x(0,i)-x(0,j); //finds the change in the x
                //this next stuff is so that you can have a boundless box, its like inifinity but not
                if(dx > L/2) { //if the change is greater than half the length
                    dx -= L; //then subtract the length from it
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
                double d = sqrt(dx*dx+dy*dy); //finds the distance although it could have just been the norm
                if (d <= r) { //compares it to the passed in radius, checks if its less than or equal to the radius, it checks to see if the two particles are close enough
                    counter++; //increments the counter
                    sumx += vectors(0,j); //sums up x components of them both
                    sumy += vectors(1,j); //sums up y components of them both
                }
            }
            
            double vx = sumx/counter; //this gives you the average vx of the swarm that clusters
            double vy = sumy/counter; //this gives you the average vy of the swarm that clusters
            
            double theta = atan2(vy, vx); //creates a theta value or angle given the x and y component, using math function
            theta = theta + (R.UniformHalf()-.5)*nu; //adds noise to the theta value
            tempvectors(0,i) = cos(theta)*v; //saves the x component to the ith particle
            tempvectors(1,i) = sin(theta)*v; //saves the y component to the ith particle
        }
        for (int k= 0; k< N; k++) { //for every particle
            x(0,k) = x(0,k) + tempvectors(0,k); //adds to the xth part of the kth point the x value of the temp
            if (x(0,k) > L){ //handles boumding to give you an infinite box.
                x(0,k) -= L;
            }
            if (x(0,k) < 0) {
                x(0,k) += L;
            }
            
            x(1,k) = x(1,k) + tempvectors(1,k); //adds to the yth part of the kth point the y value of the temp
            if (x(1,k) > L){ //handles bounding
                x(1,k) -= L;
            }
            if (x(1,k) < 0) {
                x(1,k) += L;
            }
        }
        std::swap(tempvectors,vectors); //swaps the temp and vectors with each other
        
        computed.Add(DTVectorCollection2D(DTPointCollection2D(x),vectors), timeIteration+1); //adds to the computed a created DTPointCollection2D object created with the x values
                                                                                             //then that is combined with the vector DTMutableArray to create the DTVector colleciton object.
    }
    
}
