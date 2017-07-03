// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointCollection2D.h"
#include "DTProgress.h"
#include "DTRegion2D.h"
#include "DTSeriesPointCollection2D.h"
//My imports
#include "DTRandom.h"
#include <cmath>

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTPointCollection2D &inital,const DTRegion2D &region,
                 const DTDictionary &parameters,
                 DTSeriesPointCollection2D &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTPointCollection2D inital;
    Read(inputFile,"inital",inital);
    DTRegion2D region;
    Read(inputFile,"region",region);
    DTDictionary parameters;
    Read(inputFile,"parameters",parameters);
    
    // Output series.
    DTSeriesPointCollection2D computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"inital",inital);
        WriteOne(outputFile,"region",region);
        WriteOne(outputFile,"parameters",parameters);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(inital,region,parameters,computed);
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
//    My Objects
//////////////////////////////////////////////////////////////////////////////
class Particle {
public:
    Particle() : x(0), y(0){};
    Particle(DTPoint2D singleParticle) : x(singleParticle.x), y(singleParticle.y) {};
    
    //Instance variables
    double x,y;
    //double xBef, yBef;
    double vx,vy;
    double ax, ay;
};

//////////////////////////////////////////////////////////////////////////////
//    My Function Prototypes
//////////////////////////////////////////////////////////////////////////////

DTMutableList<Particle> instantiate(DTPointCollection2D inital, int numParticles, DTRandom r, double velocity);
DTMutableDoubleArray makePointCollection(DTMutableList<Particle> particles);
void boundCheck(DTMutableList<Particle> &particles, DTRegion2D region);

//////////////////////////////////////////////////////////////////////////////
//    My Functions
//////////////////////////////////////////////////////////////////////////////

DTMutableList<Particle> instantiate(DTPointCollection2D inital, int numParticles, DTRandom r, double velocity) {
    DTMutableList<Particle> particles(numParticles);
    for (int i = 0; i < numParticles; i++) {
        DTPoint2D singleParticle = inital(i);
        particles(i) = Particle(singleParticle);
        double angle = r.UniformHalf()*2*M_PI;
        particles(i).vx = cos(angle)*velocity; //initial vx
        particles(i).vy = sin(angle)*velocity; //inital vy
    }
    return particles;
}

DTMutableDoubleArray makePointCollection(DTMutableList<Particle> particles) {
    int length = particles.Length();
    DTMutableDoubleArray points(length,2);
    for (int i = 0; i < length; i++) {
        points(i,0) = particles(i).x;
        points(i,1) = (particles(i).y);
    }
    return points;
}

void boundCheck(DTMutableList<Particle> &particles, DTRegion2D region) {
    double L = region.xmax; //only need one L because its a box
    for (int i = 0; i < particles.Length(); i++) {
        if (particles(i).x > L){ //handles boumding to give you an infinite box.
            particles(i).x -= L;
        }
        if (particles(i).x < 0) {
            particles(i).x += L;
        }
        if (particles(i).y > L){ //handles bounding
            particles(i).y -= L;
        }
        if (particles(i).y < 0) {
            particles(i).y += L;
        }
    }
}

void calculateAccelerations(DTMutableList<Particle> &particles, DTRegion2D region){
    int length = particles.Length();
    double L = region.xmax;
    
    for (int i = 0; i < length; i++) {
        particles(i).ax = 0;
        particles(i).ax = 0;
    }
    
    for (int i = 0; i < length-1; i++) {
        for (int j = i+1; j < length; j++) {
            double distXY[2];
            double dist = 0;
            
            distXY[0] = particles(i).x-particles(j).x;
            distXY[1] = particles(i).y-particles(j).y;
            for (int k = 0; k < 2; k++) {
                if (abs(distXY[k]) > 0.5 * L) {
                    if (distXY[k] > 0) {
                        distXY[k] -= L;
                    }
                    else {
                        distXY[k] += L;
                    }
                }
                dist += distXY[k]*distXY[k];
            }
            double f = 24 * (2 * pow(dist, -7) - pow(dist, -4));
            particles(i).x += distXY[0]*f;
            particles(i).y += distXY[1]*f;
            
            particles(j).x -= distXY[0]*f;
            particles(j).y -= distXY[1]*f;
        }
    }
}

void computeAccelerations() {
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            a[i][k] = 0; //all of the accelerations are set to zero to re wipe them
    for (int i = 0; i < N-1; i++) //goes over all particles except the last because it looks at pairs
        for (int j = i+1; j < N; j++) { //goes through to create pairs
            double rij[3]; //distnaces in the three dims
            double rSqd = 0; //set as zero the distance in a given dim
            for (int k = 0; k < 3; k++) { //for the three dim
                rij[k] = r[i][k] - r[j][k]; //finds the diff for a dim
                // closest image convention
                if (abs(rij[k]) > 0.5 * L) { //since its periodic particles on edges interact with opposite side
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k]; //adds the magnitdue of the three dim for a given particle
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4)); //calculates the lenard jones force for a given particle pair
            for (int k = 0; k < 3; k++) { //for the three dim
                a[i][k] += rij[k] * f; //the distances scaled by the force is added as the acceleration for the ith particle in the kth dim
                a[j][k] -= rij[k] * f; //the negation of the previous is done for the jth particle
            }
        }
}

void verlet(DTMutableList<Particle> &particles, DTRegion2D region, double dt) {
    int length = particles.Length();
    double L = region.xmax;
    
    calculateAccelerations(particles, region);
    
    for (int i = 0; i < length; i++) {
        
    }
}

void velocityVerlet(double dt) {
    computeAccelerations(); //accerations are computed for the verlet algorithm to do translations
    for (int i = 0; i < N; i++) //for N particles
        for (int k = 0; k < 3; k++) { //for the three dim
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt; //handles the displacement based on the current velocity and the current acceleration
            // use periodic boundary conditions
            if (r[i][k] < 0)
                r[i][k] += L;
            if (r[i][k] >= L)
                r[i][k] -= L;
            v[i][k] += 0.5 * a[i][k] * dt; //changes the velocity based on the acceleration
        }
    computeAccelerations(); //accelerations are calculated again for the next time step
    for (int i = 0; i < N; i++) //goes over all particles
        for (int k = 0; k < 3; k++) //goes over all dim
            v[i][k] += 0.5 * a[i][k] * dt; //creates new velocities for all particles (why is it done twice?)
}


//////////////////////////////////////////////////////////////////////////////
//    Computational routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTPointCollection2D &inital,const DTRegion2D &region,
                 const DTDictionary &parameters,
                 DTSeriesPointCollection2D &computed)
{
    computed.Add(inital, 0.0);
    
    //Read out the parameters variables
    double velocity = parameters("velocity");
    double radius = parameters("radius");
    double maxTime = parameters("maxTime");
    double nu = parameters("nu");
    double seed = parameters("seed");
    
    int numParticles = inital.NumberOfPoints();
    DTRandom r(seed);
    
    DTMutableList<Particle> particles = instantiate(inital, numParticles, r, velocity);
    
    
    for (int i = 0; i < maxTime; i++) {
        
        
        boundCheck(particles,region);
        DTMutableDoubleArray pointArray = makePointCollection(particles);
        DTPointCollection2D points = DTPointCollection2D(pointArray);
        computed.Add(points, i+1);
    }
    
    cout << particles(0).y << '\n';

    
    
    DTProgress progress;
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
}
