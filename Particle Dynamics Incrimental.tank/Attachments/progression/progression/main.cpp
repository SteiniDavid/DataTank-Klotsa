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
    Particle() : x(NAN), y(NAN), vx(NAN), vy(NAN), ax(NAN), ay(NAN) {};
    Particle(DTPoint2D singleParticle) : x(singleParticle.x), y(singleParticle.y), vx(NAN), vy(NAN), ax(NAN), ay(NAN) {};
    
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
void calculateAccelerations(DTMutableList<Particle> &particles, DTRegion2D region);
void verlet(DTMutableList<Particle> &particles, DTRegion2D region, double dt);

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
    DTMutableDoubleArray points(2,length);
    for (int i = 0; i < length; i++) {
        points(0,i) = particles(i).x;
        points(1,i) = (particles(i).y);
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
        particles(i).ay = 0;
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
//            if (f > 100) {
//                f = 10;
//            }
            cout << particles(0).x << '\n';
            cout << particles(0).y << '\n';
            
            particles(i).ax += distXY[0]*f;
            particles(i).ay += distXY[1]*f;
            
            particles(j).ax -= distXY[0]*f;
            particles(j).ay -= distXY[1]*f;
            
            cout << particles(0).x << '\n';
            cout << particles(0).y << '\n';
        }
    }
}


void calculateAccelerationsLangevin(DTMutableList<Particle> &particles, DTRegion2D region){
    int length = particles.Length();
    double L = region.xmax;
    double f[length][2];
    double rc = .1;
    double rc2 = rc*rc;
    double en = 0;
    for (int i = 0; i < length; i++) {
        particles(i).ax = 0;
        particles(i).ay = 0;
        f[i][0] = 0;
        f[i][1] = 0;
    }
    
    for (int i = 0; i < length-1; i++) {
        for (int j = i+1; j < length; j++) {
            double distXY[2];
            double r2 = 0;
            
            distXY[0] = particles(i).x-particles(j).x;
            distXY[1] = particles(i).y-particles(j).y;
            
            distXY[0] = distXY[0]-L*round(distXY[0]/L);
            distXY[1] = distXY[1]-L*round(distXY[1]/L);
            
            r2 = distXY[0]*distXY[0]+distXY[1]*distXY[1];
            //dist = dist-L*round(dist/L);
            
            if (r2 < rc2) {
                double r2i = 1/r2;
                double r6i = pow(r2i, 3.0);
                double ff = 48*r2i*r6i*(r6i-0.5);
                f[i][0] = f[i][0]+ff*distXY[0];
                f[i][1] = f[i][1]+ff*distXY[1];
                
                f[j][0] = f[j][0]-ff*distXY[0];
                f[j][1] = f[j][1]-ff*distXY[1];
                double ecut = 4*((1/pow(rc, 12.0))-(1/pow(rc, 6.0)));
                en = en + 4*r6i*(r6i-1)-ecut;
            }
        }
    }
}

void verletLangevin(DTMutableList<Particle> &particles, DTRegion2D region, double dt) {
    int length = particles.Length();
    double L = region.xmax;
    
    double sumv = 0;
    double sumv2 = 0;
    for (int i; i < length; i++) {
        
    }
    
}


void verlet(DTMutableList<Particle> &particles, DTRegion2D region, double dt) {
    int length = particles.Length();
    double L = region.xmax;
    
    calculateAccelerations(particles, region);
    
    for (int i = 0; i < length; i++) {
        particles(i).x += particles(i).vx * dt + 0.5 * particles(i).ax * dt * dt;
        particles(i).y += particles(i).vy * dt + 0.5 * particles(i).ay * dt * dt;
        
        if (particles(i).x < 0) {
            particles(i).x += L;
        }
        if (particles(i).x >= L) {
            particles(i).x -= L;
        }
        if (particles(i).y < 0) {
            particles(i).y += L;
        }
        if (particles(i).y >= L) {
            particles(i).y -= L;
        }
        
        particles(i).vx += 0.5 * particles(i).ax * dt;
        particles(i).vy += 0.5 * particles(i).ay * dt;
    }
//    cout << particles(0).x << '\n'; //its a list of particle objects
//    cout << particles(0).y << '\n';
    
    calculateAccelerations(particles, region);
    
    for (int i = 0; i < length; i++) {
        particles(i).vx += 0.5 * particles(i).ax * dt;
        particles(i).vy += 0.5 * particles(i).ay * dt;
    }
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
    double dt = parameters("dt");
    
    int numParticles = inital.NumberOfPoints();
    DTRandom r(seed);
    
    DTMutableList<Particle> particles = instantiate(inital, numParticles, r, velocity);
    
    
    for (int i = 0; i < maxTime/dt; i++) {
        verlet(particles, region, dt);
        
        boundCheck(particles,region);
        DTMutableDoubleArray pointArray = makePointCollection(particles);
        DTPointCollection2D points = DTPointCollection2D(pointArray);
        computed.Add(points, i+1); //ive got another weirder problem, the value does not change.
    }
    
    cout << particles(0).y << '\n';

    
    
    DTProgress progress;
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
}
