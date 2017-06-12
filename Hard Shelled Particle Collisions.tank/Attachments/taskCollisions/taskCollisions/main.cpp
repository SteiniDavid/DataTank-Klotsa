// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointCollection3D.h"
#include "DTProgress.h"
#include "DTRegion3D.h"
#include "DTSeriesGroup.h"
#include "DTVectorCollection3D.h"

//////////////////////////////////////////////////////////////////////////////
//    DT_RetGroup
//////////////////////////////////////////////////////////////////////////////

struct DT_RetGroup {
    DTVectorCollection3D velocities;
    
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
    cerr << pad << "velocities = "; velocities.pinfo();
}

void DT_RetGroup::WriteStructure(DTDataStorage &output,string name)
{
    output.Save("velocities",name+"_1N");
    output.Save("VectorCollection3D",name+"_1T");
    
    output.Save(1,name+"_N");
    output.Save("Group",name);
}

extern void Write(DTDataStorage &,string name,const DT_RetGroup &);

void Write(DTDataStorage &output,string name,const DT_RetGroup &var)
{
    Write(output,name+"_velocities",var.velocities);
    Write(output,name,DTDoubleArray()); // So that DataTank can see the variable.
}

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(const DTDictionary &parameters,
                 const DTPointCollection3D &initialPoints,
                 const DTRegion3D &box,DTSeriesGroup<DT_RetGroup> &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    DTDictionary parameters;
    Read(inputFile,"parameters",parameters);
    DTPointCollection3D initialPoints;
    Read(inputFile,"initialPoints",initialPoints);
    DTRegion3D box;
    Read(inputFile,"box",box);
    
    // Output series.
    DTSeriesGroup<DT_RetGroup> computed(outputFile,"Var");
    if (DTArgumentIncludesFlag("saveInput")) { // Add -saveInput to the argument list to save the input in the output file.
        WriteOne(outputFile,"parameters",parameters);
        WriteOne(outputFile,"initialPoints",initialPoints);
        WriteOne(outputFile,"box",box);
    }
    
    
    // The computation.
    clock_t t_before = clock();
    Computation(parameters,initialPoints,box,computed);
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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

const int N = 64; // number of particles
double r[N][3]; // positions
double v[N][3]; // velocities
double a[N][3]; // accelerations

double L = 10;            // linear size of cubical volume
double vMax = 0.1;        // maximum initial velocity component

void initialize() {
    // initialize positions
    int n = int(ceil(pow(N, 1.0/3)));  // number of atoms in each direction
    double a = L / n;                  // lattice spacing
    int p = 0;                         // particles placed so far
    for (int x = 0; x < n; x++)		   //for the x components of all particles
        for (int y = 0; y < n; y++)    //for the y components of all particles
            for (int z = 0; z < n; z++) { //for the z components of all of the particles
                if (p < N) {  //if p is less than N
                    r[p][0] = (x + 0.5) * a; //position at p for the x component is the x comp plus .5 times the lattice spacing
                    r[p][1] = (y + 0.5) * a;
                    r[p][2] = (z + 0.5) * a;
                }
                ++p; //increment the position
            }
    // initialize velocities
				for (int p = 0; p < N; p++) //for loop of p
                    for (int i = 0; i < 3; i++) //the three didm
                        v[p][i] = vMax * (2 * rand() / double(RAND_MAX) - 1); //
}

void computeAccelerations() {
    for (int i = 0; i < N; i++)             // set all accelerations to zero
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;
    for (int i = 0; i < N-1; i++)           // loop over all distinct pairs i,j
        for (int j = i+1; j < N; j++) {
            double rij[3];                  // position of i relative to j
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            cout << f << '\n';
            for (int k = 0; k < 3; k++) {
                a[i][k] += rij[k] * f;
                a[j][k] -= rij[k] * f;
            }
        }
}

void velocityVerlet(double dt) {
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++) {
            r[i][k] += v[i][k] * dt + 0.5 * a[i][k] * dt * dt;
            v[i][k] += 0.5 * a[i][k] * dt;
        }
    computeAccelerations();
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            v[i][k] += 0.5 * a[i][k] * dt;
}

double instantaneousTemperature() {
    double sum = 0;
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            sum += v[i][k] * v[i][k];
    return sum / (3 * (N - 1));
}


void Computation(const DTDictionary &parameters,
                 const DTPointCollection3D &initialPoints,
                 const DTRegion3D &box,DTSeriesGroup<DT_RetGroup> &computed)
{
    initialize();
    DTMutableDoubleArray position(3,N);
    DTMutableDoubleArray velocity(3,N);
    for(int i = 0; i< N; i++) {
        position(0,i) = r[i][0];
        position(1,i) = r[i][1];
        position(2,i) = r[i][2];
        velocity(0,i) = v[i][0];
        velocity(1,i) = v[i][1];
        velocity(2,i) = v[i][2];
    }
    DTVectorCollection3D vectors(DTPointCollection3D(position), velocity);
    DT_RetGroup state;
    state.velocities = vectors;
    computed.Add(state, 0.0);
    double dt = parameters("dt");
    double maxTime = parameters("maxTime");
    for (int i = 0; i < maxTime/dt; i++) {
        velocityVerlet(dt);
        cout << instantaneousTemperature() << '\n';
        for(int i = 0; i< N; i++) {
            position(0,i) = r[i][0];
            position(1,i) = r[i][1];
            position(2,i) = r[i][2];
            velocity(0,i) = v[i][0];
            velocity(1,i) = v[i][1];
            velocity(2,i) = v[i][2];
        }
        DTVectorCollection3D vectors(DTPointCollection3D(position), velocity);
        state.velocities = vectors;
        computed.Add(state, (i+1)*dt);
    }
    
    double radius = parameters("radius");
    double v0 = parameters("v0");
    double vMax = parameters("vMax");
    double seed = parameters("seed");
    
    DTRandom random(seed);
    
//    //Initializes the particles positions and velocities
//    int numPoints = initialParticles.NumberOfPoints();
//    for (int i = 0; i < numPoints; i++) {
//        DTPoint2D point = initialParticles(i);
//        particlePositions[i][0] = point.x;
//        particlePositions[i][1] = point.y;
//        
//        double angle = r.UniformHalf()*2*M_PI;
//        particleVelocities[i][0] = cos(angle)*v0; //initial vx
//        particleVelocities[i][1] = sin(angle)*v0; //inital vy
//    }
    
    DTProgress progress;
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
}
