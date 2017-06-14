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

//simulation parameters
int N = 64; // number of particles
double rho = .5; // density (number per unit volume)
double T = 1.0; // temperature

// function declarations

void initialize();        // allocates memory, calls following 2 functions
void initPositions();     // places particles on an fcc lattice
void initVelocities();    // initial Maxwell-Boltzmann velocity distribution
void rescaleVelocities(); // adjust the instanteous temperature to T
double gasdev();          // Gaussian distributed random numbers

double **r; // positions
double **v; // velocities
double **a; // accelerations

void initialize() {
    r = new double* [N];
    v = new double* [N];
    a = new double* [N];
    for (int i = 0; i < N; i++) {
        r[i] = new double [3];
        v[i] = new double [3];
        a[i] = new double [3];
    }
    initPositions();
    initVelocities();
}

double L;  // linear size of cubical volume

void initPositions() {
    
    // compute side of cube from number of particles and number density
    L = pow(N / rho, 1.0/3);
    
    // find M large enough to fit N atoms on an fcc lattice
    int M = 1;
    while (4 * M * M * M < N)
        ++M;
    double a = L / M;           // lattice constant of conventional cell
    
    // 4 atomic positions in fcc unit cell
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};   //
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};   //
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};   //
    
    int n = 0;                  // atoms placed so far
    for (int x = 0; x < M; x++) //
        for (int y = 0; y < M; y++) //
            for (int z = 0; z < M; z++) //
                for (int k = 0; k < 4; k++) 
                    if (n < N) {
                        r[n][0] = (x + xCell[k]) * a;
                        r[n][1] = (y + yCell[k]) * a;
                        r[n][2] = (z + zCell[k]) * a;
                        ++n;
                    }
}

double gasdev () {
    static bool available = false;  //set to false initially
    static double gset;             //
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        return v2*fac;
    } else {
        available = false;
        return gset;
    }
}

void initVelocities() {
    
    // Gaussian with unit variance
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev();
    
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0};
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i];
    for (int i = 0; i < 3; i++)
        vCM[i] /= N;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] -= vCM[i];
    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities();
}

void rescaleVelocities() {
    double vSqdSum = 0;
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vSqdSum += v[n][i] * v[n][i];
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum );
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda;
}

void computeAccelerations() {
    for (int i = 0; i < N; i++)
        for (int k = 0; k < 3; k++)
            a[i][k] = 0;
    for (int i = 0; i < N-1; i++)
        for (int j = i+1; j < N; j++) {
            double rij[3];
            double rSqd = 0;
            for (int k = 0; k < 3; k++) {
                rij[k] = r[i][k] - r[j][k];
                // closest image convention
                if (abs(rij[k]) > 0.5 * L) {
                    if (rij[k] > 0)
                        rij[k] -= L;
                    else
                        rij[k] += L;
                }
                rSqd += rij[k] * rij[k];
            }
            double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
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
            // use periodic boundary conditions
            if (r[i][k] < 0)
                r[i][k] += L;
            if (r[i][k] >= L)
                r[i][k] -= L;
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
        if (i % 200 == 0)
            rescaleVelocities();
        cout << instantaneousTemperature() << '\n';
        for(int j = 0; j< N; j++) {
            position(0,j) = r[j][0];
            position(1,j) = r[j][1];
            position(2,j) = r[j][2];
            velocity(0,j) = v[j][0];
            velocity(1,j) = v[j][1];
            velocity(2,j) = v[j][2];
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
