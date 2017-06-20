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
    r = new double* [N]; //the positions
    v = new double* [N]; //the velocity vectors
    a = new double* [N]; //the acceleration vectors
    for (int i = 0; i < N; i++) { //goes over for every particle
        r[i] = new double [3]; //creates an array for the three dimensions
        v[i] = new double [3]; //
        a[i] = new double [3];
    }
    initPositions();  //initializes the positions now that its done creating the arrays
    initVelocities(); //initializes the random inital velocities
}

double L;  // linear size of cubical volume

void initPositions() {
    
    // compute side of cube from number of particles and number density
    L = pow(N / rho, 1.0/3); //the side length is based on the number of particles as well as the density
    
    // find M large enough to fit N atoms on an fcc lattice
    int M = 1; //
    while (4 * M * M * M < N) //so that all of the particles fit within the lattice
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
                for (int k = 0; k < 4; k++) //
                    if (n < N) { //if the number placed so far is less than the number of particles continue
                        r[n][0] = (x + xCell[k]) * a; //places particle in x part
                        r[n][1] = (y + yCell[k]) * a; //for y part
                        r[n][2] = (z + zCell[k]) * a; //for z part
                        ++n; //the number of particles placed has increased by one
                    }
}

DTMutablePointCollection3D InitializePosition(int NumPoints,double rhoVal) { //
    
    // compute side of cube from number of particles and number density
   double cubeVol = pow(NumPoints / rhoVal, 1.0/3); //the side length is based on the number of particles as well as the density
    
    // find M large enough to fit N atoms on an fcc lattice
    int M = 1; //
    while (4 * M * M * M < N) //so that all of the particles fit within the lattice
        ++M;
    double a = cubeVol / M;           // lattice constant of conventional cell
    
    // 4 atomic positions in fcc unit cell
    double xCell[4] = {0.25, 0.75, 0.75, 0.25};   //
    double yCell[4] = {0.25, 0.75, 0.25, 0.75};   //
    double zCell[4] = {0.25, 0.25, 0.75, 0.75};   //
    
    DTMutableDoubleArray points(3,NumPoints);
    
    int n = 0;                  // atoms placed so far
    for (int x = 0; x < M; x++) {//
        for (int y = 0; y < M; y++) {//
            for (int z = 0; z < M; z++) {//
                for (int k = 0; k < 4; k++) {//
                    if (n < NumPoints) { //if the number placed so far is less than the number of particles continue
                        points(0,n) = (x + xCell[k]) * a; //places particle in x part
                        points(1,n) = (y + yCell[k]) * a; //for y part
                        points(2,n) = (z + zCell[k]) * a; //for z part
                        ++n; //the number of particles placed has increased by one
                    }
                }
            }
        }
    }
    
    return DTMutablePointCollection3D(points);
}

double gasdev () {
    static bool available = false;  //set to false initially
    static double gset;             //initialize the gas distribution number
    double fac, rsq, v1, v2; //
    if (!available) { //if available is false then do this
        do { //excecutes until while is false
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0; //v1 is 2 times a scaled random num
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0; //save for v2
            rsq = v1 * v1 + v2 * v2; //finds the square of the magnitude
        } while (rsq >= 1.0 || rsq == 0.0); //keep going until rsq is not greater than 1 or its not zero
        fac = sqrt(-2.0 * log(rsq) / rsq); //
        gset = v1 * fac; //the gset is v1 times fac
        available = true; //available set to true
        return v2*fac; //return fac times v2
    } else {
        available = false; //set available to false
        return gset; //returns gset instead
    }
}

void initVelocities() {
    
    // Gaussian with unit variance
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] = gasdev(); //calls gasdev funciton for each particle for each of the three dim
    
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0}; //the center of mass velocity is set as zero
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i]; //becomes the total velocity for a given dim
    for (int i = 0; i < 3; i++) //goes through 3 dim
        vCM[i] /= N; //devides each dims total velocity by the num of particles
    for (int n = 0; n < N; n++) //for all of the particles in the system
        for (int i = 0; i < 3; i++) //for the three dim
            v[n][i] -= vCM[i]; //subtracts the velocty average for a given dim form each particles velocity in that dim
    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities(); //rescales all of the velocities to get back to the right speed
}

DTMutablePointCollection3D initializeVelocities(int numPoints) {
    DTMutableDoubleArray points(3,numPoints);
    
    // Gaussian with unit variance
    for (int n = 0; n < numPoints; n++) {
        for (int i = 0; i < 3; i++) {
            points(i,n) = gasdev(); //calls gasdev function for each particle for each of the three dim
        }
    }
    // Adjust velocities so center-of-mass velocity is zero
    double vCM[3] = {0, 0, 0}; //the center of mass velocity is set as zero
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            vCM[i] += v[n][i]; //becomes the total velocity for a given dim
    for (int i = 0; i < 3; i++) //goes through 3 dim
        vCM[i] /= N; //devides each dims total velocity by the num of particles
    for (int n = 0; n < N; n++) //for all of the particles in the system
        for (int i = 0; i < 3; i++) //for the three dim
            points(i,n) -= vCM[i]; //subtracts the velocty average for a given dim form each particles velocity in that dim
    // Rescale velocities to get the desired instantaneous temperature
    rescaleVelocities(); //rescales all of the velocities to get back to the right speed
    
    return DTMutablePointCollection3D(points);
}

void rescaleVelocities() {
    double vSqdSum = 0; //velocity square root initialized at 0
    for (int n = 0; n < N; n++) //for all of the particles
        for (int i = 0; i < 3; i++) //for the three dim
            vSqdSum += v[n][i] * v[n][i]; //add the magnitude of all of the particles
    double lambda = sqrt( 3 * (N-1) * T / vSqdSum ); //value to scale all of the vectors by
    for (int n = 0; n < N; n++)
        for (int i = 0; i < 3; i++)
            v[n][i] *= lambda; //scales all of the vectors dims by lambda
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

double instantaneousTemperature() {
    double sum = 0; //initially sum is zero
    for (int i = 0; i < N; i++) //all particles
        for (int k = 0; k < 3; k++) //all dim
            sum += v[i][k] * v[i][k]; //adds up the magnitude of all of the dims of all of the particles
    return sum / (3 * (N - 1)); //returns the temp of the system
}




void Computation(const DTDictionary &parameters,
                 const DTPointCollection3D &initialPoints,
                 const DTRegion3D &box,DTSeriesGroup<DT_RetGroup> &computed)
{
    double numPoints = initialPoints.NumberOfPoints();
    
    //My methods version
    DTMutablePointCollection3D positions = InitializePosition(numPoints, rho);
    DTMutablePointCollection3D velocities = initializeVelocities(numPoints);
    
    
    
    
    
    
    //////////////////////////////
    
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
    int stride = parameters("stride");
    for (int i = 0; i < maxTime/dt; i++) {
        velocityVerlet(dt);
        if (i % 200 == 0)
            rescaleVelocities();
        cout << instantaneousTemperature() << '\n';
        
        if ((i+1)%stride==0) {
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
