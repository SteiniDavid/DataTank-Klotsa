// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPoint2D.h"
#include "DTPointCollection2D.h"
#include "DTRegion2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"

#include <math.h>
#include "DTRandom.h"

DTPointCollection2D Computation(const DTRegion2D &region,const DTDictionary &parameters,
                                int seed,double dt,double endTime,const DTPointCollection2D &initial,
                                const DTPoint2D &food);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTRegion2D region;
    Read(inputFile,"region",region);
    DTDictionary parameters;
    Read(inputFile,"parameters",parameters);
    int seed = int(inputFile.ReadNumber("seed"));
    double dt = inputFile.ReadNumber("dt");
    double endTime = inputFile.ReadNumber("endTime");
    DTPointCollection2D initial;
    Read(inputFile,"initial",initial);
    DTPoint2D food;
    Read(inputFile,"food",food);
    
    // The computation.
    DTPointCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(region,parameters,seed,dt,endTime,initial,food);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    
    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("PointCollection2D","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    outputFile.SaveIndex();
    
    return 0;
}

int createHomeBeacon(const DTPointCollection2D &initial, const DTDictionary &parameters, int seed);

class Ant {
public:
    enum AgentType {FOLLOWER, BEACON, LOOKING, FOUNDIT};
    
    //Constructors
    Ant() : type(FOLLOWER) {}
    Ant(DTPoint2D l, DTDictionary dict, DTPoint2D food, DTRegion2D region, DTRandom r) : x(l.x), y(l.y), dict(dict), food(food), box(region), type(FOLLOWER) {
        //Setting random vx and vy between -.1 and .1
        double angle = r.UniformHalf()*2*M_PI;
        double v0 = dict("v0");
        vx = cos(angle)*v0; //initial vx
        vy = sin(angle)*v0; //inital vy
    }
    
    //Getters
    AgentType Type(void) const {
        return type;
    }
    DTPoint2D location(void) const {
        return DTPoint2D(x,y);
    }
    
    //Setters
    void ChangeToBeacon(void) {
        type = BEACON;
    }
    void ChangeToFOUNDIT(void) {
        type = FOUNDIT;
    }
    void ChangeToLOOKING(void) {
        type = LOOKING;
    }
    void ChangeToFOLLOWER(void) {
        type = FOLLOWER;
    }
    
    //Movement functions
    void Move(DTRandom r, Ant beacon, DTMutableList<Ant> ants) {
        if (Type() == FOLLOWER) {
            advanceFOLLOWER(r, beacon, ants);
        }
        if (Type() == LOOKING) {
            advanceLOOKING(r);
        }
        if (Type() == FOUNDIT) {
            
        }
        if (Type() == BEACON) {
            //Don't move
        }
        
    }
    
    void advanceFOLLOWER(DTRandom r, Ant beacon, DTMutableList<Ant> ants) {
        
        collisionDetection(ants);
        
        bounding();
        
    }
    void advanceLOOKING(DTRandom r){
        double xy[2];
        r.Normal(xy,2); //generates a random x and y coordinate that will later be scaled
        double sqrtDT = sqrt(dt);
        double dx = sqrtDT*noise*xy[0];
        double dy = sqrtDT*noise*xy[1];
        //This handles the random walk
        x += dx;
        y += dy;
        bounding();
        
        double distanceToFood = Norm(location()-food);
        if (distanceToFood <= radius*radius) {
            ChangeToFOUNDIT(); //the agent has now found the food
        }
    }
    
    //check if there was a collision between particles and handles particle interaction
    void collisionDetection(DTMutableList<Ant> ants) {
        double rad = radius*radius;
        for (int i = 0; i < ants.Length(); i++) {
            double xDist = ants(i).x-x;
            double yDist = ants(i).y-y;
            double distance = xDist*xDist+yDist*yDist;
            if (distance < rad && xDist != 0 && yDist != 0) {
                //overlapping then lenard jones potential equation is used, in form of force
                
                double epsilon = 1;
                double sigma = radius*2;
                double force = ((24*epsilon)/(sigma))*(2*pow(sigma/distance, 13.0)-pow(sigma/distance, 7.0));
                double theta = atan2(vy, vx);
                double fx = force*cos(theta);
                double fy = force*sin(theta);
                
                //the exact same force negated impacts the other particle
                
                
                //double potential = 4*epsilon*(pow(radius*2/distance, 12.0)-pow(radius*2/distance, 6.0));
                //this force is directed in reverse of the xDist and yDist and scaled by the force
                //x -= xDist*force; //scaled by force for now since i dont know how to handle it better
                //y -= yDist*force;
            }
        }
    }
    
    
    void bounding() { //Bounds the x and y within the region
        if (x > box.xmax) {
            x = 2*box.xmax-x;
        }
        if (y > box.ymax) {
            y = 2*box.ymax-y;
        }
        if (x < box.xmin) {
            x = 2*box.xmin-x;
        }
        if (y < box.ymin) {
            y = 2*box.ymin-y;
        }
    }
    
    
private:
    double x,y; // Coordinates
    double vx, vy; //velocity's
    DTDictionary dict;
    double sight = dict("sight");
    double radius = dict("radius");
    DTPoint2D food;
    DTRegion2D box;
    double dt = dict("dt");
    double noise = dict("noise");
    double vMax = dict("vMax");
    
    AgentType type;
};


DTPointCollection2D Computation(const DTRegion2D &region,const DTDictionary &parameters,
                                int seed,double dt,double endTime,const DTPointCollection2D &initial, const DTPoint2D &food)
{
    DTMutablePointCollection2D ants = initial.Copy();
    int numAnts = ants.NumberOfPoints();
    DTMutableDoubleArray antArray = ants.Data();
    
    DTRandom r(seed);
    
    DTMutableList<Ant> antList(numAnts);
    for (int i=0;i<numAnts;i++) {
        antList(i) = Ant(DTPoint2D(ants(i)), parameters, food, region, r);
    }
    
    int beaconID = createHomeBeacon(initial, parameters, seed);
    antList(beaconID).ChangeToBeacon(); //converts the one that is picked as beacon to beacon
    
    
    double t = 0;
    while (t<endTime) {
        for (int i = 0; i<numAnts; i++) {
            antList(i).Move(r, antList(beaconID), antList);
        }
    
        t += dt;
    }
    
    
    
    
    
    
    
    
    
    
    DTPointCollection2D toReturn;

    // region.isSet - false means invalid.
    // region.xmin, xmax, ymin, ymax

    // type c = parameters("name") where type is double, string, DTDoubleArray or DTDictionary.

    // initial.Data() - DTDoubleArray - a 2xN array.
    // initial.PointNumbers() - DTIntArray - optional, list with N entries.

    return toReturn;
}

int createHomeBeacon(const DTPointCollection2D &initial, const DTDictionary &parameters, int seed){
    DTMutableIntArray swarmID(initial.NumberOfPoints());
    DTMutableIntArray swarmBiggestID(initial.NumberOfPoints());
    double radius = parameters("radius");
    DTDoubleArray antArray = initial.Data();
    int numPoints = initial.NumberOfPoints();
    DTRandom randNumber(seed);
    
    //Everyone generates their ids, assuming that none will be the same and the biggest
    for (int i = 0; i < numPoints; i++) {
        swarmID(i) = randNumber.UInteger();
    }
    
    //Attempt to do it with radius of interaction
    DTMutableDoubleArray dist(1000,1000);
    for (int i = 0; i < numPoints; i++) {
        DTPoint2D ithPoint = initial(i);
        for (int j = 0; j < numPoints; j++) {
            DTPoint2D jthPoint = initial(j);
            dist(i,j) = Norm(ithPoint-jthPoint);
        }
    }
    //Now that I have the distances between everyone I can make people only be able to talk to
    //neighbors within their specified radius
    swarmBiggestID = swarmID.Copy();
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < numPoints; j++) {
            if (dist(i,j) <= radius) {
                int ithID = swarmBiggestID(i);
                int jthID = swarmBiggestID(j);
                if(ithID > jthID) {
                    swarmBiggestID(j) = ithID;
                } else {
                    swarmBiggestID(i) = jthID;
                }
            }
        }
    }
    int beaconParticleID = -1;
    //Now we have to find which point has the biggest ID
    for (int i = 0; i < numPoints-1; i++) {
        int swarmsID = swarmID(i);
        int swarmsBiggestID = swarmBiggestID(2);
        if (swarmsBiggestID == swarmsID) {
            beaconParticleID = i;
        }
    }
    return beaconParticleID;
}

