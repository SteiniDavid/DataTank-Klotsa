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
#include "stdlib.h"
#include "DTVector2D.h"


int createHomeBeacon(const DTPointCollection2D &initial, const DTDictionary &parameters, int seed);


class Agent {
public:
    enum AgentType {FOLLOWER, BEACON, LOOKING, FOUNDIT};
    
    Agent() : type(FOLLOWER) {}
    Agent(DTPoint2D l, DTPoint2D meal, double noise, double radius) : x(l.x), y(l.y), type(FOLLOWER), u(0), v(0), food(meal), noise(noise) , radius(radius){}
    
    AgentType Type(void) const {
        return type;
    }
    
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
    
    DTPoint2D location(void) const {
        return DTPoint2D(x,y);
    }
    
    //Move command is what is called in computation
    //Its the most basic version and everything goes out of it.
    void Move(double sqrtDT, DTRegion2D box, DTRandom randNumber, const Agent &beacon) {
        //Figure out which type of advance method to use and hand in the right things.
        double xy[2];
        randNumber.Normal(xy,2);
        if (Type() == FOLLOWER) { //when looking for food its a random walk
            //Conversion rate to LOOKING type
            //Eventually it will be based on how close you are to the beacon
            //Closer meaning lower probability of conversion
            int chance = rand() % 10000 + 1;
            if (chance < 2) {
                ChangeToLOOKING();
            }
            //If not converted it moves
            advanceFOLLOWER(box, beacon);
        }
        if (Type() == LOOKING) {
           // advanceLOOKING(sqrtDT*noise*xy[0], sqrtDT*noise*xy[1], box);
        }
        if (Type() == FOUNDIT) {
            
        }
        if (Type() == BEACON) {
            //Don't move
        }
        
    }
    
    void advanceLOOKING(double dx,double dy,const DTRegion2D &box) {
        //This handles the random walk
        x += dx;
        y += dy;
        
        //add that they can't go through each other "hard shelled"
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
        //The agent object was instantiated with the location of the food
        //uses it to figure out if it found it within its given radius of search
        double distanceToFood = Norm(location()-food);
        if (distanceToFood <= radius*radius) {
            ChangeToFOUNDIT(); //the agent has now found the food
        }
    }
    
    void advanceFOLLOWER(const DTRegion2D &box, Agent beacon) {
//        //Have them headed right at the beacon, but have noise so it just swarms around it.
//        
//        //compare this distance with the manipulated dx and dy to see which way will get you closer
//        double distBeaconFollower = Norm(location()-beacon.location());
//        
//        double xFollower = beacon.x+x;
//        double yFollower = beacon.y+y;
//        //Create these as points
//        DTPoint2D xChanged(xFollower,y);
//        DTPoint2D yChanged(x,yFollower);
//        DTPoint2D xyChanged(xFollower,yFollower);
//        //the distances for the modifications
//        double distxChanged = Norm(xChanged-beacon.location());
//        double distyChanged = Norm(yChanged-beacon.location());
//        double distxyChanged = Norm(xyChanged-beacon.location());
//        
//        if (distxyChanged < distBeaconFollower) {
//            if (distxChanged) {
//                
//            }
//            x += dx;
//            y += dy;
//        } else {
//            x -= dx;
//            y -= dy;
//        }
        
        //second try, this time the closer to the beacon less variance and smaller step
//        double dx = beacon.x-x;
//        double dy = beacon.y-y;
//        double distBeaconFollower = sqrt(Norm(location()-beacon.location()))*1000;
//        if(dx > 0) {
//            x -= dx/distBeaconFollower;
//        } else {
//            x += dx/distBeaconFollower;
//        }
//        if(dy > 0) {
//            y -= dy/distBeaconFollower;
//        } else {
//            y += dy/distBeaconFollower;
//        }
        //third try
        
        
        
    }
    
    void advanceFOUNDIT(double dx,double dy,const DTRegion2D &box, Agent beaconAgent) {
        
    }
    
private:
    double x,y; // Coordinates
    double u,v; // Velocity
    double radius;
    double noise;
    DTPoint2D food;
    AgentType type;
};

void Computation(const DTPointCollection2D &initial,const DTPoint2D &food,
                 double dt,double endTime,int seed,
                 const DTDictionary &parameters,const DTRegion2D &region,
                 DTSeriesGroup<DT_RetGroup> &computed)
{
    DTMutablePointCollection2D ants = initial.Copy();
    int numAnts = ants.NumberOfPoints();
    
    DTMutableDoubleArray antArray = ants.Data();
    
    DT_RetGroup state;
    DTRandom randNumber(seed);
    double sqrtDT = sqrt(dt);
    double xy[2];
    
    double noise = parameters("noise");
    double radius = parameters("radius");
    
    DTMutableList<Agent> antList(numAnts);
    for (int i=0;i<numAnts;i++) {
        //every agent knows where it is and where the food is, although it cant directly access it
        // FIXME: Make it so they are not passed the food, I think that makes sense
        antList(i) = Agent(DTPoint2D(ants(i)), food, noise, radius);
    }
    int beaconID = createHomeBeacon(initial, parameters, seed);
    antList(beaconID).ChangeToBeacon(); //converts the one that is picked as beacon to beacon
    
    double t = 0;
    while (t<endTime) {
        for (int i = 0; i<numAnts; i++) {
            //Random Walk
            randNumber.Normal(xy,2);
            antList(i).Move(sqrtDT,region,randNumber, antList(beaconID));
        }
        
        t += dt;
        
        // Save
        // Start by converting your antList into structures that can be saved to the file
        DTMutableDoubleArray xyPoints(2,numAnts);
        DTPoint2D P;
        for (int i=0;i<numAnts;i++) {
            P = antList(i).location();
            xyPoints(0,i) = P.x;
            xyPoints(1,i) = P.y;
        }
        state.position = DTPointCollection2D(xyPoints);
        state.foodBeacon;
        state.homeBeacon;
        computed.Add(state,t);
    }
    
    
    DTProgress progress;
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
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





