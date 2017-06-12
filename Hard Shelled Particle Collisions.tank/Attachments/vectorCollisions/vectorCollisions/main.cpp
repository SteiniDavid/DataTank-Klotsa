// Include DTSource.h if you want to include all the headers.

#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointCollection2D.h"
#include "DTRegion2D.h"
#include "DTVectorCollection2D.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"

#include <math.h>

DTVectorCollection2D Computation(const DTPointCollection2D &initialParticles,const DTRegion2D &_2D_Region,
                                 const DTDictionary &parameters);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);

    DTDataFile inputFile("Input.dtbin",DTFile::ReadOnly);
    // Read in the input variables.
    DTPointCollection2D initialParticles;
    Read(inputFile,"initialParticles",initialParticles);
    DTRegion2D _2D_Region;
    Read(inputFile,"2D Region",_2D_Region);
    DTDictionary parameters;
    Read(inputFile,"parameters",parameters);

    // The computation.
    DTVectorCollection2D computed;
    clock_t t_before = clock();
    computed = Computation(initialParticles,_2D_Region,parameters);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);

    // Write the output.
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);

    // Output from computation
    Write(outputFile,"Var",computed);
    outputFile.Save("VectorCollection2D","Seq_Var");

    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");

    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");

    outputFile.SaveIndex();

    return 0;
}

DTVectorCollection2D Computation(const DTPointCollection2D &initialParticles,const DTRegion2D &_2D_Region,
                                 const DTDictionary &parameters)
{
    DTVectorCollection2D toReturn;

    // initialParticles.Data() - DTDoubleArray - a 2xN array.
    // initialParticles.PointNumbers() - DTIntArray - optional, list with N entries.

    // _2D_Region.isSet - false means invalid.
    // _2D_Region.xmin, xmax, ymin, ymax

    // type c = parameters("name") where type is double, string, DTDoubleArray or DTDictionary.

    return toReturn;
}
