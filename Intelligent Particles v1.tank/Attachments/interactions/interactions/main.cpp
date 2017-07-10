// #include "DTSource.h"
#include "DTSaveError.h"

#include "DTArguments.h"
#include "DTDataFile.h"
#include "DTDictionary.h"
#include "DTPointValueCollection2D.h"
#include "DTProgress.h"
#include "DTSeriesPointValueCollection2D.h"

//////////////////////////////////////////////////////////////////////////////
//    Main routine
//////////////////////////////////////////////////////////////////////////////

void Computation(DTSeriesPointValueCollection2D &computed);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTDataFile outputFile("Output.dtbin",DTFile::NewReadWrite);
    // Input variables.
    
    // Output series.
    DTSeriesPointValueCollection2D computed(outputFile,"Var");
    
    // The computation.
    clock_t t_before = clock();
    Computation(computed);
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

void Computation(DTSeriesPointValueCollection2D &computed)
{
    computed.Add(<#DTPointValueCollection2D v#>, <#double time#>)
    
    
    
    
    
    
    
    // Insert your code here.
    
    DTProgress progress;
    
    // Inside the loop, do
    //     progress.UpdatePercentage(fraction);
    //     computed.Add(returnStructure,time); // Call with time>=0 and strictly increasing.
    
}
