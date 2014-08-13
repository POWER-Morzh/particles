#include "tarch/logging/Log.h"
#include "tarch/tests/TestCaseRegistry.h"
#include "tarch/parallel/Node.h"

#include "peano/peano.h"

#include "particles/ExperimentSettings.h"

#include "particles/pit/runners/Runner.h"
#include "particles/pidt/runners/Runner.h"

#include <cstdlib>


#ifdef SharedTBB
#include "tarch/multicore/tbb/Core.h"
#endif


tarch::logging::Log _log("");


void usage() {
  std::cerr << "usage:"  << std::endl
            << "executable type id particle-density ppc dt dam-width [cores]" << std::endl
            << "   type             type of particle code, i.e. either pit, pidt, or rapidt (reduction-avoiding pidt)" << std::endl
            << "   id               id of experiment" << std::endl
            << "   particle-density particle density (double)" << std::endl
            << "   ppc              number of particles per cell (int)" << std::endl
            << "   dt               time step size dt=(0,1(" << std::endl
            << "   dam-width        dam width (<=1.0, double)" << std::endl
            << "   cores            number of cores to use (TBB only)" << std::endl
            << std::endl
            << "id has to be from "
            << particles::ExperimentSettings::getListOfAvailableExperiments() << std::endl;
}



int main(int argc, char** argv) {
  peano::fillLookupTables();

  int parallelSetup = peano::initParallelEnvironment(&argc,&argv);
  if ( parallelSetup!=0 ) {
    #ifdef Parallel
    // Please do not use the logging if MPI doesn't work properly.
    std::cerr << "mpi initialisation wasn't successful. Application shut down" << std::endl;
    #else
    _log.error("main()", "mpi initialisation wasn't successful. Application shut down");
    #endif
    return parallelSetup;
  }

  int sharedMemorySetup = peano::initSharedMemoryEnvironment();
  if (sharedMemorySetup!=0) {
    logError("main()", "shared memory initialisation wasn't successful. Application shut down");
    return sharedMemorySetup;
  }

  int programExitCode = 0;

  int requiredParameters = 7;
  #ifdef SharedTBB
  requiredParameters++;
  #endif

  if (argc!=requiredParameters) {
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      usage();
    }
    exit(-1);
  }

  // executable type id particle-density ppc dt dam-width [cores]
  const std::string type             = argv[1];
  const std::string experimentID     = argv[2];
  const double      particleDensity  = atof(argv[3]);
  const int         ppc              = atoi(argv[4]);
  const double      dt               = atof(argv[5]);
  const double      breakingDamWidth = atof(argv[6]);

  #ifdef SharedTBB
  const int numberOfCores = atoi(argv[7]);
  tarch::multicore::tbb::Core::getInstance().configure(numberOfCores);
  #endif


  particles::ExperimentSettings::configure(experimentID,particleDensity,ppc,dt,breakingDamWidth);
  if (particles::ExperimentSettings::chosenConfiguration==0) {
    logError( "main()", "Not a valid experiment variant" );
    usage();
    programExitCode = 1;
  }

  if (type!="pit" && type!="pidt" && type!="rapidt" ) {
    logError( "main()", "Not a valid particle code type. Use pit, pidt, or rapidt" );
    programExitCode = 2;
  }

  if (programExitCode==0 && type=="pit") {
    particles::pit::runners::Runner runner;
    switch (particles::ExperimentSettings::chosenConfiguration->runnerVariant) {
      case particles::ExperimentSettings::UnitTests:
        tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().run();
        programExitCode  = tarch::tests::TestCaseRegistry::getInstance().getTestCaseCollection().getNumberOfErrors();
        tarch::tests::TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().run();
        programExitCode += tarch::tests::TestCaseRegistry::getInstance().getIntegrationTestCaseCollection().getNumberOfErrors();
        break;
      case particles::ExperimentSettings::Simulation:
        programExitCode = runner.run();
        break;
      case particles::ExperimentSettings::MemoryBenchmark:
        programExitCode = runner.runMemoryBenchmark();
        break;
    }
  }


  if (programExitCode==0 && (type=="pidt"||type=="rapidt") ) {
    particles::pidt::runners::Runner runner;
    switch (particles::ExperimentSettings::chosenConfiguration->runnerVariant) {
      case particles::ExperimentSettings::UnitTests:
        logError("main()", "unit tests are not realised for pidt. Use pit instead");
        break;
      case particles::ExperimentSettings::Simulation:
        programExitCode = runner.run();
        break;
      case particles::ExperimentSettings::MemoryBenchmark:
        logError("main()", "memory benchmark is not realised for pidt. Use pit instead");
        break;
    }
  }

  if (programExitCode==0) {
    #ifdef Parallel
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      logInfo( "main()", "Peano terminates successfully" );
    }
    #else
    logInfo( "main()", "Peano terminates successfully" );
    #endif
  }
  else {
    logInfo( "main()", "quit with error code " << programExitCode );
  }

  std::cout << "Simulation complete" << std::endl;
  peano::shutdownParallelEnvironment();
  peano::shutdownSharedMemoryEnvironment();

  return programExitCode;
}
