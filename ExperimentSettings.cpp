#include "particles/ExperimentSettings.h"

#include "tarch/Assertions.h"

#include "tarch/logging/CommandLineLogger.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"

#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"

#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#ifdef SharedMemoryParallelisation
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"
#endif

#ifdef Parallel
#include "mpibalancing/FairNodePoolStrategy.h"
#include "mpibalancing/StaticBalancing.h"
#endif

#include <string>
#include <sstream>

#include "particles/pidt/State.h"


particles::ExperimentSettings* particles::ExperimentSettings::chosenConfiguration = 0;



void particles::ExperimentSettings::configure(
  const std::string& id,
  const double       particleDensity,
  const int          ppc,
  const double       dt,
  const double       breakingDamWidth
) {
  tarch::logging::Log _log( "particles::ExperimentSettings" );
  assertion(particles::ExperimentSettings::chosenConfiguration==0);

  if (id==particles::ExperimentSettingsRunTests::ID) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsRunTests();
  }
  else if (id==particles::ExperimentSettingsProduceVideo::ID) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsProduceVideo(dt);
  }
  else if (id.find(particles::ExperimentSettingsHomogeneousParticles::ID)!=std::string::npos) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsHomogeneousParticles(
      id.find("split-tree")!=std::string::npos,
      id.find("pipeline")!=std::string::npos,
      id.find("skip-reduce")==std::string::npos
    );
  }
  else if (id==particles::ExperimentSettingsMemoryBenchmark::ID) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsMemoryBenchmark();
  }
  else if (id==particles::ExperimentSettingsPureGrid::ID) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsPureGrid();
  }
  else if (id==particles::ExperimentSettingsStationaryParticles::ID) {
    particles::ExperimentSettings::chosenConfiguration = new particles::ExperimentSettingsStationaryParticles();
  }

  if (particles::ExperimentSettings::chosenConfiguration!=0) {
    particles::ExperimentSettings::chosenConfiguration->particleDensity  = particleDensity;
    particles::ExperimentSettings::chosenConfiguration->ppc              = ppc;
    particles::ExperimentSettings::chosenConfiguration->dt               = dt;
    particles::ExperimentSettings::chosenConfiguration->breakingDamWidth = breakingDamWidth;

    if (particleDensity<=0) {
      logWarning( "configure(...)", "no particles in system" );
    }
  }
}


std::string particles::ExperimentSettings::getListOfAvailableExperiments() {
  std::ostringstream msg;

  msg <<  "\"" << particles::ExperimentSettingsRunTests::ID << "\""
      << ",\"" << particles::ExperimentSettingsProduceVideo::ID << "\""
      << ",\"" << particles::ExperimentSettingsHomogeneousParticles::ID << ",[a],[b],[c]\""
      << ",\"" << particles::ExperimentSettingsMemoryBenchmark::ID << "\""
      << ",\"" << particles::ExperimentSettingsPureGrid::ID << "\""
      << ",\"" << particles::ExperimentSettingsStationaryParticles::ID << "\"";

  msg << "     [a=split-tree|no],[b=pipeline|no],[c=skip-reduce|no]";

  return msg.str();
}


particles::ExperimentSettingsRunTests::ExperimentSettingsRunTests() {
  runnerVariant = UnitTests;

  numberOfTimeSteps      = 1;
  maximalInitialVelocity = 1.0;

  timeSteppingAdapter = ExplicitEuler;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", false ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
}


particles::ExperimentSettingsProduceVideo::ExperimentSettingsProduceVideo(double dt) {
  runnerVariant = Simulation;

  numberOfTimeSteps      = 5.0 / dt;
  maximalInitialVelocity = 1.0;
//  maximalInitialVelocity = 0.1;

  timeSteppingAdapter = ExplicitEulerWithPlot;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();

  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );

  #if defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( tarch::logging::CommandLineLogger::FilterListEntry( "info", false ) );
  std::ostringstream traceFileName;
  traceFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.out";
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry( "debug", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
    "particles::pit",
    false
  ));
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, traceFileName.str() );
  #elif defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", false ) );
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
  #else
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", true ) );
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, false, "" );
  #endif

  peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
    new peano::datatraversal::autotuning::OracleForOnePhaseDummy(true,false,true,false,false)
  );


  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().setStrategy(
      new tarch::parallel::FCFSNodePoolStrategy()
    );
  }
  tarch::parallel::NodePool::getInstance().restart();

  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
    new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(true)
  );

  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(64);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(64);

  tarch::parallel::Node::getInstance().setDeadlockTimeOut(120);
  tarch::parallel::Node::getInstance().setTimeOutWarning(60);
}



particles::ExperimentSettingsHomogeneousParticles::ExperimentSettingsHomogeneousParticles(bool splitTree, bool pipeline, bool alwaysReduce) {
  runnerVariant = Simulation;

  numberOfTimeSteps      = 50;
//  maximalInitialVelocity = 1.0;
  maximalInitialVelocity = 0.1;

  timeSteppingAdapter = ExplicitEuler;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();

  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "", true ) );


  #if defined(Debug)
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "", false ) );
  std::ostringstream traceFileName;
  traceFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.out";
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry( "debug", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
    "particles::pit::mappings",
    false
  ));
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry( "debug", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
    "particles::pidt::mappings",
    false
  ));
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, traceFileName.str() );
  #elif defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "", false ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry(
    "info", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
    "peano::parallel", true
  ));
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
  #else
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "", false ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry(
      "info", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
      "peano::parallel", true
    ));
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry(
      "info", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
      "peano::heap", true
   ));
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry(
      "info", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
      "mpibalancing", true
    ));
  /// vorletztes false
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
  #endif

  // All these values are heuristics found by trial-and-error
  peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
    new peano::datatraversal::autotuning::OracleForOnePhaseDummy(
      true,
      false,
      splitTree ? 1 : 0,
      pipeline,              // descend
      pipeline,              // ascend
      THREE_POWER_D,         // smallestGrainSizeForAscendDescend,
      1,                     // grainSizeForAscendDescend
      tarch::la::aPowI(DIMENSIONS,9/2), //smallestGrainSizeForEnterLeaveCell (default)
      2,                     // grainSizeForEnterLeaveCell (default)
      FOUR_POWER_D,          // smallestGrainSizeForTouchFirstLast
      4                      // grainSizeForTouchFirstLast
    )
  );

  pidt::State::AlwaysReduceAlongMasterWorkers = alwaysReduce;

  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().setStrategy(
      new mpibalancing::FairNodePoolStrategy(6)
    );
  }
  #else
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().setStrategy(
      new tarch::parallel::FCFSNodePoolStrategy()
    );
  }
  #endif

  tarch::parallel::NodePool::getInstance().restart();

  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
    new mpibalancing::StaticBalancing(true)
  );
  #else
  peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
    new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(true)
  );
  #endif

  // have to be the same for all ranks
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(16);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(16);

  #if defined(Debug) || defined(Asserts)
  tarch::parallel::Node::getInstance().setDeadlockTimeOut(120*4);
  tarch::parallel::Node::getInstance().setTimeOutWarning(60*4);
  #else
  tarch::parallel::Node::getInstance().setDeadlockTimeOut(120);
  tarch::parallel::Node::getInstance().setTimeOutWarning(60);
  #endif

  #if !defined(Asserts) && !defined(Debug)
  srand(time(0));
  #endif
}




particles::ExperimentSettingsMemoryBenchmark::ExperimentSettingsMemoryBenchmark() {
  runnerVariant = MemoryBenchmark;

  numberOfTimeSteps      = 50;
  maximalInitialVelocity = 1.0;

  timeSteppingAdapter    = ExplicitEuler;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();

  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( tarch::logging::CommandLineLogger::FilterListEntry( "info", false ) );

  #if defined(Debug)
  std::ostringstream traceFileName;
  traceFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.out";
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", false ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry(
    tarch::logging::CommandLineLogger::FilterListEntry( "debug", tarch::logging::CommandLineLogger::FilterListEntry::AnyRank,
    "particles::pit::mappings",
    false
  ));
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, traceFileName.str() );
  #elif defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
  #else
  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, false, "" );
  #endif
}






particles::ExperimentSettingsPureGrid::ExperimentSettingsPureGrid() {
  runnerVariant = Simulation;

  numberOfTimeSteps      = 50;
  maximalInitialVelocity = 1.0;
  timeSteppingAdapter    = PureGrid;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "", false ) );

  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "particles::pit::mappings", true ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peano::grid::RegularGridContainer", true ) );


  //  std::ostringstream logFileName;
  //  logFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.txt";
  //  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, logFileName.str() );
    #if defined(Asserts) || defined(Debug)
    tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
    #else
    tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, false, "" );
    #endif

  #if !defined(Asserts) && !defined(Debug)
  srand(time(0));
  #endif
}


particles::ExperimentSettingsStationaryParticles::ExperimentSettingsStationaryParticles() {
  runnerVariant = Simulation;

  numberOfTimeSteps      = 50;
  maximalInitialVelocity = 0.0;
  timeSteppingAdapter    = ExplicitEuler;

  tarch::logging::CommandLineLogger::getInstance().clearFilterList();
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "info", -1, "", false ) );

  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", true ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "particles::pit::mappings", true ) );
  tarch::logging::CommandLineLogger::getInstance().addFilterListEntry( ::tarch::logging::CommandLineLogger::FilterListEntry( "debug", -1, "peano::grid::RegularGridContainer", true ) );


  //  std::ostringstream logFileName;
  //  logFileName << "rank-" << tarch::parallel::Node::getInstance().getRank() << "-trace.txt";
  //  tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, logFileName.str() );
    #if defined(Asserts) || defined(Debug)
    tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, true, "" );
    #else
    tarch::logging::CommandLineLogger::getInstance().setLogFormat( " ", false, false, true, false, false, "" );
    #endif

  #if !defined(Asserts) && !defined(Debug)
  srand(time(0));
  #endif
}




const std::string particles::ExperimentSettingsRunTests::ID = "tests";
const std::string particles::ExperimentSettingsProduceVideo::ID = "video";
const std::string particles::ExperimentSettingsHomogeneousParticles::ID = "simulation";
const std::string particles::ExperimentSettingsMemoryBenchmark::ID = "mem-benchmark";
const std::string particles::ExperimentSettingsPureGrid::ID = "pure-grid";
const std::string particles::ExperimentSettingsStationaryParticles::ID = "stationary-particles";



