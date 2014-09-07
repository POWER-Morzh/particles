#include "particles/pit/mappings/RandomParticleDistribution.h"
#include "peano/grid/aspects/VertexStateAnalysis.h"
#include "tarch/multicore/MulticoreDefinitions.h"
#include <random>


//#include "particles/pit/myfunctions/RepresentationChange.h"


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::Serial,false);
}


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::RandomParticleDistribution::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


tarch::logging::Log  particles::pit::mappings::RandomParticleDistribution::_log( "particles::pit::mappings::RandomParticleDistribution" );


//
// Non-nop events
// ==============
//
void particles::pit::mappings::RandomParticleDistribution::initCell(
  particles::pit::Cell&                     fineGridCell,
  particles::pit::Vertex * const            fineGridVertices,
  const peano::grid::VertexEnumerator&      fineGridVerticesEnumerator
) {
  const int cellIndex = ParticleHeap::getInstance().createData();
  ParticleCompressedHeap::getInstance().createData();
  fineGridCell.init(cellIndex);

//std::cout << "cellIndex: " << cellIndex << " cellIndexCompressed: " << cellIndexCompressed << " cellNumber: " << fineGridCell.getCellIndex() << std::endl;

  const bool placeParticlesHere =
    peano::grid::aspects::VertexStateAnalysis::doAllVerticesCarryRefinementFlag(
      fineGridVertices,
      fineGridVerticesEnumerator,
      Vertex::Records::Unrefined
    )
    &&
    controlPointIsInsideParticleDomain(
      fineGridVerticesEnumerator.getCellCenter(),
      _state.getBreakingDamWidth()
    );

  const double velocityScaling = _state.getMaximalInitialVelocity() / std::sqrt( static_cast<double>(DIMENSIONS) );

  static std::default_random_engine generator;
  double mu = velocityScaling/2;
  double sigma = mu/3;
  static std::normal_distribution<double> distribution(mu,sigma);
  int velocity_sign;

  if (placeParticlesHere) {
    assertion(!fineGridCell.isRoot());

    int NumberOfParticles = getNumberOfParticles(
      fineGridVerticesEnumerator.getCellSize(),
      _state.getMaximalNumberOfInitialParticlesPerUnitVolume()
    );
    NumberOfParticles = NumberOfParticles==0 ? 1 : NumberOfParticles;

    logDebug( "enterCell(...)", "insert " << NumberOfParticles << " particle(s) into cell at " << fineGridVerticesEnumerator.toString() << " for density " << _state.getMaximalNumberOfInitialParticlesPerUnitVolume() );

    for (int i=0; i<NumberOfParticles; i++) {
      particles::pit::records::Particle newParticle;
      for (int d=0; d<DIMENSIONS; d++) {
        newParticle._persistentRecords._x(d) = static_cast<double>(rand())/static_cast<double>(RAND_MAX) * fineGridVerticesEnumerator.getCellSize()(d) + fineGridVerticesEnumerator.getVertexPosition()(d);
//        if(d != 0)
//        if(true)
//          newParticle._persistentRecords._v(d) = (static_cast<double>(rand())/static_cast<double>(RAND_MAX) * 2.0 - 1.0) * velocityScaling;
//        else
//          newParticle._persistentRecords._v(d) = (static_cast<double>(rand())/static_cast<double>(RAND_MAX) * 2.0 + 1.0) * velocityScaling;
        velocity_sign = rand() % 2 == 1 ? 1 : -1;
        newParticle._persistentRecords._v(d) = velocity_sign*distribution(generator);
      }
      ParticleHeap::getInstance().getData(cellIndex).push_back(newParticle);

    }

  }
}


#ifdef Parallel
bool particles::pit::mappings::RandomParticleDistribution::prepareSendToWorker(
  particles::pit::Cell&                 fineGridCell,
  particles::pit::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  particles::pit::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  initCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator);

  return false;
}
#endif


void particles::pit::mappings::RandomParticleDistribution::enterCell(
      particles::pit::Cell&                 fineGridCell,
      particles::pit::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  initCell(fineGridCell,fineGridVertices,fineGridVerticesEnumerator);

  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void particles::pit::mappings::RandomParticleDistribution::beginIteration(
  particles::pit::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  _state = solverState;

  ParticleHeap::getInstance().setName( "particles" );
  ParticleCompressedHeap::getInstance().setName( "particlesCompressed" );

  //particles::pit::myfunctions::RepresentationChange::beginIteration();

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


#if defined(SharedMemoryParallelisation)
particles::pit::mappings::RandomParticleDistribution::RandomParticleDistribution(const RandomParticleDistribution&  masterThread):
  _state( masterThread._state ) {
}


void particles::pit::mappings::RandomParticleDistribution::mergeWithWorkerThread(const RandomParticleDistribution& workerThread) {
}
#endif

//
// Nop
// ===
//
particles::pit::mappings::RandomParticleDistribution::RandomParticleDistribution() {
}


void particles::pit::mappings::RandomParticleDistribution::endIteration(
  particles::pit::State&  solverState
) {
  //particles::pit::myfunctions::RepresentationChange::endIteration();
}


particles::pit::mappings::RandomParticleDistribution::~RandomParticleDistribution() {
}


void particles::pit::mappings::RandomParticleDistribution::createCell(
  particles::pit::Cell&                 fineGridCell,
  particles::pit::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&  fineGridVerticesEnumerator,
  particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&  coarseGridVerticesEnumerator,
  particles::pit::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void particles::pit::mappings::RandomParticleDistribution::createHangingVertex(
      particles::pit::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      particles::pit::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      particles::pit::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
}


void particles::pit::mappings::RandomParticleDistribution::createInnerVertex(
      particles::pit::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      particles::pit::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      particles::pit::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
}


void particles::pit::mappings::RandomParticleDistribution::createBoundaryVertex(
      particles::pit::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      particles::pit::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      particles::pit::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
}


void particles::pit::mappings::RandomParticleDistribution::destroyHangingVertex(
  const particles::pit::Vertex&   fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
  particles::pit::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  particles::pit::Cell&           coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void particles::pit::mappings::RandomParticleDistribution::destroyVertex(
      const particles::pit::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}




void particles::pit::mappings::RandomParticleDistribution::destroyCell(
      const particles::pit::Cell&           fineGridCell,
      particles::pit::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}

#ifdef Parallel
void particles::pit::mappings::RandomParticleDistribution::mergeWithNeighbour(
  particles::pit::Vertex&  vertex,
  const particles::pit::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
}

void particles::pit::mappings::RandomParticleDistribution::prepareSendToNeighbour(
    particles::pit::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}

void particles::pit::mappings::RandomParticleDistribution::prepareCopyToRemoteNode(
    particles::pit::Vertex&  localVertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}

void particles::pit::mappings::RandomParticleDistribution::prepareCopyToRemoteNode(
    particles::pit::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
}

void particles::pit::mappings::RandomParticleDistribution::mergeWithRemoteDataDueToForkOrJoin(
  particles::pit::Vertex&  localVertex,
  const particles::pit::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}

void particles::pit::mappings::RandomParticleDistribution::mergeWithRemoteDataDueToForkOrJoin(
  particles::pit::Cell&  localCell,
  const particles::pit::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void particles::pit::mappings::RandomParticleDistribution::prepareSendToMaster(
  particles::pit::Cell&                       localCell,
  particles::pit::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator,
  const particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const particles::pit::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void particles::pit::mappings::RandomParticleDistribution::mergeWithMaster(
    const particles::pit::Cell&                       workerGridCell,
    particles::pit::Vertex * const                    workerGridVertices,
    const peano::grid::VertexEnumerator&             workerEnumerator,
    particles::pit::Cell&                             fineGridCell,
    particles::pit::Vertex * const                    fineGridVertices,
    const peano::grid::VertexEnumerator&             fineGridVerticesEnumerator,
    particles::pit::Vertex * const                    coarseGridVertices,
    const peano::grid::VertexEnumerator&             coarseGridVerticesEnumerator,
    particles::pit::Cell&                             coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&         fineGridPositionOfCell,
    int                                              worker,
    const particles::pit::State&                      workerState,
    particles::pit::State&                            masterState
) {
}


void particles::pit::mappings::RandomParticleDistribution::receiveDataFromMaster(
    particles::pit::Cell&                        receivedCell,
    particles::pit::Vertex *                     receivedVertices,
    const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
    particles::pit::Vertex * const               receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
    particles::pit::Cell&                        receivedCoarseGridCell,
    particles::pit::Vertex * const               workersCoarseGridVertices,
    const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
    particles::pit::Cell&                        workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
}


void particles::pit::mappings::RandomParticleDistribution::mergeWithWorker(
    particles::pit::Cell&           localCell,
    const particles::pit::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
}


void particles::pit::mappings::RandomParticleDistribution::mergeWithWorker(
    particles::pit::Vertex&        localVertex,
    const particles::pit::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}
#endif


void particles::pit::mappings::RandomParticleDistribution::touchVertexFirstTime(
      particles::pit::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void particles::pit::mappings::RandomParticleDistribution::touchVertexLastTime(
      particles::pit::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}

void particles::pit::mappings::RandomParticleDistribution::leaveCell(
      particles::pit::Cell&           fineGridCell,
      particles::pit::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
//  const int cellIndex         = fineGridCell.getCellIndex();
//
//  //std::cout << "leaveCell()-------" << cellIndex << std::endl;
//
//  /*
//   * We need to check if our Cell is a Leaf Cell and if there particles in this Cell before start the computation of meanVelocity.
//   */
//  const int size = ParticleHeap::getInstance().getData(cellIndex).size();
//  if(fineGridCell.isLeaf() && (size > 0)) {
//    particles::pit::myfunctions::RepresentationChange::changeRepresentation(fineGridCell);
//  }
}


void particles::pit::mappings::RandomParticleDistribution::descend(
  particles::pit::Cell * const          fineGridCells,
  particles::pit::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  particles::pit::Cell&                 coarseGridCell
) {
}


void particles::pit::mappings::RandomParticleDistribution::ascend(
  particles::pit::Cell * const    fineGridCells,
  particles::pit::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  particles::pit::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  particles::pit::Cell&           coarseGridCell
) {
}
