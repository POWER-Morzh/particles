#include "particles/pit/mappings/PlotParticles.h"


#include "particles/pit/myfunctions/RepresentationChange.h"


tarch::logging::Log                particles::pit::mappings::PlotParticles::_log( "particles::pit::mappings::PlotParticles" );
int  particles::pit::mappings::PlotParticles::_plotCounter(0);


peano::MappingSpecification   particles::pit::mappings::PlotParticles::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::PlotParticles::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::PlotParticles::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::OnlyLeaves,peano::MappingSpecification::Serial,false);
}


peano::MappingSpecification   particles::pit::mappings::PlotParticles::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   particles::pit::mappings::PlotParticles::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


peano::MappingSpecification   particles::pit::mappings::PlotParticles::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::RunConcurrentlyOnFineGrid,false);
}


particles::pit::mappings::PlotParticles::PlotParticles():
  _vtkWriter(0),
  _particleVelocityWriter(0) {
}


particles::pit::mappings::PlotParticles::~PlotParticles() {
}


#if defined(SharedMemoryParallelisation)
particles::pit::mappings::PlotParticles::PlotParticles(const PlotParticles&  masterThread):
  _vtkWriter(masterThread._vtkWriter),
  _particleVelocityWriter(masterThread._particleVelocityWriter) {
}


void particles::pit::mappings::PlotParticles::mergeWithWorkerThread(const PlotParticles& workerThread) {
}
#endif


void particles::pit::mappings::PlotParticles::createHangingVertex(
  particles::pit::Vertex&                      fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridH,
  particles::pit::Vertex * const               coarseGridVertices,
  const peano::grid::VertexEnumerator&         coarseGridVerticesEnumerator,
  particles::pit::Cell&                        coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&     fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::destroyHangingVertex(
  const particles::pit::Vertex&                fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&  fineGridH,
  particles::pit::Vertex * const               coarseGridVertices,
  const peano::grid::VertexEnumerator&         coarseGridVerticesEnumerator,
  particles::pit::Cell&                        coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&     fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::createInnerVertex(
      particles::pit::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::createBoundaryVertex(
      particles::pit::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::destroyVertex(
      const particles::pit::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::createCell(
      particles::pit::Cell&                 fineGridCell,
      particles::pit::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void particles::pit::mappings::PlotParticles::destroyCell(
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
void particles::pit::mappings::PlotParticles::mergeWithNeighbour(
  particles::pit::Vertex&  vertex,
  const particles::pit::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
}

void particles::pit::mappings::PlotParticles::prepareSendToNeighbour(
    particles::pit::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}

void particles::pit::mappings::PlotParticles::prepareCopyToRemoteNode(
    particles::pit::Vertex&  vertex,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}

void particles::pit::mappings::PlotParticles::prepareCopyToRemoteNode(
    particles::pit::Cell&  localCell,
    int                                           toRank,
    const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
    int                                           level
) {
}

void particles::pit::mappings::PlotParticles::mergeWithRemoteDataDueToForkOrJoin(
  particles::pit::Vertex&  localVertex,
  const particles::pit::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}

void particles::pit::mappings::PlotParticles::mergeWithRemoteDataDueToForkOrJoin(
  particles::pit::Cell&  localCell,
  const particles::pit::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}

bool particles::pit::mappings::PlotParticles::prepareSendToWorker(
  particles::pit::Cell&                 fineGridCell,
  particles::pit::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  particles::pit::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}

void particles::pit::mappings::PlotParticles::prepareSendToMaster(
    particles::pit::Cell&                       localCell,
    particles::pit::Vertex *                    vertices,
    const peano::grid::VertexEnumerator&       verticesEnumerator,
    const particles::pit::Vertex * const        coarseGridVertices,
    const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
    const particles::pit::Cell&                 coarseGridCell,
    const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void particles::pit::mappings::PlotParticles::mergeWithMaster(
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


void particles::pit::mappings::PlotParticles::receiveDataFromMaster(
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


void particles::pit::mappings::PlotParticles::mergeWithWorker(
    particles::pit::Cell&           localCell,
    const particles::pit::Cell&     receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
    const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
    int                                          level
) {
}


void particles::pit::mappings::PlotParticles::mergeWithWorker(
    particles::pit::Vertex&        localVertex,
    const particles::pit::Vertex&  receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS,double>&   x,
    const tarch::la::Vector<DIMENSIONS,double>&   h,
    int                                           level
) {
}
#endif

void particles::pit::mappings::PlotParticles::touchVertexFirstTime(
      particles::pit::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      particles::pit::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      particles::pit::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::touchVertexLastTime(
      particles::pit::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void particles::pit::mappings::PlotParticles::enterCell(
  particles::pit::Cell&                     fineGridCell,
  particles::pit::Vertex * const            fineGridVertices,
  const peano::grid::VertexEnumerator&      fineGridVerticesEnumerator,
  particles::pit::Vertex * const            coarseGridVertices,
  const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
  particles::pit::Cell&                     coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&  fineGridPositionOfCell
) {
  logTraceInWith4Arguments( "enterCell(...)", fineGridCell, fineGridVerticesEnumerator.toString(), coarseGridCell, fineGridPositionOfCell );

  const int cellIndex         = fineGridCell.getCellIndex();
  if (!fineGridCell.isRefined() && cellIndex>=0) {
    const int numberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
    for (int i=0; i<numberOfParticles; i++) {
      const particles::pit::records::Particle& currentParticle = ParticleHeap::getInstance().getData(cellIndex).at(i);
      const int particleNumber = _vtkWriter->plotPoint(currentParticle._persistentRecords._x);
      _particleVelocityWriter->plot(particleNumber,currentParticle._persistentRecords._v);
    }
  }

  logTraceOutWith1Argument( "enterCell(...)", fineGridCell );
}


void particles::pit::mappings::PlotParticles::leaveCell(
      particles::pit::Cell&           fineGridCell,
      particles::pit::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  const int cellIndex         = fineGridCell.getCellIndex();

  /*
   * We need to check if our Cell is a Leaf Cell and if there particles in this Cell before start the computation of meanVelocity.
   */
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  if(fineGridCell.isLeaf() && (NumberOfParticles > 0)) {
    particles::pit::myfunctions::RepresentationChange::leaveCell(fineGridCell);
  }
}


void particles::pit::mappings::PlotParticles::beginIteration(
  particles::pit::State&  solverState
) {
  logTraceInWith1Argument( "beginIteration(State)", solverState );

  assertion( _vtkWriter==0 );

  _vtkWriter = new tarch::plotter::pointdata::vtk::VTKTextFileWriter();

  _particleVelocityWriter = _vtkWriter->createPointDataWriter("Velocity", DIMENSIONS);

  particles::pit::myfunctions::RepresentationChange::beginIteration(solverState);

  logTraceOutWith1Argument( "beginIteration(State)", solverState);
}


void particles::pit::mappings::PlotParticles::endIteration(
  particles::pit::State&  solverState
) {
  logTraceInWith1Argument( "endIteration(State)", solverState );


  _vtkWriter->close();
  _particleVelocityWriter->close();

  delete _particleVelocityWriter;

  _particleVelocityWriter  = 0;

  std::ostringstream snapshotFileName;
  snapshotFileName << "particles"
                   #ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
                   #endif
                   << "-" << _plotCounter
                   << ".vtk";
  _vtkWriter->writeToFile( snapshotFileName.str() );

  _plotCounter++;

  delete _vtkWriter;
  _vtkWriter = 0;

  particles::pit::myfunctions::RepresentationChange::endIteration();

  logTraceOutWith1Argument( "endIteration(State)", solverState);
}



void particles::pit::mappings::PlotParticles::descend(
  particles::pit::Cell * const          fineGridCells,
  particles::pit::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  particles::pit::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  particles::pit::Cell&                 coarseGridCell
) {
}


void particles::pit::mappings::PlotParticles::ascend(
  particles::pit::Cell * const    fineGridCells,
  particles::pit::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  particles::pit::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  particles::pit::Cell&           coarseGridCell
) {
	particles::pit::myfunctions::RepresentationChange::ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
}
