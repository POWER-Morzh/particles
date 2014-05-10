#ifndef PARTICLES_PIT_MYFUNCTIONS_RepresentationChange_H_
#define PARTICLES_PIT_MYFUNCTIONS_RepresentationChange_H_

#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"

#include "peano/grid/VertexEnumerator.h"
#include "peano/MappingSpecification.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "particles/pit/Vertex.h"
#include "particles/pit/Cell.h"
#include "particles/pit/State.h"

#include "particles/pit/records/Particle.h"
#include "particles/pit/records/ParticleCompressed.h"


namespace particles {
  namespace pit {
    namespace myfunctions {
      class RepresentationChange;
    }
  }
}


class particles::pit::myfunctions::RepresentationChange {
  public:
	static double _globalMaxRelativeError;
	static double _globalMaxL2ErrorNorm;
	static int _globalNormAdditions;
	static tarch::la::Vector<DIMENSIONS, double> _globalL2ErrorNorm;
	static tarch::la::Vector<DIMENSIONS, double> _globalL2OffsetNorm;
	static int _iteration;

	static bool _outputInConsole;

	static std::ostringstream _maxRelativeErrorOut;
	static std::ostringstream _maxOffsetOut;
	static std::ostringstream _maxErrorOut;
	static std::ostringstream _RMSDOut;
	static std::ostringstream _L2NormOut;
	static std::ostringstream _L2ErrorNormOut;
	static std::ostringstream _meanVelocityOut;

	static void printParticlesInfo( const particles::pit::Cell& fineGridCell, const std::string normName, const tarch::la::Vector<DIMENSIONS, double> norm );

    static void leaveCell(particles::pit::Cell& fineGridCell);

    static tarch::la::Vector<DIMENSIONS, double> computeMeanVelocity( const ParticleHeap::HeapEntries& currentParticles, const int& size );

    /*
     * Compute norms
     */
    static tarch::la::Vector<DIMENSIONS, double> computeRMSD( const particles::pit::Cell& fineGridCell );
    static double computeMaxRelativeError( const particles::pit::Cell& fineGridCell );
    static double computeMaxOffset( const particles::pit::Cell& fineGridCell );
    static double computeMaxError( const particles::pit::Cell& fineGridCell );
    static tarch::la::Vector<DIMENSIONS, double> computeL2ErrorNorm( const particles::pit::Cell& fineGridCell );
    static tarch::la::Vector<DIMENSIONS, double> computeL2Norm( const particles::pit::Cell& fineGridCell );

    static void writeInCompressedHeap(
      const ParticleHeap::HeapEntries&                        currentParticles,
      const int                                               cellIndex,
      const tarch::la::Vector<DIMENSIONS, double>&            meanVelocity,
      const tarch::la::Vector<DIMENSIONS, double>&            meanCoordinate);

    /*
     * Skip value of globalL2Norm to 0 to start new iteration.
     */
    static void beginIteration();

    /*
     * Increment number of iterations and write all Data in file.
     */
    static void endIteration();

    /*
     * Write all Norms in files
     */
    static void writeAllInFile();

    /*
     * Write globalL2Norm in file, then it could be read from Matlab to plot it.
     */
    static void writeGlobalNorm( const std::string& filename, const tarch::la::Vector<DIMENSIONS,double>& Norm, const bool& writeFirstTime );

    /*
     * Write L2Norm in file, then it could be read from Matlab to plot it.
     */
    static void writeNorm( const std::string& normName, const std::ostringstream& normData );

    /*
     * Writes coordinates and values of norms for each cell in every iteration.
     */
    static void ascend(
      particles::pit::Cell * const    fineGridCells,
      particles::pit::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      particles::pit::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      particles::pit::Cell&           coarseGridCell
    );
};


#endif
