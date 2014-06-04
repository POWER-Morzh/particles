#include "particles/pit/myfunctions/RepresentationChange.h"
#include "particles/pit/myfunctions/CoordinatesRepresentationChange.h"
#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <string>
#include "peano/utils/Loop.h"

double particles::pit::myfunctions::CoordinatesRepresentationChange::_global_max_error = 0;
double particles::pit::myfunctions::CoordinatesRepresentationChange::_globalMaxRelativeError = 0;
double particles::pit::myfunctions::CoordinatesRepresentationChange::_globalMaxL2ErrorNorm = 0;
int particles::pit::myfunctions::CoordinatesRepresentationChange::_globalNormAdditions = 0;
tarch::la::Vector<DIMENSIONS, double> particles::pit::myfunctions::CoordinatesRepresentationChange::_globalL2ErrorNorm(0);

int particles::pit::myfunctions::CoordinatesRepresentationChange::_iteration = 0;

bool particles::pit::myfunctions::CoordinatesRepresentationChange::_outputInConsole = 0;

std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_maxRelativeErrorOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_maxErrorOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_maxOffsetOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_RMSDOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_L2ErrorNormOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_L2NormOut;
std::ostringstream particles::pit::myfunctions::CoordinatesRepresentationChange::_MeanCoordinateOut;


void particles::pit::myfunctions::CoordinatesRepresentationChange::printParticlesInfo(const particles::pit::Cell& fineGridCell, const std::string normName, const tarch::la::Vector<DIMENSIONS, double> norm ) {
  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();

  if(NumberOfParticles > 0) {
    const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
    const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);
    const tarch::la::Vector<DIMENSIONS, double> MeanCoordinate = fineGridCell.getMeanCoordinate();

    // Setuping view of output
    std::cout.setf( std::ios::fixed, std:: ios::floatfield );
    std::cout.precision(20);
    //typedef std::numeric_limits< double > dbl;
    //std::cout << "Precision for output: " << dbl::digits10 << std::endl;

    std::cout << "\n************ Data for Cell [" << cellIndex << "];" << " Number of Particles: " << NumberOfParticles << " **************\n";
    std::cout << "************ Mean coordinates saved in Cell: [";
    for(int d=0; d<DIMENSIONS; d++) {
      std::cout << MeanCoordinate[d] << " ";
    }
    std::cout << "]\n";

    // Norm
    std::cout << "\n////" << normName << ": ( " ;
    for(int d=0; d<DIMENSIONS; d++) {
      std::cout << norm[d] << " ";
    }
    std::cout << " );\n\n";

    // Coordinates
    std::cout << "---Coordinates---------------------------:" << std::endl;
    for (int i=0; i<NumberOfParticles; i++) {
      for(int d=0; d<DIMENSIONS; d++) {
    	std::cout << ( currentParticles.at(i)._persistentRecords._x(d) ) << " ";
      }
      std::cout << std::endl;
    }

    // Offsets of coordinates
    std::cout << "---Offsets of coordinates---------------------------:" << std::endl;
    for (int i=0; i<NumberOfParticles; i++) {
      for(int d=0; d<DIMENSIONS; d++) {
    	std::cout << ( currentParticles.at(i)._persistentRecords._x(d) - MeanCoordinate[d] ) << " ";
      }
      std::cout << std::endl;
    }

    // Compressed offsets of coordinates
    std::cout << "---Offsets of coordinates compressed---------------------------:" << std::endl;
    for(int d=0; d<DIMENSIONS; d++) {
      for (int i=0; i<NumberOfParticles; i++) {
    	std::cout << std::abs( compressedParticles.at(i).getX()[d] ) << " + ";
      }
      std::cout << std::endl;
    }

    // Error
    std::cout << "---Errors for coordinates---------------------------:" << std::endl;
    for(int d=0; d<DIMENSIONS; d++) {
      for (int i=0; i<NumberOfParticles; i++) {
    	std::cout << std::abs( (currentParticles.at(i)._persistentRecords._x[d] - MeanCoordinate[d] - compressedParticles.at(i).getX()[d]) ) << " + ";
      }
      std::cout << std::endl;
    }

    std::cout << "********End of output for Cell [" << cellIndex << "];" << " Number of Particles: " << NumberOfParticles << " ***********\n\n";
  }
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::myfunctions::CoordinatesRepresentationChange::computeMeanCoordinate(
  const ParticleHeap::HeapEntries& currentParticles,
  const int& NumberOfParticles
) {
  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate(0);
  for (int i=0; i<NumberOfParticles; i++) {
	MeanCoordinate += currentParticles.at(i)._persistentRecords._x;
  }
  for (int d =0; d<DIMENSIONS; d++) {
    MeanCoordinate[d] /= NumberOfParticles;
  }


  return MeanCoordinate;
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::myfunctions::CoordinatesRepresentationChange::computeRMSD( const particles::pit::Cell& fineGridCell ) {
  tarch::la::Vector<DIMENSIONS, double> rmsd(0);

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
  const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);

  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate(0);

  //Here we check if there more than one particle to prevent division by zero
  if ( NumberOfParticles > 1 ) {
    // Compute the mean value of the coordinate for each axes
    MeanCoordinate = computeMeanCoordinate(currentParticles, NumberOfParticles);

    double preciseOffset = 0;

    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	    preciseOffset = currentParticles.at(i)._persistentRecords._x[d] - MeanCoordinate[d];

	    rmsd[d] += std::pow ( ((preciseOffset - compressedParticles.at(i).getX()[d]) / preciseOffset), 2);
	  }
    }
    for(int d=0; d < DIMENSIONS; d++) {
      rmsd[d] = std::sqrt( rmsd[d] / NumberOfParticles );
    }
  }

  return rmsd;
}


double particles::pit::myfunctions::CoordinatesRepresentationChange::computeMaxRelativeError( const particles::pit::Cell& fineGridCell ) {

  double maxRelativeError = 0;

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
  const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);

  tarch::la::Vector<DIMENSIONS, double> meanCoordinate(0);

  // Here we check if there more than one particle to prevent division by zero
  if ( NumberOfParticles > 1 ) {
    // Compute the mean value of the coordinate for each axes
    meanCoordinate = computeMeanCoordinate(currentParticles, NumberOfParticles);

    double precise_offset = 0;

    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	      precise_offset = currentParticles.at(i)._persistentRecords._x[d] - meanCoordinate[d];

	    if( maxRelativeError <  std::abs( (precise_offset - compressedParticles.at(i).getX()[d]) / precise_offset ) ) {
	      maxRelativeError = std::abs( (precise_offset - compressedParticles.at(i).getX()[d]) / precise_offset);
	    }
	  }
    }
  }

  return maxRelativeError;
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::myfunctions::CoordinatesRepresentationChange::computeL2ErrorNorm( const particles::pit::Cell& fineGridCell ) {

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
  const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);

  tarch::la::Vector<DIMENSIONS, double> l2ErrorNorm(0);
  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate(0);

  if ( NumberOfParticles > 0 ) {
    // Compute the mean value of the coordinate for each axes
    MeanCoordinate = computeMeanCoordinate(currentParticles, NumberOfParticles);

    double precise_offset = 0;
    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	    precise_offset = currentParticles.at(i)._persistentRecords._x[d] - MeanCoordinate[d];
	    l2ErrorNorm[d] += std::abs(MUL_FACTOR*precise_offset - compressedParticles.at(i).getX()[d]) / MUL_FACTOR;
	  }
    }
    for (int d =0; d<DIMENSIONS; d++) {
      l2ErrorNorm[d] /= NumberOfParticles;
    }
  }


  return l2ErrorNorm;
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::myfunctions::CoordinatesRepresentationChange::computeL2Norm( const particles::pit::Cell& fineGridCell ) {

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
  const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);

  tarch::la::Vector<DIMENSIONS, double> l2Norm(0);
  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate(0);

  if ( NumberOfParticles > 0 ) {
    // Compute the mean value of the coordinate for each axes
    MeanCoordinate = computeMeanCoordinate(currentParticles, NumberOfParticles);


    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	    l2Norm[d] += std::abs(compressedParticles.at(i).getX()[d]);
	  }
    }
    for (int d =0; d<DIMENSIONS; d++) {
      l2Norm[d] /= NumberOfParticles;
    }
  }


  return l2Norm;
}


double particles::pit::myfunctions::CoordinatesRepresentationChange::computeMaxOffset( const particles::pit::Cell& fineGridCell ) {

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);

  double maxOffset = 0;
  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate = fineGridCell.getMeanCoordinate();

  if ( NumberOfParticles > 0 ) {

    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	    if( maxOffset < std::abs(currentParticles.at(i)._persistentRecords._x[d] - MeanCoordinate[d]) )
	      maxOffset = std::abs(currentParticles.at(i)._persistentRecords._x[d] - MeanCoordinate[d]);
	  }
    }

  }


  return maxOffset;
}


double particles::pit::myfunctions::CoordinatesRepresentationChange::computeMaxError( const particles::pit::Cell& fineGridCell ) {

  const int cellIndex = fineGridCell.getCellIndex();
  const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();
  const ParticleHeap::HeapEntries& currentParticles = ParticleHeap::getInstance().getData(cellIndex);
  const ParticleCompressedHeap::HeapEntries& compressedParticles = ParticleCompressedHeap::getInstance().getData(cellIndex);

  double maxOffset = 0;
  tarch::la::Vector<DIMENSIONS, double> MeanCoordinate = fineGridCell.getMeanCoordinate();

  if ( NumberOfParticles > 0 ) {
    double real_error = 0;

    for (int i=0; i<NumberOfParticles; i++) {
	  for(int d=0; d < DIMENSIONS; d++) {
	    real_error = MUL_FACTOR*(std::abs(currentParticles.at(i)._persistentRecords._x[d])
	                 - MeanCoordinate[d])
	                 - compressedParticles.at(i).getX()[d];
	    if( maxOffset < std::abs(real_error) / MUL_FACTOR )
	      maxOffset = std::abs(real_error) / MUL_FACTOR;
	  }
    }

  }


  return maxOffset;
}


void particles::pit::myfunctions::CoordinatesRepresentationChange::beginIteration() {
  _global_max_error = 0;
  _globalMaxRelativeError = 0;
  _globalMaxL2ErrorNorm = 0;
  tarch::la::Vector<DIMENSIONS, double> zeroVector(0);
  _globalL2ErrorNorm = zeroVector;
  _globalNormAdditions = 0;

  _maxRelativeErrorOut.str("");
  _maxRelativeErrorOut.clear();
  _maxErrorOut.str("");
  _maxErrorOut.clear();
  _maxOffsetOut.str("");
  _maxOffsetOut.clear();
  _RMSDOut.str("");
  _RMSDOut.clear();
  _L2ErrorNormOut.str("");
  _L2ErrorNormOut.clear();
  _L2NormOut.str("");
  _L2NormOut.clear();
  _MeanCoordinateOut.str("");
  _MeanCoordinateOut.clear();
}


void particles::pit::myfunctions::CoordinatesRepresentationChange::endIteration() {

    if ( _globalNormAdditions > 0 ) {
        for(int d = 0; d<DIMENSIONS; d++) {
            _globalL2ErrorNorm[d] /= _globalNormAdditions;
        }
    }

  writeAllInFile();

  ++_iteration;

}


void particles::pit::myfunctions::CoordinatesRepresentationChange::writeAllInFile() {
  // Write maxRelativeError
  //writeNorm( "CoordinatemaxRelativeError", _maxRelativeErrorOut );

  // Write maxError
  //writeNorm( "coordinatemaxError", _maxErrorOut );

  // Write maxOffset
  //writeNorm( "coordinatemaxOffset", _maxOffsetOut );

  // Write maxRelativeError
  //writeNorm( "coordinateRMSDOut", _RMSDOut );

  // Write l2ErrorNorm
  //writeNorm( "coordinateL2ErrorNorm", _L2ErrorNormOut );

  // Write l2Norm
  //writeNorm( "coordinateL2Norm", _L2NormOut );

  // Write MeanCoordinate
  //writeNorm( "coordinateMeanCoordinate", _MeanCoordinateOut );

  static bool writeFirstTime = 1;
  // Write globalL2Norm
  particles::pit::myfunctions::RepresentationChange::writeGlobalNorm( "coordinateGlobalL2ErrorNorm", _globalL2ErrorNorm, writeFirstTime );

  // Write _globalMaxL2ErrorNorm
  particles::pit::myfunctions::RepresentationChange::writeGlobalNorm( "coordinateGlobalMaxL2ErrorNorm.dat", _globalMaxL2ErrorNorm, writeFirstTime );

  // Write _globalMaxRelativeError
  //particles::pit::myfunctions::RepresentationChange::writeGlobalNorm( "coordinateGlobalMaxRelativeError", _globalMaxRelativeError, writeFirstTime );

  // Write _globalMaxRelativeError
  particles::pit::myfunctions::RepresentationChange::writeGlobalNorm( "coordinate_global_max_error", _global_max_error, writeFirstTime );
  writeFirstTime = 0;
}


void particles::pit::myfunctions::CoordinatesRepresentationChange::writeGlobalNorm( const std::string& filename, const tarch::la::Vector<DIMENSIONS,double>& Norm, const bool& writeFirstTime ) {
  std::ofstream out;
  if( writeFirstTime ) {
	out.open( filename.c_str() );
  } else {
    out.close();
    out.open( filename.c_str(), std::ofstream::app );
  }
  if ( (!out.fail()) && out.is_open() && !writeFirstTime ) {
	for(int d = 0; d<DIMENSIONS; d++) {
      out << Norm[d] << " ";
	}
	out << std::endl;
  }
}


void particles::pit::myfunctions::CoordinatesRepresentationChange::writeNorm( const std::string& normName, const std::ostringstream& normData ) {
  std::ostringstream normFileName;
  normFileName << normName
   		       << "-" << _iteration
    	       << ".dat";
  std::ofstream out;
  out.open( normFileName.str().c_str() );
  if ( (!out.fail()) && out.is_open() ) {
    out << normData.str() << std::endl;
  }
}


void particles::pit::myfunctions::CoordinatesRepresentationChange::ascend(
  particles::pit::Cell * const    fineGridCells,
  particles::pit::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  particles::pit::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  particles::pit::Cell&           coarseGridCell
) {

  dfor3(k)
	particles::pit::Cell fineGridCell = fineGridCells[ fineGridVerticesEnumerator.cell(k) ];
    const tarch::la::Vector<DIMENSIONS,double> cellOffset     = fineGridVerticesEnumerator.getVertexPosition(k);
    const tarch::la::Vector<DIMENSIONS,double> MeanCoordinate = fineGridCell.getMeanCoordinate();
    bool isLeaf = fineGridCell.isLeaf();
    const int cellIndex = fineGridCell.getCellIndex();
    const int NumberOfParticles = ParticleHeap::getInstance().getData(cellIndex).size();

    if( isLeaf && NumberOfParticles>1 ) {
      // Compute Max-Norm
      double maxRelativeError = computeMaxRelativeError( fineGridCell );
      // Save maximal maxRelativeError in _globalMaxRelativeError
      if (_globalMaxRelativeError < maxRelativeError) {
        _globalMaxRelativeError = maxRelativeError;
      }
      double maxError = computeMaxError( fineGridCell );
      if(_global_max_error < maxError) {
          _global_max_error = maxError;
      }
      double maxOffset = computeMaxOffset( fineGridCell );
      // Compute RMSD
      tarch::la::Vector<DIMENSIONS,double> rmsd = computeRMSD( fineGridCell );
      // Computer L2-Norm
      tarch::la::Vector<DIMENSIONS,double> l2ErrorNorm = computeL2ErrorNorm( fineGridCell );
      tarch::la::Vector<DIMENSIONS,double> l2Norm = computeL2Norm( fineGridCell );

      // Save maximal l2ErrorNorm in _globalMaxL2ErrorNorm
      for(int d = 0; d<DIMENSIONS; d++) {
        if(_globalMaxL2ErrorNorm < l2ErrorNorm[d]) {
          _globalMaxL2ErrorNorm = l2ErrorNorm[d];
        }
      }
      // Add l2Norm to _globalL2Norm
      _globalL2ErrorNorm += l2ErrorNorm;
      // Don't forget to increment _globalNormAdditions to divide _globalL2Norm by it
      // at the end of iteration before writing it in the file!
      ++_globalNormAdditions;
      std::cout << "_globalNormAdditions: " << _globalNormAdditions << std::endl;

      // Output for checking
      //printParticlesInfo( fineGridCell, "l2ErrorNorm", l2ErrorNorm );

      _maxRelativeErrorOut << maxRelativeError << " ";
      _maxErrorOut << maxError << " ";
      _maxOffsetOut << maxOffset << " ";
      for(int d=0; d<DIMENSIONS; d++) {
        _RMSDOut << rmsd[d] << " ";
        _L2ErrorNormOut << l2ErrorNorm[d] << " ";
        _L2NormOut << l2Norm[d] << " ";
        _MeanCoordinateOut << MeanCoordinate[d] << " ";
      }
      for(int d=0; d<DIMENSIONS; d++) {
        _maxRelativeErrorOut << cellOffset[d] << " ";
    	_maxErrorOut << cellOffset[d] << " ";
        _maxOffsetOut << cellOffset[d] << " ";
        _RMSDOut << cellOffset[d] << " ";
        _L2ErrorNormOut << cellOffset[d] << " ";
        _L2NormOut << cellOffset[d] << " ";
        _MeanCoordinateOut << cellOffset[d] << " ";
      }
      _maxRelativeErrorOut << std::endl;
      _maxErrorOut << std::endl;
      _maxOffsetOut << std::endl;
      _RMSDOut << std::endl;
      _L2ErrorNormOut << std::endl;
      _L2NormOut << std::endl;
      _MeanCoordinateOut << std::endl;
    }
  enddforx
}
