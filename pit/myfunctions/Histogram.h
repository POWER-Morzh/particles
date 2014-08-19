/*
 * Histogram.h
 *
 *  Created on: Aug 19, 2014
 *      Author: power-morzh
 */

#ifndef PARTICLES_PIT_MYFUNCTIONS_HISTOGRAM_H_
#define PARTICLES_PIT_MYFUNCTIONS_HISTOGRAM_H_

#include "peano/utils/Dimensions.h"
#include "tarch/la/Vector.h"

#define N_INTERVALS_HISTOGRAM 21

namespace particles {
  namespace pit {
    namespace myfunctions {
      class Histogram;
    } /* namespace myfunctions */
  } /* namespace pit */
} /* namespace particles */

class particles::pit::myfunctions::Histogram {
 public:
  Histogram(
    const double small_boundary=1e-7,
    const double big_boundary=1e-3);
  virtual ~Histogram();

  void processHistogram(const tarch::la::Vector<DIMENSIONS,double>& norm);
  void writeHistogramData(
    const std::string& filename,
    const bool writeFirstTime=true);
 private:
  const double small_boundary_;
  const double big_boundary_;
  tarch::la::Vector<N_INTERVALS_HISTOGRAM, int> histogram_data_;
};

#endif /* HISTOGRAM_H_ */
