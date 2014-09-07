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
#include "particles/pit/State.h"

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

  /**
   * If number_intervals==0 then it will create a histogram with distances of 10.
   */
  Histogram(
      particles::pit::State&  solverState,
    const double small_boundary=1e-7,
    const double big_boundary=1e-3,
    const int number_intervals=0);
  Histogram(
    particles::pit::State&  solverState,
    const int small_boundary=-7,
    const int big_boundary=-3,
    const int number_intervals=0);
  virtual ~Histogram();

  void processHistogram(const tarch::la::Vector<DIMENSIONS,double>& norm);
  void processHistogram10(const tarch::la::Vector<DIMENSIONS,double>& norm);
  void writeHistogramData(
    const std::string& filename,
    const bool writeFirstTime=true);
 private:
  const double small_boundary_;
  const double big_boundary_;
  State state_;
  std::vector<int> histogram_data_;
};

#endif /* HISTOGRAM_H_ */
