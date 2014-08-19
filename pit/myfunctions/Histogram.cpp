/*
 * Histogram.cpp
 *
 *  Created on: Aug 19, 2014
 *      Author: Denys Korzh
 */

#include <iostream>
#include <cmath>
#include <limits>
#include <fstream>
#include <string>

#include "tarch/Assertions.h"
#include "particles/pit/myfunctions/Histogram.h"


particles::pit::myfunctions::Histogram::Histogram(
  const double small_boundary /*=1e-7*/,
  const double big_boundary /*=1e-3*/)
  : small_boundary_(small_boundary),
    big_boundary_(big_boundary)
{
  tarch::la::Vector<N_INTERVALS_HISTOGRAM, int> zeroVector_N_INTERVALS_HISTOGRAM(0);
  histogram_data_ = zeroVector_N_INTERVALS_HISTOGRAM;
}

particles::pit::myfunctions::Histogram::~Histogram()
{
  // TODO Auto-generated destructor stub
}


void particles::pit::myfunctions::Histogram::processHistogram(
  const tarch::la::Vector<DIMENSIONS,double>& norm
) {
  assertion(small_bound < big_bound);
  double check = small_boundary_;
  double oneStep = pow(big_boundary_/small_boundary_, 1.0/(N_INTERVALS_HISTOGRAM-3)); // We subtract 3 because two places already reserved
//  std::cout << "oneStep: " << oneStep << std::endl;
  //std::cout << "processHistogram(), oneStep=" << oneStep << std::endl;
  //std::cout << "processHistogram(), n_intervals=" << N_INTERVALS_HISTOGRAM << std::endl;

  // We save the data in increasing order.

  for(int d = 0; d<DIMENSIONS; d++) {
    check = small_boundary_;

    if(norm[d] < small_boundary_) {
        histogram_data_[0]++;
    } else if(norm[d] >= big_boundary_) {
        histogram_data_[N_INTERVALS_HISTOGRAM - 1]++;
    } else {

      for(int i = 1; i< N_INTERVALS_HISTOGRAM-1; i++) {
        check *= oneStep;

        if(norm[d] < check) {
            histogram_data_[i]++;
          break;
        }
      }
//std::cout << "CHECK:\n";
//      check = small_bound;
//      for(int i = 1; i< N_INTERVALS_HISTOGRAM-1; i++) {
//        check *= oneStep;
//        std::cout << check << " ";
//
//      }
//std::cout << std::endl;
      if(norm[d] >= check && norm[d] <big_boundary_) {
        // If we came here it means that we did not save data
        std::cout << "processHistogram() - ERROR -  we didn't save data "
          "about Norm in _histogramData for: " << norm[d] << std::endl;
      }

    }
  }
}


void particles::pit::myfunctions::Histogram::writeHistogramData(
  const std::string& filename,
  const bool writeFirstTime /*=true*/
) {
  std::ofstream out;
  if( writeFirstTime ) {
    out.open( filename.c_str() );
  } else {
    out.open( filename.c_str(), std::ofstream::app );
  }
  if ( (!out.fail()) && out.is_open()) {
    if( writeFirstTime ) {
      // Information about boundarys
      out << small_boundary_ << " " << big_boundary_ << " ";
      // Fill  till the end of the row to have the same number of elements in
      // each row (needed for matlab)
      for(int d = 2; d<histogram_data_.size(); d++) {
        out << histogram_data_[d] << " ";
      }
      out << std::endl;
    } else {
      for(int d = 0; d<histogram_data_.size(); d++) {
        out << histogram_data_[d] << " ";
      }
      out << std::endl;
    }
  }
  out.close();
}
