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
  particles::pit::State&  solverState,
  const double small_boundary /*=1e-7*/,
  const double big_boundary, /*=1e-3*/
  const int number_intervals /*=6*/)
  : small_boundary_(small_boundary),
    big_boundary_(big_boundary)
{
  assertion(small_boundary < big_boundary);
  state_ = solverState;
  if(number_intervals) {
    histogram_data_.resize(number_intervals,0);
  } else {
    histogram_data_.resize(log10(big_boundary)-log10(small_boundary)+2,0);
  }
}


particles::pit::myfunctions::Histogram::Histogram(
  particles::pit::State&  solverState,
  const int small_boundary /*=-7*/,
  const int big_boundary, /*=-3*/
  const int number_intervals /*=6*/)
  : small_boundary_(pow(10,small_boundary)),
    big_boundary_(pow(10,big_boundary))
{
  assertion(small_boundary < big_boundary);
  state_ = solverState;
  if(number_intervals) {
    histogram_data_.resize(number_intervals,0);
  } else {
    histogram_data_.resize(big_boundary-small_boundary+2,0);
  }
}


particles::pit::myfunctions::Histogram::~Histogram()
{
}


void particles::pit::myfunctions::Histogram::processHistogram(
  const tarch::la::Vector<DIMENSIONS,double>& norm
) {
  double check = small_boundary_;
  double oneStep = pow(big_boundary_/small_boundary_, 1.0/(histogram_data_.size()-2)); // We subtract 2 because two places already reserved

  // We save the data in increasing order.

  for(int d = 0; d<DIMENSIONS; d++) {
    check = small_boundary_;

    if(norm[d] < small_boundary_) {
        histogram_data_[0]++;
    } else if(norm[d] >= big_boundary_) {
        histogram_data_[histogram_data_.size() - 1]++;
    } else {

      for(int i = 1; i< histogram_data_.size()-1; i++) {
        check *= oneStep;

        if(norm[d] < check) {
            histogram_data_[i]++;
          break;
        }
      }
//std::cout << "CHECK:\n";
//      check = small_boundary_;
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


void particles::pit::myfunctions::Histogram::processHistogram10(
  const tarch::la::Vector<DIMENSIONS,double>& norm
) {
  double check = small_boundary_;
  double oneStep = 10; // We subtract 3 because two places already reserved
//  std::cout << "oneStep: " << oneStep << std::endl;
  //std::cout << "processHistogram(), oneStep=" << oneStep << std::endl;
  //std::cout << "processHistogram(), n_intervals=" << N_INTERVALS_HISTOGRAM << std::endl;

  // We save the data in increasing order.

  for(int d = 0; d<DIMENSIONS; d++) {
    check = small_boundary_;

    if(norm[d] < small_boundary_) {
        histogram_data_[0]++;
    } else if(norm[d] >= big_boundary_) {
        histogram_data_[histogram_data_.size() - 1]++;
    } else {

      for(int i = 1; i< histogram_data_.size()-1; i++) {
        check *= oneStep;

        if(norm[d] < check) {
            histogram_data_[i]++;
          break;
        }
      }
std::cout << "CHECK:\n";
      check = small_boundary_;
      for(int i = 1; i< histogram_data_.size()-1; i++) {
        check *= oneStep;
        std::cout << check << " ";

      }
std::cout << std::endl;
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
  std::ostringstream full_file_name;
  full_file_name << filename << "-" << state_.getMinimalNumberOfParticlesPerCell();
  std::ofstream out;
  if( writeFirstTime ) {
//    out.open( filename.c_str() );
    out.open( full_file_name.str().c_str() );
  } else {
//    out.open( filename.c_str(), std::ofstream::app );
    out.open( full_file_name.str().c_str(), std::ofstream::app );
  }
  if ( (!out.fail()) && out.is_open()) {
    if( writeFirstTime ) {
      // Information about boundaries
      out << small_boundary_ << " " << big_boundary_<< " "
        << state_.getMinimalNumberOfParticlesPerCell() << " "
        << state_.getMaximalInitialVelocity()/std::sqrt( static_cast<double>(DIMENSIONS) ) << " ";
      // Fill  till the end of the row to have the same number of elements in
      // each row (needed for matlab)
      for(int d = 4; d<histogram_data_.size(); d++) {
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
