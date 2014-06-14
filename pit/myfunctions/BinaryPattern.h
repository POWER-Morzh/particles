/*
 * BinaryPattern.h
 *
 *  Created on: Jun 10, 2014
 *      Author: Denys Korzh
 */

#ifndef PARTICLES_PIT_MYFUNCTIONS_BINARYPATTERN_H_
#define PARTICLES_PIT_MYFUNCTIONS_BINARYPATTERN_H_

namespace particles {
  namespace pit {
    namespace myfunctions {
      class BinaryPattern;
    }
  }
}

class particles::pit::myfunctions::BinaryPattern {
 public:
  static void PrintBinaryDouble(const double& d);
  static void PrintBinaryDoubleCompressed(const double& d, const int mantissa);
};

#endif /* PARTICLES_PIT_MYFUNCTIONS_BINARYPATTERN_H_ */
