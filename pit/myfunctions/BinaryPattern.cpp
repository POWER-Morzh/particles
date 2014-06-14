/*
 * BinaryPattern.cpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Denys Korzh
 */

#include "particles/pit/myfunctions/BinaryPattern.h"
#include <bitset>
#include <iostream>


void particles::pit::myfunctions::BinaryPattern::PrintBinaryDouble(const double& d){
    int tmp = 0;
    const char* double_ptr = reinterpret_cast<const char*> (&d) + 4;
    char* int_ptr = reinterpret_cast<char*> (&tmp);
    *int_ptr = *double_ptr;
    *(int_ptr + 1) = *(double_ptr + 1);
    *(int_ptr + 2) = *(double_ptr + 2);
    *(int_ptr + 3) = *(double_ptr + 3);
    std::bitset<32> bX(tmp);
    std::cout << bX << " ";
    int tmp2 = 0;
    const char* double_ptr2 = reinterpret_cast<const char*> (&d);
    char* int_ptr2 = reinterpret_cast<char*> (&tmp2);
    *int_ptr2 = *double_ptr2;
    *(int_ptr2 + 1) = *(double_ptr2 + 1);
    *(int_ptr2 + 2) = *(double_ptr2 + 2);
    *(int_ptr2 + 3) = *(double_ptr2 + 3);
    std::bitset<32> bX2(tmp2);
    std::cout << bX2 << " ";
}


void particles::pit::myfunctions::BinaryPattern::PrintBinaryDoubleCompressed(const double& d, const int mantissa) {
   short int _x;
   int tmp = 0;
   const char* double_ptr = reinterpret_cast<const char*> (&d) + 4;
   char* int_ptr = reinterpret_cast<char*> (&tmp);
   *int_ptr = *double_ptr;
   *(int_ptr + 1) = *(double_ptr + 1);
   *(int_ptr + 2) = *(double_ptr + 2);
   *(int_ptr + 3) = *(double_ptr + 3);
   if (tmp == 0) {
       _x = 0;
   } else {

   short int exponent = (tmp & 0x7ff00000) >> 20;

   tmp = tmp >> (20 - mantissa);


   _x = tmp & (0xffff >> (8+8-mantissa));
   exponent = exponent - 1023 + (0xfff >> (mantissa - 2));
   std::bitset<32> bExponent(exponent);

   _x |= (exponent << mantissa);
   _x |= (tmp >> (mantissa - 4)) & 0x8000;
   }
   std::bitset<16> b_X(_x);
   std::cout << b_X << std::endl;

}
