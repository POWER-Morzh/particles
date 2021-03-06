#ifndef _PARTICLES_PIT_RECORDS_PARTICLECOMPRESSED_H
#define _PARTICLES_PIT_RECORDS_PARTICLECOMPRESSED_H

#include "peano/utils/Globals.h"
#include "tarch/compiler/CompilerSpecificSettings.h"
#include "peano/utils/PeanoOptimisations.h"
#ifdef Parallel
	#include "tarch/parallel/Node.h"
#endif
#ifdef Parallel
	#include <mpi.h>
#endif
#include "tarch/logging/Log.h"
#include "tarch/la/Vector.h"
#include <bitset>
#include <complex>
#include <string>
#include <iostream>

namespace particles {
   namespace pit {
      namespace records {
         class ParticleCompressed;
         class ParticleCompressedPacked;
      }
   }
}

/**
 * @author This class is generated by DaStGen
 * 		   DataStructureGenerator (DaStGen)
 * 		   2007-2009 Wolfgang Eckhardt
 * 		   2012      Tobias Weinzierl
 *
 * 		   build date: 09-02-2014 14:40
 *
 * @date   28/04/2014 17:38
 */
class particles::pit::records::ParticleCompressed { 
   
   public:
      
      typedef particles::pit::records::ParticleCompressedPacked Packed;
      
      struct PersistentRecords {
         /**
          * Generated
          */
         PersistentRecords();
         
         
      };
      
   private: 
      PersistentRecords _persistentRecords;
      public:

      #ifdef UseManualAlignment
      tarch::la::Vector<DIMENSIONS,double> _x __attribute__((aligned(VectorisationAlignment)));
      #else
      tarch::la::Vector<DIMENSIONS,double> _x;
      #endif
      private:

      public:

      #ifdef UseManualAlignment
      tarch::la::Vector<DIMENSIONS,double> _v __attribute__((aligned(VectorisationAlignment)));
      #else
      tarch::la::Vector<DIMENSIONS,double> _v;
      #endif
      private:

      
   public:
      /**
       * Generated
       */
      ParticleCompressed();
      
      /**
       * Generated
       */
      ParticleCompressed(const PersistentRecords& persistentRecords);
      
      /**
       * Generated
       */
      ParticleCompressed(const tarch::la::Vector<DIMENSIONS,double>& x, const tarch::la::Vector<DIMENSIONS,double>& v);
      
      /**
       * Generated
       */
      ~ParticleCompressed();
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       tarch::la::Vector<DIMENSIONS,double> getX() const ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       void setX(const tarch::la::Vector<DIMENSIONS,double>& x) ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       double getX(int elementIndex) const ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       void setX(int elementIndex, const double& x) ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       tarch::la::Vector<DIMENSIONS,double> getV() const ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       void setV(const tarch::la::Vector<DIMENSIONS,double>& v) ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       double getV(int elementIndex) const ;
      
      /**
       * Generated and optimized
       * 
       * If you realise a for loop using exclusively arrays (vectors) and compile 
       * with -DUseManualAlignment you may add 
       * \code
       #pragma vector aligned
       #pragma simd
       \endcode to this for loop to enforce your compiler to use SSE/AVX.
       * 
       * The alignment is tied to the unpacked records, i.e. for packed class
       * variants the machine's natural alignment is switched off to recude the  
       * memory footprint. Do not use any SSE/AVX operations or 
       * vectorisation on the result for the packed variants, as the data is misaligned. 
       * If you rely on vectorisation, convert the underlying record 
       * into the unpacked version first. 
       * 
       * @see convert()
       */
       void setV(int elementIndex, const double& v) ;
      
      /**
       * Generated
       */
      std::string toString() const;
      
      /**
       * Generated
       */
      void toString(std::ostream& out) const;
      
      
      PersistentRecords getPersistentRecords() const;
      /**
       * Generated
       */
      ParticleCompressedPacked convert() const;
      
      
   #ifdef Parallel
      protected:
         static tarch::logging::Log _log;
         
      public:
         
         /**
          * Global that represents the mpi datatype.
          * There are two variants: Datatype identifies only those attributes marked with
          * parallelise. FullDatatype instead identifies the whole record with all fields.
          */
         static MPI_Datatype Datatype;
         static MPI_Datatype FullDatatype;
         
         /**
          * Initializes the data type for the mpi operations. Has to be called
          * before the very first send or receive operation is called.
          */
         static void initDatatype();
         
         static void shutdownDatatype();
         
         void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
         
         static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
         
         #endif
            
         };
         
         /**
          * @author This class is generated by DaStGen
          * 		   DataStructureGenerator (DaStGen)
          * 		   2007-2009 Wolfgang Eckhardt
          * 		   2012      Tobias Weinzierl
          *
          * 		   build date: 09-02-2014 14:40
          *
          * @date   28/04/2014 17:38
          */
         class particles::pit::records::ParticleCompressedPacked { 
            
            public:
               
               struct PersistentRecords {
                  /**
                   * Generated
                   */
                  PersistentRecords();
                  
                  
               };
               
            private: 
               PersistentRecords _persistentRecords;
               short int _x[DIMENSIONS];
            short int _v[DIMENSIONS];
         
      public:
         /**
          * Generated
          */
         ParticleCompressedPacked();
         
         /**
          * Generated
          */
         ParticleCompressedPacked(const PersistentRecords& persistentRecords);
         
         /**
          * Generated
          */
         ParticleCompressedPacked(const tarch::la::Vector<DIMENSIONS,double>& x, const tarch::la::Vector<DIMENSIONS,double>& v);
         
         /**
          * Generated
          */
         ~ParticleCompressedPacked();
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          tarch::la::Vector<DIMENSIONS,double> getX() const ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          void setX(const tarch::la::Vector<DIMENSIONS,double>& x) ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          double getX(int elementIndex) const ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          void setX(int elementIndex, const double& x) ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          tarch::la::Vector<DIMENSIONS,double> getV() const ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          void setV(const tarch::la::Vector<DIMENSIONS,double>& v) ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          double getV(int elementIndex) const ;
         
         /**
          * Generated and optimized
          * 
          * If you realise a for loop using exclusively arrays (vectors) and compile 
          * with -DUseManualAlignment you may add 
          * \code
          #pragma vector aligned
          #pragma simd
          \endcode to this for loop to enforce your compiler to use SSE/AVX.
          * 
          * The alignment is tied to the unpacked records, i.e. for packed class
          * variants the machine's natural alignment is switched off to recude the  
          * memory footprint. Do not use any SSE/AVX operations or 
          * vectorisation on the result for the packed variants, as the data is misaligned. 
          * If you rely on vectorisation, convert the underlying record 
          * into the unpacked version first. 
          * 
          * @see convert()
          */
          void setV(int elementIndex, const double& v) ;
         
         /**
          * Generated
          */
         std::string toString() const;
         
         /**
          * Generated
          */
         void toString(std::ostream& out) const;
         
         
         PersistentRecords getPersistentRecords() const;
         /**
          * Generated
          */
         ParticleCompressed convert() const;
         
         
      #ifdef Parallel
         protected:
            static tarch::logging::Log _log;
            
         public:
            
            /**
             * Global that represents the mpi datatype.
             * There are two variants: Datatype identifies only those attributes marked with
             * parallelise. FullDatatype instead identifies the whole record with all fields.
             */
            static MPI_Datatype Datatype;
            static MPI_Datatype FullDatatype;
            
            /**
             * Initializes the data type for the mpi operations. Has to be called
             * before the very first send or receive operation is called.
             */
            static void initDatatype();
            
            static void shutdownDatatype();
            
            void send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
            
            void receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking);
            
            static bool isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise);
            
            #endif
               
            };
            
            #endif
            
