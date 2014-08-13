#include "particles/pit/records/ParticleCompressed.h"

#define MANTISSA 10
#if defined MANTISSA
#warning "MANTISSA DEFINED as " MANTISSA
#endif

particles::pit::records::ParticleCompressed::PersistentRecords::PersistentRecords() {
   
}


particles::pit::records::ParticleCompressed::ParticleCompressed() {
   
}


particles::pit::records::ParticleCompressed::ParticleCompressed(const PersistentRecords& persistentRecords):
_persistentRecords() {
   
}


particles::pit::records::ParticleCompressed::ParticleCompressed(const tarch::la::Vector<DIMENSIONS,double>& x, const tarch::la::Vector<DIMENSIONS,double>& v):
_x(x),
_v(v) {
   
}

particles::pit::records::ParticleCompressed::~ParticleCompressed() { }


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
 tarch::la::Vector<DIMENSIONS,double> particles::pit::records::ParticleCompressed::getX() const  {
   return _x;
}



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
 void particles::pit::records::ParticleCompressed::setX(const tarch::la::Vector<DIMENSIONS,double>& x)  {
   _x = (x);
}



 double particles::pit::records::ParticleCompressed::getX(int elementIndex) const  {
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   return _x[elementIndex];
   
}



 void particles::pit::records::ParticleCompressed::setX(int elementIndex, const double& x)  {
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   _x[elementIndex]= x;
   
}



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
 tarch::la::Vector<DIMENSIONS,double> particles::pit::records::ParticleCompressed::getV() const  {
   return _v;
}



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
 void particles::pit::records::ParticleCompressed::setV(const tarch::la::Vector<DIMENSIONS,double>& v)  {
   _v = (v);
}



 double particles::pit::records::ParticleCompressed::getV(int elementIndex) const  {
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   return _v[elementIndex];
   
}



 void particles::pit::records::ParticleCompressed::setV(int elementIndex, const double& v)  {
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   _v[elementIndex]= v;
   
}




std::string particles::pit::records::ParticleCompressed::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void particles::pit::records::ParticleCompressed::toString (std::ostream& out) const {
   out << "("; 
   out << "x:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getX(i) << ",";
   }
   out << getX(DIMENSIONS-1) << "]";
   out << ",";
   out << "v:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getV(i) << ",";
   }
   out << getV(DIMENSIONS-1) << "]";
   out <<  ")";
}


particles::pit::records::ParticleCompressed::PersistentRecords particles::pit::records::ParticleCompressed::getPersistentRecords() const {
   return _persistentRecords;
}

particles::pit::records::ParticleCompressedPacked particles::pit::records::ParticleCompressed::convert() const{
   return ParticleCompressedPacked(
      getX(),
      getV()
   );
}

#ifdef Parallel
   tarch::logging::Log particles::pit::records::ParticleCompressed::_log( "particles::pit::records::ParticleCompressed" );
   
   MPI_Datatype particles::pit::records::ParticleCompressed::Datatype = 0;
   MPI_Datatype particles::pit::records::ParticleCompressed::FullDatatype = 0;
   
   
   void particles::pit::records::ParticleCompressed::initDatatype() {
      {
         ParticleCompressed dummyParticleCompressed[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //x
            MPI_DOUBLE,		 //v
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            DIMENSIONS,		 //x
            DIMENSIONS,		 //v
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]._x[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]._v[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyParticleCompressed[1]._x[0])), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ParticleCompressed::Datatype );
         MPI_Type_commit( &ParticleCompressed::Datatype );
         
      }
      {
         ParticleCompressed dummyParticleCompressed[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_DOUBLE,		 //x
            MPI_DOUBLE,		 //v
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            DIMENSIONS,		 //x
            DIMENSIONS,		 //v
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]._x[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressed[0]._v[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyParticleCompressed[1]._x[0])), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ParticleCompressed::FullDatatype );
         MPI_Type_commit( &ParticleCompressed::FullDatatype );
         
      }
      
   }
   
   
   void particles::pit::records::ParticleCompressed::shutdownDatatype() {
      MPI_Type_free( &ParticleCompressed::Datatype );
      MPI_Type_free( &ParticleCompressed::FullDatatype );
      
   }
   
   void particles::pit::records::ParticleCompressed::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message particles::pit::records::ParticleCompressed "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Isend(
               this, 1, Datatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         else {
            result = MPI_Isend(
               this, 1, FullDatatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message particles::pit::records::ParticleCompressed "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished send task for particles::pit::records::ParticleCompressed "
               << toString()
               << " sent to node " << destination
               << " failed: " << tarch::parallel::MPIReturnValueToString(result);
               _log.error("send(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "particles::pit::records::ParticleCompressed",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::ParticleCompressed",
               "send(int)", destination,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         #ifdef Debug
         _log.debug("send(int,int)", "sent " + toString() );
         #endif
         
      }
      
   }
   
   
   
   void particles::pit::records::ParticleCompressed::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive particles::pit::records::ParticleCompressed from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Irecv(
               this, 1, Datatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         else {
            result = MPI_Irecv(
               this, 1, FullDatatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive particles::pit::records::ParticleCompressed from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished receive task for particles::pit::records::ParticleCompressed failed: "
               << tarch::parallel::MPIReturnValueToString(result);
               _log.error("receive(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "particles::pit::records::ParticleCompressed",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::ParticleCompressed",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool particles::pit::records::ParticleCompressed::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   
#endif


particles::pit::records::ParticleCompressedPacked::PersistentRecords::PersistentRecords() {
   
}


particles::pit::records::ParticleCompressedPacked::ParticleCompressedPacked() {
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   
}


particles::pit::records::ParticleCompressedPacked::ParticleCompressedPacked(const PersistentRecords& persistentRecords):
_persistentRecords() {
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   
}


particles::pit::records::ParticleCompressedPacked::ParticleCompressedPacked(const tarch::la::Vector<DIMENSIONS,double>& x, const tarch::la::Vector<DIMENSIONS,double>& v) {
   setX(x);
   setV(v);
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   assertion((sizeof(short int) >= 2));
   assertion((sizeof(double) == 8));
   
}

particles::pit::records::ParticleCompressedPacked::~ParticleCompressedPacked() { }


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
 tarch::la::Vector<DIMENSIONS,double> particles::pit::records::ParticleCompressedPacked::getX() const  {
   tarch::la::Vector<DIMENSIONS,double> result;
   for (int i = 0; i < DIMENSIONS; i++) {
      result[i] = getX(i);
   }
   return result;
}



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
 void particles::pit::records::ParticleCompressedPacked::setX(const tarch::la::Vector<DIMENSIONS,double>& x)  {
   for (int i = 0; i < DIMENSIONS; i++) {
      setX(i, x[i]);
   }
}



 double particles::pit::records::ParticleCompressedPacked::getX(int elementIndex) const  {
   
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);

   if( _x[elementIndex] == 0 ) {
        return 0;
   }
   int tmp = 0;
   tmp |= ((_x[elementIndex] & 0x8000) << 16);
   int exponent = ((_x[elementIndex] & (0x000fff >> (MANTISSA-3)) << MANTISSA) >> MANTISSA);
   exponent = exponent - (0xfff >> (MANTISSA - 2)) + 1023;
   tmp |= (exponent << 20);
   tmp |= ((_x[elementIndex] & (0xffff >> (16 -MANTISSA))) << (20 - MANTISSA));
   double result = 0;
   char* double_ptr = reinterpret_cast<char*> (&result) + 4;
   char* int_ptr = reinterpret_cast<char*> (&tmp);
   *double_ptr = *int_ptr;
   *(double_ptr + 1) = *(int_ptr + 1);
   *(double_ptr + 2) = *(int_ptr + 2);
   *(double_ptr + 3) = *(int_ptr + 3);
   
   return result;
   
}



 void particles::pit::records::ParticleCompressedPacked::setX(int elementIndex, const double& x)  {
   
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   int tmp = 0;
   const char* double_ptr = reinterpret_cast<const char*> (&x) + 4;
   char* int_ptr = reinterpret_cast<char*> (&tmp);
   *int_ptr = *double_ptr;
   *(int_ptr + 1) = *(double_ptr + 1);
   *(int_ptr + 2) = *(double_ptr + 2);
   *(int_ptr + 3) = *(double_ptr + 3);
   if( tmp == 0 ) {
      _x[elementIndex] = 0;
   } else {
      short int exponent = (tmp & 0x7ff00000) >> 20;
      tmp = tmp >> (20 - MANTISSA);
   
      _x[elementIndex] = tmp & (0xffff >> (8+8-MANTISSA));
      exponent = exponent - 1023 + (0xfff >> (MANTISSA - 2));
      _x[elementIndex] |= (exponent << MANTISSA);
      _x[elementIndex] |= (tmp >> (MANTISSA - 4)) & 0x8000;
   }
   
}



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
 tarch::la::Vector<DIMENSIONS,double> particles::pit::records::ParticleCompressedPacked::getV() const  {
   tarch::la::Vector<DIMENSIONS,double> result;
   for (int i = 0; i < DIMENSIONS; i++) {
      result[i] = getV(i);
   }
   return result;
}



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
 void particles::pit::records::ParticleCompressedPacked::setV(const tarch::la::Vector<DIMENSIONS,double>& v)  {
   for (int i = 0; i < DIMENSIONS; i++) {
      setV(i, v[i]);
   }
}



 double particles::pit::records::ParticleCompressedPacked::getV(int elementIndex) const  {
   
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);

   if( _v[elementIndex] == 0 ) {
	  return 0;
   }
   int tmp = 0;
   tmp |= ((_v[elementIndex] & 0x8000) << 16);
   int exponent = ((_v[elementIndex] & (0x000fff >> (MANTISSA-3)) << MANTISSA) >> MANTISSA);
   exponent = exponent - (0xfff >> (MANTISSA - 2)) + 1023;
   tmp |= (exponent << 20);
   tmp |= ((_v[elementIndex] & (0xffff >> (16 -MANTISSA))) << (20 - MANTISSA));
   double result = 0;
   char* double_ptr = reinterpret_cast<char*> (&result) + 4;
   char* int_ptr = reinterpret_cast<char*> (&tmp);
   *double_ptr = *int_ptr;
   *(double_ptr + 1) = *(int_ptr + 1);
   *(double_ptr + 2) = *(int_ptr + 2);
   *(double_ptr + 3) = *(int_ptr + 3);
   
   return result;
   
}



 void particles::pit::records::ParticleCompressedPacked::setV(int elementIndex, const double& v)  {
   
   assertion(elementIndex>=0);
   assertion(elementIndex<DIMENSIONS);
   int tmp = 0;
   const char* double_ptr = reinterpret_cast<const char*> (&v) + 4;
   char* int_ptr = reinterpret_cast<char*> (&tmp);
   *int_ptr = *double_ptr;
   *(int_ptr + 1) = *(double_ptr + 1);
   *(int_ptr + 2) = *(double_ptr + 2);
   *(int_ptr + 3) = *(double_ptr + 3);
   if( tmp == 0 ) {
      _v[elementIndex] = 0;
   } else {
      short int exponent = (tmp & 0x7ff00000) >> 20;
      tmp = tmp >> (20 - MANTISSA);
   
      _v[elementIndex] = tmp & (0xffff >> (8+8-MANTISSA));
      exponent = exponent - 1023 + (0xfff >> (MANTISSA - 2));
      _v[elementIndex] |= (exponent << MANTISSA);
      _v[elementIndex] |= (tmp >> (MANTISSA - 4)) & 0x8000;
   }
   
}




std::string particles::pit::records::ParticleCompressedPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void particles::pit::records::ParticleCompressedPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "x:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getX(i) << ",";
   }
   out << getX(DIMENSIONS-1) << "]";
   out << ",";
   out << "v:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getV(i) << ",";
   }
   out << getV(DIMENSIONS-1) << "]";
   out <<  ")";
}


particles::pit::records::ParticleCompressedPacked::PersistentRecords particles::pit::records::ParticleCompressedPacked::getPersistentRecords() const {
   return _persistentRecords;
}

particles::pit::records::ParticleCompressed particles::pit::records::ParticleCompressedPacked::convert() const{
   return ParticleCompressed(
      getX(),
      getV()
   );
}

#ifdef Parallel
   tarch::logging::Log particles::pit::records::ParticleCompressedPacked::_log( "particles::pit::records::ParticleCompressedPacked" );
   
   MPI_Datatype particles::pit::records::ParticleCompressedPacked::Datatype = 0;
   MPI_Datatype particles::pit::records::ParticleCompressedPacked::FullDatatype = 0;
   
   
   void particles::pit::records::ParticleCompressedPacked::initDatatype() {
      {
         ParticleCompressedPacked dummyParticleCompressedPacked[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //x
            MPI_SHORT,		 //v
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            DIMENSIONS,		 //x
            DIMENSIONS,		 //v
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]._x[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]._v[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyParticleCompressedPacked[1]._x[0])), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ParticleCompressedPacked::Datatype );
         MPI_Type_commit( &ParticleCompressedPacked::Datatype );
         
      }
      {
         ParticleCompressedPacked dummyParticleCompressedPacked[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //x
            MPI_SHORT,		 //v
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            DIMENSIONS,		 //x
            DIMENSIONS,		 //v
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]._x[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyParticleCompressedPacked[0]._v[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyParticleCompressedPacked[1]._x[0])), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ParticleCompressedPacked::FullDatatype );
         MPI_Type_commit( &ParticleCompressedPacked::FullDatatype );
         
      }
      
   }
   
   
   void particles::pit::records::ParticleCompressedPacked::shutdownDatatype() {
      MPI_Type_free( &ParticleCompressedPacked::Datatype );
      MPI_Type_free( &ParticleCompressedPacked::FullDatatype );
      
   }
   
   void particles::pit::records::ParticleCompressedPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message particles::pit::records::ParticleCompressedPacked "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Isend(
               this, 1, Datatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         else {
            result = MPI_Isend(
               this, 1, FullDatatype, destination,
               tag, tarch::parallel::Node::getInstance().getCommunicator(),
               sendRequestHandle
            );
            
         }
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message particles::pit::records::ParticleCompressedPacked "
            << toString()
            << " to node " << destination
            << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "send(int)",msg.str() );
         }
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished send task for particles::pit::records::ParticleCompressedPacked "
               << toString()
               << " sent to node " << destination
               << " failed: " << tarch::parallel::MPIReturnValueToString(result);
               _log.error("send(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "particles::pit::records::ParticleCompressedPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::ParticleCompressedPacked",
               "send(int)", destination,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         #ifdef Debug
         _log.debug("send(int,int)", "sent " + toString() );
         #endif
         
      }
      
   }
   
   
   
   void particles::pit::records::ParticleCompressedPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive particles::pit::records::ParticleCompressedPacked from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
      }
      else {
      
         MPI_Request* sendRequestHandle = new MPI_Request();
         MPI_Status   status;
         int          flag = 0;
         int          result;
         
         clock_t      timeOutWarning   = -1;
         clock_t      timeOutShutdown  = -1;
         bool         triggeredTimeoutWarning = false;
         
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            result = MPI_Irecv(
               this, 1, Datatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         else {
            result = MPI_Irecv(
               this, 1, FullDatatype, source, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
            );
            
         }
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive particles::pit::records::ParticleCompressedPacked from node "
            << source << ": " << tarch::parallel::MPIReturnValueToString(result);
            _log.error( "receive(int)", msg.str() );
         }
         
         result = MPI_Test( sendRequestHandle, &flag, &status );
         while (!flag) {
            if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
            if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
            result = MPI_Test( sendRequestHandle, &flag, &status );
            if (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "testing for finished receive task for particles::pit::records::ParticleCompressedPacked failed: "
               << tarch::parallel::MPIReturnValueToString(result);
               _log.error("receive(int)", msg.str() );
            }
            
            // deadlock aspect
            if (
               tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
               (clock()>timeOutWarning) &&
               (!triggeredTimeoutWarning)
            ) {
               tarch::parallel::Node::getInstance().writeTimeOutWarning(
               "particles::pit::records::ParticleCompressedPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::ParticleCompressedPacked",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool particles::pit::records::ParticleCompressedPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   
#endif



