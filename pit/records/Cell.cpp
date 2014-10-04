#include "particles/pit/records/Cell.h"

#if !defined(Debug) && !defined(Parallel) && defined(SharedMemoryParallelisation)
   particles::pit::records::Cell::PersistentRecords::PersistentRecords() {
      
   }
   
   
   particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
   _cellIndex(cellIndex),
   _meanVelocity(meanVelocity),
   _meanCoordinate(meanCoordinate),
   _myNorm(myNorm),
   _isInside(isInside),
   _state(state),
   _evenFlags(evenFlags),
   _accessNumber(accessNumber),
   _numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
   _numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
      
   }
   
   particles::pit::records::Cell::Cell() {
      
   }
   
   
   particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
      
   }
   
   
   particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
   _persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
      
   }
   
   
   particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
   _persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
      
   }
   
   particles::pit::records::Cell::~Cell() { }
   
   std::string particles::pit::records::Cell::toString(const State& param) {
      switch (param) {
         case Leaf: return "Leaf";
         case Refined: return "Refined";
         case Root: return "Root";
      }
      return "undefined";
   }
   
   std::string particles::pit::records::Cell::getStateMapping() {
      return "State(Leaf=0,Refined=1,Root=2)";
   }
   
   
   std::string particles::pit::records::Cell::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void particles::pit::records::Cell::toString (std::ostream& out) const {
      out << "("; 
      out << "cellIndex:" << getCellIndex();
      out << ",";
      out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
      out << ",";
      out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
      out << ",";
      out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
      out << ",";
      out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
      out << ",";
      out << "isInside:" << getIsInside();
      out << ",";
      out << "state:" << toString(getState());
      out << ",";
      out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
      out << ",";
      out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
      out << ",";
      out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
      out << ",";
      out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
      out <<  ")";
   }
   
   
   particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
      return CellPacked(
         getCellIndex(),
         getNumberOfParticlesInChildren(),
         getMeanVelocity(),
         getMeanCoordinate(),
         getMyNorm(),
         getIsInside(),
         getState(),
         getEvenFlags(),
         getAccessNumber(),
         getNumberOfLoadsFromInputStream(),
         getNumberOfStoresToOutputStream()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );
      
      MPI_Datatype particles::pit::records::Cell::Datatype = 0;
      MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;
      
      
      void particles::pit::records::Cell::initDatatype() {
         {
            Cell dummyCell[2];
            
            const int Attributes = 4;
            MPI_Datatype subtypes[Attributes] = {
               MPI_CHAR,		 //isInside
               MPI_INT,		 //state
               MPI_INT,		 //numberOfParticlesInChildren
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //isInside
               1,		 //state
               1,		 //numberOfParticlesInChildren
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[3] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
            MPI_Type_commit( &Cell::Datatype );
            
         }
         {
            Cell dummyCell[2];
            
            const int Attributes = 12;
            MPI_Datatype subtypes[Attributes] = {
               MPI_INT,		 //cellIndex
               MPI_DOUBLE,		 //meanVelocity
               MPI_DOUBLE,		 //meanCoordinate
               MPI_DOUBLE,		 //myNorm
               MPI_CHAR,		 //isInside
               MPI_INT,		 //state
               MPI_INT,		 //evenFlags
               MPI_SHORT,		 //accessNumber
               MPI_INT,		 //numberOfLoadsFromInputStream
               MPI_INT,		 //numberOfStoresToOutputStream
               MPI_INT,		 //numberOfParticlesInChildren
               MPI_UB		 // end/displacement flag
            };
            
            int blocklen[Attributes] = {
               1,		 //cellIndex
               DIMENSIONS,		 //meanVelocity
               DIMENSIONS,		 //meanCoordinate
               DIMENSIONS,		 //myNorm
               1,		 //isInside
               1,		 //state
               DIMENSIONS,		 //evenFlags
               DIMENSIONS_TIMES_TWO,		 //accessNumber
               1,		 //numberOfLoadsFromInputStream
               1,		 //numberOfStoresToOutputStream
               1,		 //numberOfParticlesInChildren
               1		 // end/displacement flag
            };
            
            MPI_Aint     disp[Attributes];
            
            MPI_Aint base;
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[6] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[7] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[8] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[9] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[10] );
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[11] );
            
            for (int i=1; i<Attributes; i++) {
               assertion1( disp[i] > disp[i-1], i );
            }
            for (int i=0; i<Attributes; i++) {
               disp[i] -= base;
            }
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
            MPI_Type_commit( &Cell::FullDatatype );
            
         }
         
      }
      
      
      void particles::pit::records::Cell::shutdownDatatype() {
         MPI_Type_free( &Cell::Datatype );
         MPI_Type_free( &Cell::FullDatatype );
         
      }
      
      void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
         _senderDestinationRank = destination;
         
         if (communicateBlocking) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message particles::pit::records::Cell "
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
            msg << "was not able to send message particles::pit::records::Cell "
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
               msg << "testing for finished send task for particles::pit::records::Cell "
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
               "particles::pit::records::Cell",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::Cell",
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
   
   
   
   void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      if (communicateBlocking) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         _senderDestinationRank = status.MPI_SOURCE;
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive particles::pit::records::Cell from node "
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
            msg << "failed to start to receive particles::pit::records::Cell from node "
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
               msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
               "particles::pit::records::Cell",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "particles::pit::records::Cell",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
         }
         
         delete sendRequestHandle;
         
         _senderDestinationRank = status.MPI_SOURCE;
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int particles::pit::records::Cell::getSenderRank() const {
      assertion( _senderDestinationRank!=-1 );
      return _senderDestinationRank;
      
   }
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
   setIsInside(isInside);
   setState(state);
   setEvenFlags(evenFlags);
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}

particles::pit::records::CellPacked::CellPacked() {
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
   assertion((DIMENSIONS+3 < (8 * sizeof(short int))));
   
}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
   return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
   return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "cellIndex:" << getCellIndex();
   out << ",";
   out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
   out << ",";
   out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
   out << ",";
   out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
   out << ",";
   out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
   out << ",";
   out << "isInside:" << getIsInside();
   out << ",";
   out << "state:" << toString(getState());
   out << ",";
   out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
   out << ",";
   out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
   out << ",";
   out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
   out << ",";
   out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
   out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
   return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
   return Cell(
      getCellIndex(),
      getNumberOfParticlesInChildren(),
      getMeanVelocity(),
      getMeanCoordinate(),
      getMyNorm(),
      getIsInside(),
      getState(),
      getEvenFlags(),
      getAccessNumber(),
      getNumberOfLoadsFromInputStream(),
      getNumberOfStoresToOutputStream()
   );
}

#ifdef Parallel
   tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );
   
   MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
   MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;
   
   
   void particles::pit::records::CellPacked::initDatatype() {
      {
         CellPacked dummyCellPacked[2];
         
         const int Attributes = 3;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //_packedRecords0
            MPI_INT,		 //numberOfParticlesInChildren
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //_packedRecords0
            1,		 //numberOfParticlesInChildren
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._packedRecords0))), 		&disp[2] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
         MPI_Type_commit( &CellPacked::Datatype );
         
      }
      {
         CellPacked dummyCellPacked[2];
         
         const int Attributes = 10;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //cellIndex
            MPI_DOUBLE,		 //meanVelocity
            MPI_DOUBLE,		 //meanCoordinate
            MPI_DOUBLE,		 //myNorm
            MPI_SHORT,		 //accessNumber
            MPI_INT,		 //numberOfLoadsFromInputStream
            MPI_INT,		 //numberOfStoresToOutputStream
            MPI_SHORT,		 //_packedRecords0
            MPI_INT,		 //numberOfParticlesInChildren
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //cellIndex
            DIMENSIONS,		 //meanVelocity
            DIMENSIONS,		 //meanCoordinate
            DIMENSIONS,		 //myNorm
            DIMENSIONS_TIMES_TWO,		 //accessNumber
            1,		 //numberOfLoadsFromInputStream
            1,		 //numberOfStoresToOutputStream
            1,		 //_packedRecords0
            1,		 //numberOfParticlesInChildren
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[9] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
         MPI_Type_commit( &CellPacked::FullDatatype );
         
      }
      
   }
   
   
   void particles::pit::records::CellPacked::shutdownDatatype() {
      MPI_Type_free( &CellPacked::Datatype );
      MPI_Type_free( &CellPacked::FullDatatype );
      
   }
   
   void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
      _senderDestinationRank = destination;
      
      if (communicateBlocking) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message particles::pit::records::CellPacked "
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
         msg << "was not able to send message particles::pit::records::CellPacked "
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
            msg << "testing for finished send task for particles::pit::records::CellPacked "
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
            "particles::pit::records::CellPacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
   if (communicateBlocking) {
   
      MPI_Status  status;
      const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
      _senderDestinationRank = status.MPI_SOURCE;
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
         msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
            msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
            "particles::pit::records::CellPacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "particles::pit::records::CellPacked",
            "receive(int)", source,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
      }
      
      delete sendRequestHandle;
      
      _senderDestinationRank = status.MPI_SOURCE;
      #ifdef Debug
      _log.debug("receive(int,int)", "received " + toString() ); 
      #endif
      
   }
   
}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
   assertion( _senderDestinationRank!=-1 );
   return _senderDestinationRank;
   
}
#endif



#elif !defined(Parallel) && defined(Debug) && !defined(SharedMemoryParallelisation)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
   case Leaf: return "Leaf";
   case Refined: return "Refined";
   case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
   getCellIndex(),
   getNumberOfParticlesInChildren(),
   getMeanVelocity(),
   getMeanCoordinate(),
   getMyNorm(),
   getIsInside(),
   getState(),
   getLevel(),
   getEvenFlags(),
   getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
   {
      Cell dummyCell[2];
      
      const int Attributes = 5;
      MPI_Datatype subtypes[Attributes] = {
         MPI_CHAR,		 //isInside
         MPI_INT,		 //state
         MPI_INT,		 //level
         MPI_INT,		 //numberOfParticlesInChildren
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //isInside
         1,		 //state
         1,		 //level
         1,		 //numberOfParticlesInChildren
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[4] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
      MPI_Type_commit( &Cell::Datatype );
      
   }
   {
      Cell dummyCell[2];
      
      const int Attributes = 11;
      MPI_Datatype subtypes[Attributes] = {
         MPI_INT,		 //cellIndex
         MPI_DOUBLE,		 //meanVelocity
         MPI_DOUBLE,		 //meanCoordinate
         MPI_DOUBLE,		 //myNorm
         MPI_CHAR,		 //isInside
         MPI_INT,		 //state
         MPI_INT,		 //level
         MPI_INT,		 //evenFlags
         MPI_SHORT,		 //accessNumber
         MPI_INT,		 //numberOfParticlesInChildren
         MPI_UB		 // end/displacement flag
      };
      
      int blocklen[Attributes] = {
         1,		 //cellIndex
         DIMENSIONS,		 //meanVelocity
         DIMENSIONS,		 //meanCoordinate
         DIMENSIONS,		 //myNorm
         1,		 //isInside
         1,		 //state
         1,		 //level
         DIMENSIONS,		 //evenFlags
         DIMENSIONS_TIMES_TWO,		 //accessNumber
         1,		 //numberOfParticlesInChildren
         1		 // end/displacement flag
      };
      
      MPI_Aint     disp[Attributes];
      
      MPI_Aint base;
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[6] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[7] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[8] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[9] );
      MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[10] );
      
      for (int i=1; i<Attributes; i++) {
         assertion1( disp[i] > disp[i-1], i );
      }
      for (int i=0; i<Attributes; i++) {
         disp[i] -= base;
      }
      MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
      MPI_Type_commit( &Cell::FullDatatype );
      
   }
   
}


void particles::pit::records::Cell::shutdownDatatype() {
   MPI_Type_free( &Cell::Datatype );
   MPI_Type_free( &Cell::FullDatatype );
   
}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
   _senderDestinationRank = destination;
   
   if (communicateBlocking) {
   
      const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message particles::pit::records::Cell "
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
      msg << "was not able to send message particles::pit::records::Cell "
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
         msg << "testing for finished send task for particles::pit::records::Cell "
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
         "particles::pit::records::Cell",
         "send(int)", destination,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

   MPI_Status  status;
   const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
   _senderDestinationRank = status.MPI_SOURCE;
   if ( result != MPI_SUCCESS ) {
      std::ostringstream msg;
      msg << "failed to start to receive particles::pit::records::Cell from node "
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
      msg << "failed to start to receive particles::pit::records::Cell from node "
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
         msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
         "particles::pit::records::Cell",
         "receive(int)", source,tag,1
         );
         triggeredTimeoutWarning = true;
      }
      if (
         tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
         (clock()>timeOutShutdown)
      ) {
         tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
         "particles::pit::records::Cell",
         "receive(int)", source,tag,1
         );
      }
      tarch::parallel::Node::getInstance().receiveDanglingMessages();
   }
   
   delete sendRequestHandle;
   
   _senderDestinationRank = status.MPI_SOURCE;
   #ifdef Debug
   _log.debug("receive(int,int)", "received " + toString() ); 
   #endif
   
}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_level(level),
_accessNumber(accessNumber) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
   CellPacked dummyCellPacked[2];
   
   const int Attributes = 4;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //level
      MPI_SHORT,		 //_packedRecords0
      MPI_INT,		 //numberOfParticlesInChildren
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //level
      1,		 //_packedRecords0
      1,		 //numberOfParticlesInChildren
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._level))), 		&disp[3] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
   MPI_Type_commit( &CellPacked::Datatype );
   
}
{
   CellPacked dummyCellPacked[2];
   
   const int Attributes = 9;
   MPI_Datatype subtypes[Attributes] = {
      MPI_INT,		 //cellIndex
      MPI_DOUBLE,		 //meanVelocity
      MPI_DOUBLE,		 //meanCoordinate
      MPI_DOUBLE,		 //myNorm
      MPI_INT,		 //level
      MPI_SHORT,		 //accessNumber
      MPI_SHORT,		 //_packedRecords0
      MPI_INT,		 //numberOfParticlesInChildren
      MPI_UB		 // end/displacement flag
   };
   
   int blocklen[Attributes] = {
      1,		 //cellIndex
      DIMENSIONS,		 //meanVelocity
      DIMENSIONS,		 //meanCoordinate
      DIMENSIONS,		 //myNorm
      1,		 //level
      DIMENSIONS_TIMES_TWO,		 //accessNumber
      1,		 //_packedRecords0
      1,		 //numberOfParticlesInChildren
      1		 // end/displacement flag
   };
   
   MPI_Aint     disp[Attributes];
   
   MPI_Aint base;
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[4] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[6] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[7] );
   MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[8] );
   
   for (int i=1; i<Attributes; i++) {
      assertion1( disp[i] > disp[i-1], i );
   }
   for (int i=0; i<Attributes; i++) {
      disp[i] -= base;
   }
   MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
   MPI_Type_commit( &CellPacked::FullDatatype );
   
}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

   const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
   if  (result!=MPI_SUCCESS) {
      std::ostringstream msg;
      msg << "was not able to send message particles::pit::records::CellPacked "
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
   msg << "was not able to send message particles::pit::records::CellPacked "
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
      msg << "testing for finished send task for particles::pit::records::CellPacked "
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
      "particles::pit::records::CellPacked",
      "send(int)", destination,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
   std::ostringstream msg;
   msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
   msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
      msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
      "particles::pit::records::CellPacked",
      "receive(int)", source,tag,1
      );
      triggeredTimeoutWarning = true;
   }
   if (
      tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
      (clock()>timeOutShutdown)
   ) {
      tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
      "particles::pit::records::CellPacked",
      "receive(int)", source,tag,1
      );
   }
   tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && !defined(Debug) && !defined(SharedMemoryParallelisation)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 18;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[17] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate()) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 9;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._subtreeHoldsWorker))), 		&disp[8] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 15;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[14] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif !defined(Debug) && !defined(Parallel) && !defined(SharedMemoryParallelisation)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 4;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[3] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_accessNumber(accessNumber) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 3;
MPI_Datatype subtypes[Attributes] = {
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._packedRecords0))), 		&disp[2] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 8;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_SHORT,		 //accessNumber
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[7] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && defined(SharedMemoryParallelisation) && defined(Debug)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 21;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[18] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[19] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[20] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_level(level),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate(), persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._level))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 18;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[17] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && defined(Debug) && !defined(SharedMemoryParallelisation)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 19;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[18] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_level(level),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate()) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //level
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //level
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._level))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 16;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[15] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif defined(Parallel) && !defined(Debug) && defined(SharedMemoryParallelisation)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_cellIsAForkCandidate(cellIsAForkCandidate),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords._cellIsAForkCandidate, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 10;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[9] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 20;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_CHAR,		 //cellIsAForkCandidate
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //cellIsAForkCandidate
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._responsibleRank))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._nodeWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._localWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._totalWorkload))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._maxWorkload))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._minWorkload))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIsAForkCandidate))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[16] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[17] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[18] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[19] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_accessNumber(accessNumber),
_responsibleRank(responsibleRank),
_subtreeHoldsWorker(subtreeHoldsWorker),
_nodeWorkload(nodeWorkload),
_localWorkload(localWorkload),
_totalWorkload(totalWorkload),
_maxWorkload(maxWorkload),
_minWorkload(minWorkload),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
setCellIsAForkCandidate(cellIsAForkCandidate);
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._responsibleRank, persistentRecords._subtreeHoldsWorker, persistentRecords._nodeWorkload, persistentRecords._localWorkload, persistentRecords._totalWorkload, persistentRecords._maxWorkload, persistentRecords._minWorkload, persistentRecords.getCellIsAForkCandidate(), persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& responsibleRank, const bool& subtreeHoldsWorker, const double& nodeWorkload, const double& localWorkload, const double& totalWorkload, const double& maxWorkload, const double& minWorkload, const bool& cellIsAForkCandidate, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, evenFlags, accessNumber, responsibleRank, subtreeHoldsWorker, nodeWorkload, localWorkload, totalWorkload, maxWorkload, minWorkload, cellIsAForkCandidate, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+4 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "responsibleRank:" << getResponsibleRank();
out << ",";
out << "subtreeHoldsWorker:" << getSubtreeHoldsWorker();
out << ",";
out << "nodeWorkload:" << getNodeWorkload();
out << ",";
out << "localWorkload:" << getLocalWorkload();
out << ",";
out << "totalWorkload:" << getTotalWorkload();
out << ",";
out << "maxWorkload:" << getMaxWorkload();
out << ",";
out << "minWorkload:" << getMinWorkload();
out << ",";
out << "cellIsAForkCandidate:" << getCellIsAForkCandidate();
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getEvenFlags(),
getAccessNumber(),
getResponsibleRank(),
getSubtreeHoldsWorker(),
getNodeWorkload(),
getLocalWorkload(),
getTotalWorkload(),
getMaxWorkload(),
getMinWorkload(),
getCellIsAForkCandidate(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 9;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._subtreeHoldsWorker))), 		&disp[8] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 17;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_SHORT,		 //accessNumber
MPI_INT,		 //responsibleRank
MPI_CHAR,		 //subtreeHoldsWorker
MPI_DOUBLE,		 //nodeWorkload
MPI_DOUBLE,		 //localWorkload
MPI_DOUBLE,		 //totalWorkload
MPI_DOUBLE,		 //maxWorkload
MPI_DOUBLE,		 //minWorkload
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //responsibleRank
1,		 //subtreeHoldsWorker
1,		 //nodeWorkload
1,		 //localWorkload
1,		 //totalWorkload
1,		 //maxWorkload
1,		 //minWorkload
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._responsibleRank))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._subtreeHoldsWorker))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._nodeWorkload))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._localWorkload))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._totalWorkload))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._maxWorkload))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._minWorkload))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[12] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[13] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[14] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[15] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[16] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#elif !defined(Parallel) && defined(SharedMemoryParallelisation) && defined(Debug)
particles::pit::records::Cell::PersistentRecords::PersistentRecords() {

}


particles::pit::records::Cell::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_isInside(isInside),
_state(state),
_level(level),
_evenFlags(evenFlags),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {

}

particles::pit::records::Cell::Cell() {

}


particles::pit::records::Cell::Cell(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords._isInside, persistentRecords._state, persistentRecords._level, persistentRecords._evenFlags, persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {

}


particles::pit::records::Cell::Cell(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {

}

particles::pit::records::Cell::~Cell() { }

std::string particles::pit::records::Cell::toString(const State& param) {
switch (param) {
case Leaf: return "Leaf";
case Refined: return "Refined";
case Root: return "Root";
}
return "undefined";
}

std::string particles::pit::records::Cell::getStateMapping() {
return "State(Leaf=0,Refined=1,Root=2)";
}


std::string particles::pit::records::Cell::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::Cell::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::Cell::PersistentRecords particles::pit::records::Cell::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::CellPacked particles::pit::records::Cell::convert() const{
return CellPacked(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::Cell::_log( "particles::pit::records::Cell" );

MPI_Datatype particles::pit::records::Cell::Datatype = 0;
MPI_Datatype particles::pit::records::Cell::FullDatatype = 0;


void particles::pit::records::Cell::initDatatype() {
{
Cell dummyCell[2];

const int Attributes = 5;
MPI_Datatype subtypes[Attributes] = {
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //isInside
1,		 //state
1,		 //level
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._isInside))), 		&disp[4] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::Datatype );
MPI_Type_commit( &Cell::Datatype );

}
{
Cell dummyCell[2];

const int Attributes = 13;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_CHAR,		 //isInside
MPI_INT,		 //state
MPI_INT,		 //level
MPI_INT,		 //evenFlags
MPI_SHORT,		 //accessNumber
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //isInside
1,		 //state
1,		 //level
DIMENSIONS,		 //evenFlags
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._isInside))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._state))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._level))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._evenFlags))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._accessNumber[0]))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[10] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[0]._numberOfParticlesInChildren))), 		&disp[11] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCell[1]._persistentRecords._cellIndex))), 		&disp[12] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &Cell::FullDatatype );
MPI_Type_commit( &Cell::FullDatatype );

}

}


void particles::pit::records::Cell::shutdownDatatype() {
MPI_Type_free( &Cell::Datatype );
MPI_Type_free( &Cell::FullDatatype );

}

void particles::pit::records::Cell::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "was not able to send message particles::pit::records::Cell "
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
msg << "testing for finished send task for particles::pit::records::Cell "
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
"particles::pit::records::Cell",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
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



void particles::pit::records::Cell::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "failed to start to receive particles::pit::records::Cell from node "
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
msg << "testing for finished receive task for particles::pit::records::Cell failed: "
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
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::Cell",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::Cell::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::Cell::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::PersistentRecords::PersistentRecords(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_cellIndex(cellIndex),
_meanVelocity(meanVelocity),
_meanCoordinate(meanCoordinate),
_myNorm(myNorm),
_level(level),
_accessNumber(accessNumber),
_numberOfLoadsFromInputStream(numberOfLoadsFromInputStream),
_numberOfStoresToOutputStream(numberOfStoresToOutputStream) {
setIsInside(isInside);
setState(state);
setEvenFlags(evenFlags);
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::CellPacked() {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._cellIndex, persistentRecords._meanVelocity, persistentRecords._meanCoordinate, persistentRecords._myNorm, persistentRecords.getIsInside(), persistentRecords.getState(), persistentRecords._level, persistentRecords.getEvenFlags(), persistentRecords._accessNumber, persistentRecords._numberOfLoadsFromInputStream, persistentRecords._numberOfStoresToOutputStream) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}


particles::pit::records::CellPacked::CellPacked(const int& cellIndex, const int& numberOfParticlesInChildren, const tarch::la::Vector<DIMENSIONS,double>& meanVelocity, const tarch::la::Vector<DIMENSIONS,double>& meanCoordinate, const tarch::la::Vector<DIMENSIONS,double>& myNorm, const bool& isInside, const State& state, const int& level, const std::bitset<DIMENSIONS>& evenFlags, const tarch::la::Vector<DIMENSIONS_TIMES_TWO,short int>& accessNumber, const int& numberOfLoadsFromInputStream, const int& numberOfStoresToOutputStream):
_persistentRecords(cellIndex, meanVelocity, meanCoordinate, myNorm, isInside, state, level, evenFlags, accessNumber, numberOfLoadsFromInputStream, numberOfStoresToOutputStream),_numberOfParticlesInChildren(numberOfParticlesInChildren) {
assertion((DIMENSIONS+3 < (8 * sizeof(short int))));

}

particles::pit::records::CellPacked::~CellPacked() { }

std::string particles::pit::records::CellPacked::toString(const State& param) {
return particles::pit::records::Cell::toString(param);
}

std::string particles::pit::records::CellPacked::getStateMapping() {
return particles::pit::records::Cell::getStateMapping();
}



std::string particles::pit::records::CellPacked::toString() const {
std::ostringstream stringstr;
toString(stringstr);
return stringstr.str();
}

void particles::pit::records::CellPacked::toString (std::ostream& out) const {
out << "("; 
out << "cellIndex:" << getCellIndex();
out << ",";
out << "numberOfParticlesInChildren:" << getNumberOfParticlesInChildren();
out << ",";
out << "meanVelocity:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanVelocity(i) << ",";
   }
   out << getMeanVelocity(DIMENSIONS-1) << "]";
out << ",";
out << "meanCoordinate:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMeanCoordinate(i) << ",";
   }
   out << getMeanCoordinate(DIMENSIONS-1) << "]";
out << ",";
out << "myNorm:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMyNorm(i) << ",";
   }
   out << getMyNorm(DIMENSIONS-1) << "]";
out << ",";
out << "isInside:" << getIsInside();
out << ",";
out << "state:" << toString(getState());
out << ",";
out << "level:" << getLevel();
out << ",";
out << "evenFlags:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getEvenFlags(i) << ",";
   }
   out << getEvenFlags(DIMENSIONS-1) << "]";
out << ",";
out << "accessNumber:[";
   for (int i = 0; i < DIMENSIONS_TIMES_TWO-1; i++) {
      out << getAccessNumber(i) << ",";
   }
   out << getAccessNumber(DIMENSIONS_TIMES_TWO-1) << "]";
out << ",";
out << "numberOfLoadsFromInputStream:" << getNumberOfLoadsFromInputStream();
out << ",";
out << "numberOfStoresToOutputStream:" << getNumberOfStoresToOutputStream();
out <<  ")";
}


particles::pit::records::CellPacked::PersistentRecords particles::pit::records::CellPacked::getPersistentRecords() const {
return _persistentRecords;
}

particles::pit::records::Cell particles::pit::records::CellPacked::convert() const{
return Cell(
getCellIndex(),
getNumberOfParticlesInChildren(),
getMeanVelocity(),
getMeanCoordinate(),
getMyNorm(),
getIsInside(),
getState(),
getLevel(),
getEvenFlags(),
getAccessNumber(),
getNumberOfLoadsFromInputStream(),
getNumberOfStoresToOutputStream()
);
}

#ifdef Parallel
tarch::logging::Log particles::pit::records::CellPacked::_log( "particles::pit::records::CellPacked" );

MPI_Datatype particles::pit::records::CellPacked::Datatype = 0;
MPI_Datatype particles::pit::records::CellPacked::FullDatatype = 0;


void particles::pit::records::CellPacked::initDatatype() {
{
CellPacked dummyCellPacked[2];

const int Attributes = 4;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //level
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //level
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._level))), 		&disp[3] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::Datatype );
MPI_Type_commit( &CellPacked::Datatype );

}
{
CellPacked dummyCellPacked[2];

const int Attributes = 11;
MPI_Datatype subtypes[Attributes] = {
MPI_INT,		 //cellIndex
MPI_DOUBLE,		 //meanVelocity
MPI_DOUBLE,		 //meanCoordinate
MPI_DOUBLE,		 //myNorm
MPI_INT,		 //level
MPI_SHORT,		 //accessNumber
MPI_INT,		 //numberOfLoadsFromInputStream
MPI_INT,		 //numberOfStoresToOutputStream
MPI_SHORT,		 //_packedRecords0
MPI_INT,		 //numberOfParticlesInChildren
MPI_UB		 // end/displacement flag
};

int blocklen[Attributes] = {
1,		 //cellIndex
DIMENSIONS,		 //meanVelocity
DIMENSIONS,		 //meanCoordinate
DIMENSIONS,		 //myNorm
1,		 //level
DIMENSIONS_TIMES_TWO,		 //accessNumber
1,		 //numberOfLoadsFromInputStream
1,		 //numberOfStoresToOutputStream
1,		 //_packedRecords0
1,		 //numberOfParticlesInChildren
1		 // end/displacement flag
};

MPI_Aint     disp[Attributes];

MPI_Aint base;
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]))), &base);
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._cellIndex))), 		&disp[0] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanVelocity[0]))), 		&disp[1] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._meanCoordinate[0]))), 		&disp[2] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._myNorm[0]))), 		&disp[3] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._level))), 		&disp[4] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._accessNumber[0]))), 		&disp[5] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfLoadsFromInputStream))), 		&disp[6] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._numberOfStoresToOutputStream))), 		&disp[7] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._persistentRecords._packedRecords0))), 		&disp[8] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[0]._numberOfParticlesInChildren))), 		&disp[9] );
MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellPacked[1]._persistentRecords._cellIndex))), 		&disp[10] );

for (int i=1; i<Attributes; i++) {
assertion1( disp[i] > disp[i-1], i );
}
for (int i=0; i<Attributes; i++) {
disp[i] -= base;
}
MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellPacked::FullDatatype );
MPI_Type_commit( &CellPacked::FullDatatype );

}

}


void particles::pit::records::CellPacked::shutdownDatatype() {
MPI_Type_free( &CellPacked::Datatype );
MPI_Type_free( &CellPacked::FullDatatype );

}

void particles::pit::records::CellPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
_senderDestinationRank = destination;

if (communicateBlocking) {

const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
if  (result!=MPI_SUCCESS) {
std::ostringstream msg;
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "was not able to send message particles::pit::records::CellPacked "
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
msg << "testing for finished send task for particles::pit::records::CellPacked "
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
"particles::pit::records::CellPacked",
"send(int)", destination,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
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



void particles::pit::records::CellPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, bool communicateBlocking) {
if (communicateBlocking) {

MPI_Status  status;
const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
_senderDestinationRank = status.MPI_SOURCE;
if ( result != MPI_SUCCESS ) {
std::ostringstream msg;
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "failed to start to receive particles::pit::records::CellPacked from node "
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
msg << "testing for finished receive task for particles::pit::records::CellPacked failed: "
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
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
triggeredTimeoutWarning = true;
}
if (
tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
(clock()>timeOutShutdown)
) {
tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
"particles::pit::records::CellPacked",
"receive(int)", source,tag,1
);
}
tarch::parallel::Node::getInstance().receiveDanglingMessages();
}

delete sendRequestHandle;

_senderDestinationRank = status.MPI_SOURCE;
#ifdef Debug
_log.debug("receive(int,int)", "received " + toString() ); 
#endif

}

}



bool particles::pit::records::CellPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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

int particles::pit::records::CellPacked::getSenderRank() const {
assertion( _senderDestinationRank!=-1 );
return _senderDestinationRank;

}
#endif




#endif


