#include "particles/pit/Cell.h"


particles::pit::Cell::Cell():
  Base() { 
  #ifdef Asserts
  _cellData.setCellIndex(-1);
  _cellData.setNumberOfParticlesInChildren(0);
  #endif
}


particles::pit::Cell::Cell(const Base::DoNotCallStandardConstructor& value):
  Base(value) { 
  // Please do not insert anything here
}

particles::pit::Cell::Cell(const Base::PersistentCell& argument):
  Base(argument) {
}


void particles::pit::Cell::init(int cellIndex) {
  assertion( cellIndex>=0 );
  _cellData.setCellIndex(cellIndex);
}


int particles::pit::Cell::getCellIndex() const {
  assertion1( _cellData.getCellIndex()!=-1,toString() );
  return _cellData.getCellIndex();
}


void particles::pit::Cell::clearNumberOfParticlesInChildren(int particles) {
  assertion2(particles>=0, particles, toString());
  _cellData.setNumberOfParticlesInChildren(particles);
}


void particles::pit::Cell::incNumberOfParticlesInChildren(const Cell& childCell) {
  assertion2(childCell._cellData.getNumberOfParticlesInChildren()>=0, childCell.toString(), toString());
  assertion2(_cellData.getNumberOfParticlesInChildren()>=0, childCell.toString(), toString());
  _cellData.setNumberOfParticlesInChildren(
    _cellData.getNumberOfParticlesInChildren() +
    childCell._cellData.getNumberOfParticlesInChildren()
  );
}


int particles::pit::Cell::getNumberOfParticlesInChildren() const {
  assertion1(_cellData.getNumberOfParticlesInChildren()>=0, toString());
  return _cellData.getNumberOfParticlesInChildren();
}


#ifdef Asserts
bool particles::pit::Cell::isInitialised() const {
  return _cellData.getCellIndex() >= 0;
}
#endif


tarch::la::Vector<DIMENSIONS, double> particles::pit::Cell::getMeanVelocity() const {
  return _cellData.getMeanVelocity();
}


void particles::pit::Cell::setMeanVelocity(const tarch::la::Vector<DIMENSIONS, double>& meanVelocity) {
	_cellData.setMeanVelocity(meanVelocity);
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::Cell::getMeanCoordinate() const {
  return _cellData.getMeanCoordinate();
}


void particles::pit::Cell::setMeanCoordinate(const tarch::la::Vector<DIMENSIONS, double>& meanVelocity) {
	_cellData.setMeanCoordinate(meanVelocity);
}


tarch::la::Vector<DIMENSIONS, double> particles::pit::Cell::getMyNorm() const {
  return _cellData.getMyNorm();
}


void particles::pit::Cell::setMyNorm(const tarch::la::Vector<DIMENSIONS, double>& myNorm) {
        _cellData.setMyNorm(myNorm);
}
