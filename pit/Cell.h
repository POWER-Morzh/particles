// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _PARTICLES_PIT_CELL_H_ 
#define _PARTICLES_PIT_CELL_H_


#include "particles/pit/records/Cell.h"
#include "peano/grid/Cell.h"
#include "peano/heap/Heap.h"
#include "particles/pit/records/Particle.h"
#include "particles/pit/records/ParticleCompressed.h"


namespace particles { 
  namespace pit { 
    class Cell;

    /**
     * The heap we should use has critical impact on the runtime, in particular
     * its (algorithmic) latency. Therefore, using a runtime-encoded version
     * compressing zeros is not a good idea - also zero messages have to be
     * passed to mpi immediately to get rid at least of algorithmic latency.
     *
     * While the arguing above in principle is correct, it lacks one imporant
     * fact: The heap uses three different buffers (fork/join, boundary,
     * master/worker). And the choice of the predefined heap just is a choice
     * of different boundary buffers. Those however are not used to transport
     * particles at all in pit.
     */
    typedef peano::heap::PlainHeap<particles::pit::records::Particle>  ParticleHeap;
    typedef peano::heap::PlainHeap<particles::pit::records::ParticleCompressedPacked>  ParticleCompressedHeap;
  }
}


/**
 * PiT Cell
 * 
 * A cell in the PiT code is pretty simple: It basically only holds an index
 * for the heap. This index tells the heap which particles are assigned to the
 * cell.
 *
 * Besides that simple attribute, there is also a counter of particles. For
 * the leaves, this counter (set by the lift mechanism) equals the number of
 * particles assigned to the cell. For refined cells, the counter holds the
 * reduced values of all the children.
 *
 * @author Tobias Weinzierl
 */
class particles::pit::Cell: public peano::grid::Cell< particles::pit::records::Cell > { 
  private: 
    typedef class peano::grid::Cell< particles::pit::records::Cell >  Base;

  public:
    /**
     * Default constructor
     *
     * If we are working in assert mode, we set the particle pointer and the
     * accumulated data to -1. This is not required to make the algorithm work,
     * but it simplifies debugging to have well-defined attribute values here.
     */
    Cell();

    /**
     * Prescribed by framework.
     */
    Cell(const Base::DoNotCallStandardConstructor&);

    /**
     * Copy constructor
     */
    Cell(const Base::PersistentCell& argument);

    /**
     * Assign the cell an index
     */
    void init(int cellIndex);

    #ifdef Asserts
    bool isInitialised() const;
    #endif

    /**
     * @return Cell index for the heap. Has to be set before due to init().
     */
    int getCellIndex() const;

    /**
     * Clear number of particles in children
     *
     * This transition function is used to analyse the number of particles
     * bottom-up along the tree. It resets the number of children particles.
     * If the cell is refined, this reset basically passes 0. If the
     * operation is called for an unrefined cell, the argument typically
     * equals the number of particles contained in the leaf cell.
     */
    void clearNumberOfParticlesInChildren(int particles);

    /**
     * Reduction
     *
     * For the passed child cell, tell the current object to take the
     * particles contained there into account. This is the operation
     * realising the accumulative reduction.
     */
    void incNumberOfParticlesInChildren(const Cell& childCell);

    /**
     * Return number of particles
     *
     * For a leaf, this operation returns the particles contained within
     * the leaf. For a refined cell, it returns the particles contained
     * within all children. The latter is a recursive definition.
     *
     * @return Number of particles contained in the present cell volume.
     */
    int getNumberOfParticlesInChildren() const;

    /*
     * Get the mean velocity of Particles in this sell
     */
    tarch::la::Vector<DIMENSIONS, double> getMeanVelocity() const;

    /**
     * Set the mean velocity of Particles in this sell
     *
     * This operation sets the mean velocity only for leaf cells with particles.
     */
    void setMeanVelocity(const tarch::la::Vector<DIMENSIONS, double>& meanVelocity);

    /*
     * Get the mean Coordinate of Particles in this sell
     */
    tarch::la::Vector<DIMENSIONS, double> getMeanCoordinate() const;

    /**
     * Set the mean Coordinate of Particles in this sell
     *
     * This operation sets the mean velocity only for leaf cells with particles.
     */
    void setMeanCoordinate(const tarch::la::Vector<DIMENSIONS, double>& meanVelocity);

    /*
     * Get the norm which we saved in this sell
     */
    tarch::la::Vector<DIMENSIONS, double> getMyNorm() const;

    /**
     * Set some norm(which we want) in this sell
     *
     * This operation sets the mean velocity only for leaf cells with particles.
     */
    void setMyNorm(const tarch::la::Vector<DIMENSIONS, double>& myNorm);
};


#endif
