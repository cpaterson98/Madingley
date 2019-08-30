#ifndef CONTINUOUSDISPERSAL_H
#define CONTINUOUSDISPERSAL_H

#include "Dispersal.h"
#include "UtilityFunctions.h"

class ContinuousDispersal : public Dispersal{
public:
    /** \brief Assigns all parameter values for continuous dispersal */
    ContinuousDispersal( );

    /** \brief Run dispersal
    @param gridForDispersal The model grid to run dispersal for 
    @param cohortToDisperse The cohort for which to apply the dispersal process 
    @param currentMonth The current model month */
    void Run( Grid&, Cohort*, const unsigned& ) override;

    /** \brief Calculate the average diffusive dispersal speed of individuals in a cohort given their body mass
    @param bodyMass The current body mass of an individual in the cohort 
    @return The (average) dispersal speed in kilometres per month */
    double CalculateDispersalSpeed( double );
    
    private:
    
    /** \brief Get new cell 
    @param madingleyGrid Model grid 
    @param newIndex 
    @param c Working cohort  
    @param newLocation   */
    //should be overridden - come back to this later
    void NewCell( Grid&, unsigned int, Cohort*, Location* );

    /** \brief The time units associated with this implementation of dispersal */
    std::string mTimeUnitImplementation;
    /** \brief Scalar relating dispersal speed to individual body mass */
    double mDispersalSpeedBodyMassScalar;
    /** \brief Body-mass exponent of the relationship between disperal speed and individual body mass */
    double mDispersalSpeedBodyMassExponent;
    /** \brief Scalar to convert from the time step units used by this formulation of dispersal to global model time step units */
    double mDeltaT;

};

#endif /* CONTINUOUSDISPERSAL_H */

