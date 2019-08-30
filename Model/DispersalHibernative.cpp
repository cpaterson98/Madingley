#include "DispersalHibernative.h"
#include "UtilityFunctions.h"
#include "Constants.h"

 DispersalHibernative::DispersalHibernative() {
    mTimeUnitImplementation = "month";
    mDispersalSpeedBodyMassScalar = 0.0278;
    mDispersalSpeedBodyMassExponent = 0.48;

    // Calculate the scalar to convert from the time step units used by this implementation of dispersal to the global model time step units
    mDeltaT = mUtilities.ConvertTimeUnits( Parameters::Get( )->GetTimeStepUnits( ), mTimeUnitImplementation );
   

    // Set the seed for the random number generator
    mRandomNumberA.SetSeedFromSystemTime();
  
}

void DispersalHibernative::Run ( Grid& gridForDispersal, Cohort* cohortToDisperse, const unsigned& currentMonth ) {
   //Calculate movement speed - correct/most accurate speeds need to be figured out later
    double adultMass = cohortToDisperse->mAdultMass;
    double dispersalSpeed = this->CalculateDispersalSpeed(adultMass); //* mDeltaT; - currently no change with time step
    
    
    
    //Moves towards preferred direction with 50% chance, otherwise choose randomly
    double bearing = 0;
   
        bearing = mRandomNumberA.GetUniform( ) * 2 * acos( -1. );
   
    
    //Creating new location
    std::pair<double, double> coords = mUtilities.NewLocationGivenDistance(cohortToDisperse->mCurrentLocation, dispersalSpeed, bearing);
    unsigned newIndex = mUtilities.FindGridCellGivenLongandLat(coords.first, coords.second);
    
    
    //Dispersing cohort
    if (newIndex == cohortToDisperse->mCurrentCell->GetIndex()){
        //If cohort hasn't changed cells, update the new location
        Location* newLocation = new Location ( ); 
    newLocation->SetCoordinates(coords.first, coords.second);
    newLocation->SetIndices(gridForDispersal.GetACell(newIndex).GetLatitudeIndex(), gridForDispersal.GetACell(newIndex).GetLongitudeIndex());
        
       
        // this is changing the cell this->NewCell(gridForDispersal, newIndex, cohortToDisperse, newLocation);
        cohortToDisperse->mCurrentLocation = (*newLocation);
    }
}   // If the change would leave the cell, ignore it as hibernating animals cannot venture outside the grid



// units? m/s?
double DispersalHibernative::CalculateDispersalSpeed( double bodyMass ) {
    double speed =  mDispersalSpeedBodyMassScalar * pow( bodyMass, mDispersalSpeedBodyMassExponent )*Constants::cHibernationSlowdown;
    return speed;
}
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

