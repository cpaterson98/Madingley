#include "Activity.h"
#include "Environment.h"
#include "MadingleyInitialisation.h"
#include "Maths.h"

Activity::Activity( ) {
    mWarmingTolerance = 0;
    mThermalSafetyMargin = 0;
    mTopt = 0;
    mCTmax = 0;
    mCTmin = 0;
    mAmbientTemp = 0;
    mDTR = 0;
    // Initialise ecological parameters for predation
    InitialiseActivityParameters( );
}

void Activity::InitialiseActivityParameters( ) {
    // Source: Deutsch et al (2008), Impacts of climate warming on terrestrial ecototherms across latitude, PNAS.
    mTerrestrialWarmingToleranceIntercept = 6.61;
    mTerrestrialWarmingToleranceSlope = 1.6;
    mTerrestrialTSMIntercept = 1.51;
    mTerrestrialTSMSlope = 1.53;

}

void Activity::AssignProportionTimeActive( GridCell& gcl, Cohort* actingCohort, unsigned currentTimestep, unsigned currentMonth, MadingleyInitialisation& params ) {
    //Only work on heterotroph cohorts
    if( params.mCohortFunctionalGroupDefinitions.GetTraitNames( "Heterotroph/Autotroph", actingCohort->mFunctionalGroupIndex ) == "heterotroph" ) {
        //Check if this is an endotherm or ectotherm
        bool Endotherm = params.mCohortFunctionalGroupDefinitions.GetTraitNames( "Endo/Ectotherm", actingCohort->mFunctionalGroupIndex ) == "endotherm";
        if( Endotherm ) {
            //Assumes the whole timestep is suitable for endotherms to be active - actual time active is therefore the proportion specified for this functional group.
            actingCohort->mProportionTimeActive = params.mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "proportion suitable time active", actingCohort->mFunctionalGroupIndex );
            if (actingCohort->mHibernating == true){
                actingCohort->mProportionTimeActive *= Constants::cWinterSlowDown; // Hibernators are not as active (but do wake up and feed occasionally
            }
        } else {
            //If ectotherm then use realm specific function
            if( !gcl.IsMarine( ) ) {
                actingCohort->mProportionTimeActive = CalculateProportionTimeSuitableTerrestrial( gcl, currentMonth, Endotherm ) *
                        params.mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "proportion suitable time active", actingCohort->mFunctionalGroupIndex );
            } else {
                //update from PLOS biology code which just had 1.0 here
                actingCohort->mProportionTimeActive = CalculateProportionTimeSuitableMarine( gcl, currentMonth, Endotherm ) *
                        params.mCohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup( "proportion suitable time active", actingCohort->mFunctionalGroupIndex );
            }
        }
    }
}

double Activity::CalculateProportionTimeSuitableTerrestrial( GridCell& gcl, unsigned currentMonth, bool endotherm ) {

    mAmbientTemp = Environment::Get( "Temperature", gcl );
    mDTR = Environment::Get( "DiurnalTemperatureRange", gcl );

    //Calculate the Warming tolerance and thermal safety margin given standard deviation of monthly temperature
    mWarmingTolerance = mTerrestrialWarmingToleranceSlope * Environment::Get( "SDTemperature", gcl ) + mTerrestrialWarmingToleranceIntercept;
    mThermalSafetyMargin = mTerrestrialTSMSlope * Environment::Get( "SDTemperature", gcl ) + mTerrestrialTSMIntercept;

    mTopt = mThermalSafetyMargin + Environment::Get( "AnnualTemperature", gcl );
    mCTmax = mWarmingTolerance + Environment::Get( "AnnualTemperature", gcl );


    double PerformanceStandardDeviation = ( mCTmax - mTopt ) / 12;

    mCTmin = mTopt - 4 * PerformanceStandardDeviation;

    return ProportionDaySuitable( );

}

double Activity::CalculateProportionTimeSuitableMarine( GridCell& gcl, unsigned currentMonth, bool endotherm ) {
    return 1.0;
}

double Activity::ProportionDaySuitable( ) {

    //Calculate the diurnal maximum in the current month
    double DTmax = mAmbientTemp + ( 0.5 * mDTR );
    double DTmin = mAmbientTemp - ( 0.5 * mDTR );

    double temp = 2 * ( mCTmax - mAmbientTemp ) / mDTR;
    if( mCTmax - DTmax > 0.0 ) {
        temp = 1.0;
    } else if( mCTmax - DTmin < 0.0 ) {
        temp = -1.0;
    }

    //Proportion of time for which ambient temperatures are greater than the critical upper temperature
    double POver = ( ( Maths::Get( )->Pi( ) / 2.0 ) - asin( temp ) ) / Maths::Get( )->Pi( );

    if( mCTmin - DTmax > 0.0 ) {
        temp = 1.0;
    } else if( mCTmin - DTmin < 0.0 ) {
        temp = -1.0;
    } else {
        temp = 2 * ( mCTmin - mAmbientTemp ) / mDTR;
    }
    //Proportion of time for which ambient temperatures are below the critical lower temperature
    double PBelow = 1 - ( ( Maths::Get( )->Pi( ) / 2.0 ) - asin( temp ) ) / Maths::Get( )->Pi( );

    return 1 - ( POver + PBelow );
}