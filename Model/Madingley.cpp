
#include "Madingley.h"
#include "Environment.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <list>

//# new
#include <dirent.h>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>
#include "GridCell.h"
#include <chrono>

//# end new

Madingley::Madingley( ) {
    // Set up list of global diagnostics
    SetUpGlobalDiagnosticsList( );

   
    // Initialise the cohort ID to zero
    mNextCohortID = 0;
    bool UseMadingleySpinUp = Parameters::Get( )->GetApplyModelSpinup( ); // read yes(1)/no(0) from csv 
    mParams = MadingleyInitialisation( mNextCohortID, mGlobalDiagnosticVariables["NumberOfCohortsInModel"],
      mGlobalDiagnosticVariables["NumberOfStocksInModel"], mModelGrid, UseMadingleySpinUp );


    mDispersalSet = new DispersalSet( );
    
    mStockLeafStrategy = mParams.mStockFunctionalGroupDefinitions.mTraitLookupFromIndex[ "leaf strategy" ];
    mCohortNutritionSource = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "nutrition source" ];
    mCohortThermoregulation = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "endo/ectotherm" ];
    mCohortReproductiveStrategy = mParams.mCohortFunctionalGroupDefinitions.mTraitLookupFromIndex[ "reproductive strategy" ];
   
    //############## new output for csv files
   // string NameOfOutputFile;
    //std::cin >> NameOfOutputFile;
    mOutputDirectory = "C://Users/science.intern2/Desktop/output/";
            //+ NameOfOutputFile + "/";
    std::cout << "Cohort outputs: " << mOutputDirectory << std::endl;
    //############## end new

}

void Madingley::Run( ) {
    auto start = std::chrono::system_clock::now();
    //# Console colors
    std::string default_console = "\033[0m"; //# default
    std::string clr1 = "\033[0;31m"; //# red
    std::string clr2 = "\033[0;35m"; //# magenta
    std::string clr3 = "\033[0;33m"; //# yellow
    std::string clr4 = "\033[0;32m"; //# green

    std::cout << clr1 << std::endl;
    std::cout << "Extinction scenario..." << std::endl;
    std::cout << "TimeStepStartExtinction: " << Parameters::Get( )->GetTimeStepStartExtinction( ) << std::endl;
    std::cout << "StartBodyMass: "  << Parameters::Get( )->GetStartBodyMass( ) << std::endl;
    std::cout << "EndBodyMass: " << Parameters::Get( )->GetEndBodyMass( ) << std::endl;
    std::cout << "StepBodyMass: " << Parameters::Get( )->GetStepBodyMass( ) << std::endl;
    std::cout << "SelectCarnivores: "  << Parameters::Get( )->GetSelectCarnivores( ) << std::endl;
    std::cout << "SelectOmnivores: " << Parameters::Get( )->GetSelectOmnivores( ) << std::endl;
    std::cout << "SelectHerbivores: " << Parameters::Get( )->GetSelectHerbivores( ) << std::endl;

    // Write out model run details to the console
    std::cout << clr4 << std::endl;
    std::cout << "Starting model run" << std::endl;
    unsigned RunParallel = 0;

    if(Parameters::Get( )->GetRunParallel( )==1){
        std::cout << "Running in parallel, using " << Parameters::Get( )->GetThreadNumber( ) << " Threads" << std::endl;
        RunParallel = 1;
    }else{
        std::cout << "Running in serial" << std::endl;
        RunParallel = 0;
    }
    std::cout<<default_console<<std::endl;

    // Store EcoTime
    double EcoTime = 0;
    mDispersals = 0;
    int SimulationInMonths_print = Parameters::Get( )->GetLengthOfSimulationInMonths( );

    //############## new create vector with relevant grid indices (terrestrial only)
    std::vector<int> TerrestrialGridcellIndices; 
    unsigned cellCounter = 0;
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gcl ) {
        if ( !gcl.IsMarine( ) ) {
            TerrestrialGridcellIndices.push_back( gcl.GetIndex( ) );
            cellCounter++;
           
        } 
    } );

    std::cout << "Running terrestrial only ( "<< cellCounter << " of the " << 
        Parameters::Get( )->GetNumberOfGridCells( ) << " gridcells )"<< std::endl;
    std::cout<< " " <<std::endl;

    // write initialization grid properties
    OutputGridCSV( 0 );

    // write initialization cohort properties
    OutputCSV( 9999999 );
    
    /// Run the model
    for( unsigned timeStep = 0; timeStep < Parameters::Get( )->GetLengthOfSimulationInMonths( ); timeStep += 1 ) {

        TimeStep::Get( )->SetMonthly( timeStep );

        std::cout << clr2<< "Running time step " << timeStep + 1 << " / " << SimulationInMonths_print << default_console << std::endl;

        // Get current time step and month
        mCurrentTimeStep = timeStep;
        
        
        mCurrentMonth = mUtilities.GetCurrentMonth( timeStep );
       
        mEcologyTimer.Start( );

        Environment::Update( mCurrentMonth );


        // Run in gridcell ecology in parallel or serial depending on configuration
        if(RunParallel==1){
            RunWithinCellsInParallel( cellCounter, TerrestrialGridcellIndices );
        }else{
            RunWithinCells( cellCounter, TerrestrialGridcellIndices );
        }

        mEcologyTimer.Stop( );
        std::cout << "Within grid ecology took: " << mEcologyTimer.GetElapsedTimeSecs( ) << std::endl;
        
        EcoTime += mEcologyTimer.GetElapsedTimeSecs( );

        mDispersalTimer.Start( );
        
        RunCrossGridCellEcology( mDispersals );
        
        mDispersalTimer.Stop( );
        std::cout << "Across grid ecology took: " << mDispersalTimer.GetElapsedTimeSecs( ) << std::endl;

        

        mOutputTimer.Start( );

        // Netcdf basic and grid outputs
        Output( timeStep );

        // Output cohort specifics
        if( timeStep  % 50 == 0 ) {
            std::cout << clr4 << "Writing cohort specifics..." << default_console << std::endl;
            OutputCSV( timeStep );
        }

        mOutputTimer.Stop( );
        std::cout << "Global Outputs took: " << mOutputTimer.GetElapsedTimeSecs( ) << std::endl;

        // Write the results of dispersal to the console
        std::cout << "Total Cohorts remaining " << mGlobalDiagnosticVariables["NumberOfCohortsInModel"] << std::endl;
        std::cout << default_console << std::endl;

        // Write model outputs which can be utilized for initialisation of new model run (use current model run as spin-up)
        if( timeStep == Parameters::Get( )->GetLengthOfSimulationInMonths( ) - 1 ) {
            CohortSpinUpOutput( );
            StockSpinUpOutput( );
        }
         
         

    }
    std::cout << "Total EcoTime: " << EcoTime << std::endl;
    std::cout << "Mean EcoTime: " << EcoTime/SimulationInMonths_print << std::endl;
}

void Madingley::RunWithinCells( unsigned cellCounter, std::vector<int> TerrestrialGridcellIndices ) {
    // Instantiate a class to hold thread locked global diagnostic variables
    ThreadVariables singleThreadDiagnostics( 0, 0, 0, mNextCohortID );
    cellCounter = 0;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gcl ) {

        RunWithinCellStockEcology( gcl );
        RunWithinCellCohortEcology( gcl, singleThreadDiagnostics );
        cellCounter++;
        
    } );
    // Update the variable tracking cohort unique IDs
    mNextCohortID = singleThreadDiagnostics.mNextCohortID;

    //std::cout << mNextCohortID << std::endl;

    // Take the results from the thread local variables and apply to the global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = singleThreadDiagnostics.mExtinctions - singleThreadDiagnostics.mCombinations;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = singleThreadDiagnostics.mProductions;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = mGlobalDiagnosticVariables["NumberOfCohortsInModel"] + singleThreadDiagnostics.mProductions - singleThreadDiagnostics.mExtinctions;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = singleThreadDiagnostics.mCombinations;
}

void Madingley::RunWithinCellsInParallel( unsigned cellCounter, std::vector<int> TerrestrialGridcellIndices ) {
    // Instantiate a class to hold thread locked global diagnostic variables


    #ifdef _OPENMP
    std::cout<<"Running RunWithinCellsInParallel..."<<endl;
    double startTimeTest = omp_get_wtime( );
    #endif

    list<ThreadVariables> partialsDiagnostics;
    unsigned gridCellIndexS;

    #pragma omp parallel num_threads(Parameters::Get( )->GetThreadNumber( )) shared(partialsDiagnostics)
    {
        ThreadVariables singleThreadDiagnostics( 0, 0, 0, mNextCohortID );

        #pragma omp for schedule(dynamic)
        for( unsigned gridCellIndex = 0; gridCellIndex < Parameters::Get( )->GetNumberOfGridCells( ); gridCellIndex++ )
        {   
            RunWithinCellStockEcology( mModelGrid.GetACell( gridCellIndex) );
            RunWithinCellCohortEcology( mModelGrid.GetACell( gridCellIndex ), singleThreadDiagnostics ); 
        }
        partialsDiagnostics.push_back(singleThreadDiagnostics);
    }//end parallel

    ThreadVariables globalDiagnostics( 0, 0, 0, mNextCohortID);
    for (list<ThreadVariables>::iterator it=partialsDiagnostics.begin(); it != partialsDiagnostics.end(); it++)
    {
        ThreadVariables tmp=*it;
        globalDiagnostics.mProductions+=tmp.mProductions;
        globalDiagnostics.mExtinctions+=tmp.mExtinctions;
        globalDiagnostics.mCombinations+=tmp.mCombinations;
    }

    // Update the variable tracking cohort unique IDs
    mNextCohortID = globalDiagnostics.mNextCohortID;

    // Take the results from the thread local variables and apply to the global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = globalDiagnostics.mExtinctions - globalDiagnostics.mCombinations;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = globalDiagnostics.mProductions;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = mGlobalDiagnosticVariables["NumberOfCohortsInModel"] + globalDiagnostics.mProductions - globalDiagnostics.mExtinctions;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = globalDiagnostics.mCombinations;

    #ifdef _OPENMP
    double endTimeTest = omp_get_wtime( );
    //std::cout << "RunWithinCellsInParallel( ) took: " << endTimeTest - startTimeTest << endl;
    #endif

}

void Madingley::RunWithinCellStockEcology( GridCell& gcl ) {

    if ( !gcl.IsMarine( ) ) {
    // Create a local instance of the stock ecology class
    EcologyStock MadingleyEcologyStock;
    // Get the list of functional group indices for autotroph stocks
    std::vector<int> AutotrophStockFunctionalGroups = mParams.mStockFunctionalGroupDefinitions.GetFunctionalGroupIndex( "Heterotroph/Autotroph", "Autotroph", false );
    // Loop over autotroph functional groups
    for( unsigned FunctionalGroup: AutotrophStockFunctionalGroups ) {
        for( auto& ActingStock: gcl.mStocks[FunctionalGroup] ) {

            // Run stock ecology
            MadingleyEcologyStock.RunWithinCellEcology( gcl, ActingStock, mCurrentTimeStep, mCurrentMonth, mParams );
        }
    }
    }

}

void Madingley::RunWithinCellCohortEcology( GridCell& gcl, ThreadVariables& partial ) {
    // Local instances of classes
    // Initialize ecology for stocks and cohorts - needed fresh every timestep?
    if ( !gcl.IsMarine( ) ) { 
       
             EcologyCohort mEcologyCohort;
    
    mEcologyCohort.InitialiseEating( gcl, mParams );
    
    Activity CohortActivity;
    //clock_t begin = clock();
    std::vector< std::vector<int> > SortedCohortIndices;
    
    SortedCohortIndices = DetSortIndicesCohorts( gcl, false );
    
    //clock_t end = clock();
    //std::cout << double(end - begin) / CLOCKS_PER_SEC << std::endl;
    //std::cout << "##############" << std::endl;

     //std::vector< std::vector<int> > SortedCohortIndices;
   // SortedCohortIndices = mGridCell.GetSortedCohortIndices();
    // for(unsigned i = 0; i < SortedCohortIndices.size(); i++) std::cout << SortedCohortIndices[i].size() << std::endl;
//[&]( Cohort* c )
    // Loop over randomly ordered gridCellCohorts to implement biological functions
     
        // Perform all biological functions except dispersal (which is cross grid cell)
           
       
        std::vector<unsigned > RandomCohortOrder;
        std::vector<std::pair<int, int> > indexedList;
        unsigned TotalCohorts = 0;
        int totalcarnivores = 0;
        for( int functionalTypeIndex = 0; functionalTypeIndex < gcl.mCohorts.size( ); functionalTypeIndex++ ) {
            // Work through the list of cohorts and create a list of pairs so as to be able to lookup cohorts easily
            //for( int cohortNum = 0; cohortNum < gcl.mCohorts.size( ); cohortNum++ ) {
                for( int cohortNum = 0; cohortNum < gcl.mCohorts[ functionalTypeIndex ].size( ); cohortNum++ ) {
                indexedList.push_back( std::make_pair( functionalTypeIndex, cohortNum ) );
                TotalCohorts++;
            }
        }
        
               //despite the name of this utility, it actually returns a random list, but with a given fixed seed
        RandomCohortOrder = mUtilities.NonRandomlyOrderedCohorts( TotalCohorts, mCurrentTimeStep );

        for( int i = 0; i < RandomCohortOrder.size( ); i++ ) {
          
        
            Cohort* c = gcl.mCohorts[indexedList[RandomCohortOrder[i]].first][indexedList[RandomCohortOrder[i]].second];
            if (c->mFunctionalGroupIndex == 1 || c->mFunctionalGroupIndex == 4){
             totalcarnivores +=1;
           }
            
         //Start hibernation if necessary
            if (Constants::cHibernationEnabled == true){
            if(c->mCurrentLocation.mLatitude > Constants::cminLatforHib) {
                c->mGoodPlaceToHibernate = true;
            }else{c->mGoodPlaceToHibernate = false;}

        bool timeisgood;
        bool tempisgood;
       
           if (mCurrentMonth == 0 % 12 || mCurrentMonth == 1 % 12 || mCurrentMonth == 11 % 12 ){
                 timeisgood = true;
           }else{ timeisgood = false;}
           
        if (Environment::Get( "Temperature", c->GetCurrentCell( ) ) < 2){
            tempisgood = true;}
        else{tempisgood = false;}
        
        if ( timeisgood == true && tempisgood == true) {
            c->mGoodTimeToHibernate = true;}
           else{c->mGoodTimeToHibernate = false;}
            
            if (c->mWinterSurvivalTechnique ==1 && c->mGoodPlaceToHibernate == true && c->mGoodTimeToHibernate == true){
                c->mHibernating = true;
               //std::cout << "In this cell, this cohort is Hibernating." << endl;
           } else{
              c->mHibernating = false;
               // std::cout << "In this cell, this cohort is very much awake." << endl;
            }
       
                                      }
         // end of hibernation code: 
    
        
               if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) != 0 && c->mCohortAbundance > Parameters::Get( )->GetExtinctionThreshold( ) ) {
            
            CohortActivity.AssignProportionTimeActive( gcl, c, mCurrentTimeStep, mCurrentMonth, mParams );

            // Run ecology
            mEcologyCohort.RunWithinCellEcology( gcl, c, mCurrentTimeStep, partial, mCurrentMonth, mParams, SortedCohortIndices);
            // Update the properties of the acting cohort
           
            mEcologyCohort.UpdateEcology( gcl, c, mCurrentTimeStep );
            Cohort::ResetMassFluxes( );
            // Check that the mass of individuals in this cohort is still >= 0 after running ecology
           assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );
        }
        // Check that the mass of individuals in this cohort is still >= 0 after running ecology
       if( gcl.mCohorts[c->mFunctionalGroupIndex].size( ) > 0 )
        { assert( c->mIndividualBodyMass >= 0.0 && "Biomass < 0 for this cohort" );}
        }
     
        
    //for(unsigned i = 0; i < SortedCohortIndices[10].size(); i++) std::cout << SortedCohortIndices[10][i] << "  ";
    //SortedCohortIndices.clear();
    //std::cout << "##" << std::endl;

        
        
        
         for( auto c: GridCell::mNewCohorts ) {
        gcl.InsertCohort( c ); //HELP - ERROR IS HERE SOMEWHERE!!!
       

        if( c->mDestinationCell != &gcl ) std::cout << "what? wrong cell?" << std::endl;
        
    }
    partial.mProductions += GridCell::mNewCohorts.size( );
    GridCell::mNewCohorts.clear( );
    RunExtinction( gcl, partial );
    // Merge cohorts, if necessary
    if( gcl.GetNumberOfCohorts( ) > Parameters::Get( )->GetMaximumNumberOfCohorts( ) ) {
        mCohortMerger.ResetRandom( );
        partial.mCombinations += mCohortMerger.MergeToReachThresholdFast( gcl, mParams );

        //Run extinction a second time to remove those cohorts that have been set to zero abundance when merging
        RunExtinction( gcl, partial );
    }
    //# Megafauna extinction
    int TimeStepStartExtinction = Parameters::Get( )->GetTimeStepStartExtinction( ); //# start (time step in # months)
    int ExtMassThresholdStart = Parameters::Get( )->GetStartBodyMass( ); //# (g)
    int ExtMassThresholdEnd = Parameters::Get( )->GetEndBodyMass( );  //# (g)
    int ExtMassThresholdStep = Parameters::Get( )->GetStepBodyMass( ); //# (g)
    int ExtAdultMass = 9999999; //#temp placeholder
    if( mCurrentTimeStep >= TimeStepStartExtinction ) {

        int ExtAdultMass = ExtMassThresholdStart - ( ( mCurrentTimeStep - TimeStepStartExtinction ) * ExtMassThresholdStep );
        if( ExtAdultMass <= ExtMassThresholdEnd ) { ExtAdultMass = ExtMassThresholdEnd; }

        mGlobalDiagnosticVariables["ExtAdultMass"] = ExtAdultMass;
        MegaFaunaExtinction_v1( gcl, partial, ExtAdultMass );
    }
    //# end megafauna extinction 
        

       // if( totalcarnivores > 0){
            //std::cout << totalcarnivores << std::endl;
      //  }
        }
   
  
}

void Madingley::RunExtinction( GridCell& gcl, ThreadVariables& partial ) {

    // Loop over cohorts and remove any whose abundance is below the extinction threshold
    std::vector<Cohort*> CohortsToRemove; 
    gcl.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
        if( c->mCohortAbundance - Parameters::Get( )->GetExtinctionThreshold( ) <= 0 || c->mIndividualBodyMass <= 0 ) {
            CohortsToRemove.push_back( c );
            partial.mExtinctions += 1;
        }
    } );

    // Code to add the biomass to the biomass pool and dispose of the cohort
    for( auto c: CohortsToRemove ) {

        // Add biomass of the extinct cohort to the organic matter pool
        double deadMatter = ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance;
        if( deadMatter < 0 ) std::cout << "Dead " << deadMatter << std::endl;
        Environment::Get( "Organic Pool", c->GetCurrentCell( ) ) += deadMatter;
        assert( Environment::Get( "Organic Pool", c->GetCurrentCell( ) ) >= 0 && "Organic pool < 0" );

        // Remove the extinct cohort from the list of cohorts
        gcl.RemoveCohort( c );
    }
    for( auto c: CohortsToRemove ) {delete(c);}
}

//# MegaFaunaExtinction: new function, incorporates megafauna extinctions
void Madingley::MegaFaunaExtinction_v1( GridCell& gcl, ThreadVariables&, int ExtAdultMass ) {

    std::string color1 = "\033[0;31m"; //# red
    std::string color2 = "\033[0;35m"; //# magenta
    std::string default_console = "\033[0m"; //# default color

    // Loop over cohorts and remove
    std::vector<Cohort*> MegaFaunaCohortsToRemove;
    gcl.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
        if( c->mAdultMass > ExtAdultMass &&
           ( (mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", c->mFunctionalGroupIndex) == "carnivore"
              && Parameters::Get( )->GetSelectCarnivores( ) == 1 ) ||
            (mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", c->mFunctionalGroupIndex) == "omnivore"
             && Parameters::Get( )->GetSelectOmnivores( ) == 1 ) ||
            (mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", c->mFunctionalGroupIndex) == "herbivore"
             && Parameters::Get( )->GetSelectHerbivores( ) == 1 ) )
           ) {
            //###### check properties of cohorts going extinct
            if( true ) {
                std::cout << color2 << "Cohort extinct: " <<
                mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Nutrition source", c->mFunctionalGroupIndex ) << ", " <<
                mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Realm", c->mFunctionalGroupIndex ) << ", " <<
                mParams.mCohortFunctionalGroupDefinitions.GetTraitNames( "Endo/Ectotherm", c->mFunctionalGroupIndex ) <<
                ", c->mIndividualBodyMass = " << c->mIndividualBodyMass << ", c->mAdultMass = " << c->mAdultMass << ", c->mJuvenileMass = " << c->mJuvenileMass <<
                ", GridCell = " << &gcl << ", lat = " << gcl.GetLatitudeIndex( ) << ", long = " << gcl.GetLongitudeIndex( ) << default_console << std::endl;
            }
            //##### end check properties of cohorts going extinct
            MegaFaunaCohortsToRemove.push_back( c );
        }
    } );

    // Remove cohorts
    for( auto c: MegaFaunaCohortsToRemove ) {

        // Remove the extinct cohort from the list of cohorts
        gcl.RemoveCohort( c );
    }
    for( auto c: MegaFaunaCohortsToRemove ) {delete(c);}

    mGlobalDiagnosticVariables["NumberOfCohortForcedToExtinction"] += MegaFaunaCohortsToRemove.size( );

}//# end MegaFaunaExtinction

void Madingley::RunCrossGridCellEcology( unsigned& dispersals ) {
    // Loop through each grid cell, and run dispersal for each.
    // In the original model a new dispersal object is made every timestep - this resets the random number generators
    mDispersalSet->ResetRandoms( );
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & c ) {
        mDispersalSet->RunCrossGridCellEcologicalProcess( c, mModelGrid, mParams, mCurrentMonth );
    } );

    // Apply the changes from dispersal
    mDispersalSet->UpdateCrossGridCellEcology( dispersals );
    
}

void Madingley::SetUpGlobalDiagnosticsList( ) {
    // Add global diagnostic variables
    mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsProduced"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsCombined"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfCohortsInModel"] = 0.0;
    mGlobalDiagnosticVariables["NumberOfStocksInModel"] = 0.0;
    mGlobalDiagnosticVariables["ExtAdultMass"] = 999999.0; //#
    mGlobalDiagnosticVariables["NumberOfCohortForcedToExtinction"] = 0.0; //#
}


//##
std::vector< std::vector<int> > Madingley::DetSortIndicesCohorts( GridCell& gcl, bool print ) {

    // Vectors of mIndividualBodyMass, stored per function group
    // ("columns" are functional group index, "rows" are filled with the corresponding mIndividualBodyMass)
    std::vector< std::vector<double> > x(20,std::vector<double>(0));
    
  //  gcl.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
  //      x[c->mFunctionalGroupIndex].push_back(c->mIndividualBodyMass);
  //  } );

    for (auto& cell : gcl.mMooreNeighbourhood){
            for (vector<Cohort*> vc : cell->mCohorts){
                for (Cohort* c : vc){
                    x[c->mFunctionalGroupIndex].push_back(c->mIndividualBodyMass);
                }
            }
        }

    // Empty vector for storing sorted indices
    std::vector< std::vector<int> > y2(20,std::vector<int>(0));

    // Determine and store sorted indices per functional group
    for( int ii = 0; ii < 20; ii++ ) {
      if( x[ii].size()>0 ){
        std::vector<double> x2 = x[ii];

        std::vector<int> y(x2.size());

        std::size_t n(0);
        std::generate(y.begin(), y.end(), [&]{ return n++; });
         std::sort( y.begin(), y.end(), [&](int i1, int i2) { return x2[i1] < x2[i2]; } );

        for( unsigned yy = 0; yy < y.size(); yy++) y2[ii].push_back(y[yy]);
        y.clear();
      }
    }

    x.clear();

    // Print y2 to console
    if (print == true){
      for ( unsigned yy = 0; yy < y2.size(); yy++ ) {
        for ( unsigned tt = 0; tt < y2[yy].size(); tt++ ) std::cout << y2[yy][tt] << ' ';
        std::cout << "|||" <<std::endl;
      }
    }
    return y2;
}
//##

//##
std::string Madingley::ReturnOutputfolder( ) {
    std::string last_dir; const char* PATH = "C://Users/science.intern2/Desktop/output/";
    
    DIR *dir = opendir(PATH); struct dirent *entry = readdir(dir);

    while (entry != NULL)
   {
      //if (entry->d_type == DT_DIR)
         //  if (S_ISDIR(entry->st_mode))
       // {
            last_dir = entry->d_name;
       // }
       // entry = readdir(dir);
    }

    closedir(dir);

    return last_dir;
}

//##

//##
void Madingley::OutputCSV( unsigned step ) {
        stringstream ss;
    ss << step;
    string number = ss.str();
    std::string csvName = mOutputDirectory + number + "_month.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "c.mFunctionalGroupIndex" << "," << "c.mWinterSurvivalTechnique" << "," << "c.mAdultMass" << "," << "c.mIndividualBodyMass" << "," << "c.mCohortAbundance" << "," <<
        "c.mLogOptimalPreyBodySizeRatio" << "," << "step" << "," <<  "gridCell" << "," << "lat" << "," << "long" << std::endl;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                //std::cout << exp(c->mLogOptimalPreyBodySizeRatio) << std::endl;
                CSV << c->mFunctionalGroupIndex << "," << c->mWinterSurvivalTechnique << "," << c->mAdultMass << "," << c->mIndividualBodyMass << "," << c->mCohortAbundance << "," <<
                c->mLogOptimalPreyBodySizeRatio << "," << step << "," <<  gridCell.GetIndex( ) << "," <<
                gridCell.GetLatitudeIndex( ) << "," << gridCell.GetLongitudeIndex( ) << std::endl;

        } );
    } );
   CSV.close();
}
//##


//##
void Madingley::OutputGridCSV( unsigned step ) {
    stringstream ss;
    ss << step;
    string number = ss.str();
    std::string csvName = mOutputDirectory + number + "_grid.csv";
    ofstream CSV;
    CSV.open (csvName);
    CSV << "Index" << "," << "Lat_Index" << "," << "Long_Index" << "," << "Cell_Area" << "," <<
        "Cell_Width" << "," << "Cell_Height" << std::endl;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
                CSV << gridCell.GetIndex( ) << "," <<
                gridCell.GetLatitudeIndex( ) << "," << 
                gridCell.GetLongitudeIndex( ) << "," << 
                gridCell.GetCellArea( ) << "," << 
                gridCell.GetCellWidth( ) << "," << 
                gridCell.GetCellHeight( ) << std::endl;

    } );
   CSV.close();
}
//##

//##
void Madingley::CohortSpinUpOutput( ) {
    std::string csvName = mOutputDirectory + "CohortSpinUpOutput" + ".csv";
    ofstream CSV;
    CSV.open (csvName);
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                CSV <<
                gridCell.GetIndex( ) << "," <<
                c->mFunctionalGroupIndex << "," <<
                c->mJuvenileMass << "," <<
                c->mAdultMass << "," <<
                c->mIndividualBodyMass << "," <<
                c->mCohortAbundance << "," <<
                c->mLogOptimalPreyBodySizeRatio << "," <<
                0 << "," << // mBirthTimeStep
                c->mProportionTimeActive << "," <<
                c->mWinterSurvivalTechnique <<        
                std::endl;
        } );
    } );
   CSV.close();
}
//##

//##
void Madingley::StockSpinUpOutput( ) {
    std::string csvName2 = mOutputDirectory + "StockSpinUpOutput" + ".csv";
    ofstream CSV2;
    CSV2.open (csvName2);
    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        gridCell.ApplyFunctionToAllStocks( [&]( Stock & s ) {
              CSV2 <<
              gridCell.GetIndex( ) << "," <<
              s.mFunctionalGroupIndex << "," <<
              s.mTotalBiomass <<
              std::endl;
        } );
    } );
}
//##

void Madingley::Output( unsigned step ) {
    double totalLivingBiomass = 0;
    double totalBiomass = 0;

    double organicMatterPool = 0;
    double respiratoryPool = 0;

    double totalStockBiomass = 0;

    double totalCohortBiomass = 0;
    long totalCohorts = 0;
    double totalCohortAbundance = 0;

    mModelGrid.ApplyFunctionToAllCells( [&]( GridCell & gridCell ) {
        double organicMatterThisCell = Environment::Get( "Organic Pool", gridCell ) / 1000.;
        double respirationThisCell = Environment::Get( "Respiratory CO2 Pool", gridCell ) / 1000.;

                double cohortBiomassThisCell = 0;
                double stockBiomassThisCell = 0;

                double phytoplanktonBiomassThisCell = 0;
                double evergreenBiomassThisCell = 0;
                double deciduousBiomassThisCell = 0;

                double cohortAbundanceThisCell = 0;
                double herbivoreBiomassThisCell = 0;
                double herbivoreAbundanceThisCell = 0;
                double omnivoreBiomassThisCell = 0;
                double omnivoreAbundanceThisCell = 0;
                double carnivoreBiomassThisCell = 0;
                double carnivoreAbundanceThisCell = 0;

                double ectothermBiomassThisCell = 0;
                double ectothermAbundanceThisCell = 0;
                double endothermBiomassThisCell = 0;
                double endothermAbundanceThisCell = 0;

                double iteroparousBiomassThisCell = 0;
                double iteroparousAbundanceThisCell = 0;
                double semelparousBiomassThisCell = 0;
                double semelparousAbundanceThisCell = 0;

                //###### custom madingley outputs (turn on/off in OutputControlParameters.csv)

                // split carnivores (endo/ecto)
                double carnivoreEndothermBiomassThisCell = 0;
                double carnivoreEndothermAbundanceThisCell = 0;
                double carnivoreEctothermBiomassThisCell = 0;
                double carnivoreEctothermAbundanceThisCell = 0;

                // split omnivores (endo/ecto)
                double omnivoreEndothermBiomassThisCell = 0;
                double omnivoreEndothermAbundanceThisCell = 0;
                double omnivoreEctothermBiomassThisCell = 0;
                double omnivoreEctothermAbundanceThisCell = 0;

                // split herbivores (endo/ecto)
                double herbivoreEndothermBiomassThisCell = 0;
                double herbivoreEndothermAbundanceThisCell = 0;
                double herbivoreEctothermBiomassThisCell = 0;
                double herbivoreEctothermAbundanceThisCell = 0;

                //##### end new code


                //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                // Mega carnivores endotherm
                double MCarnivoreEndothermBMThisCell = 0;
                double MCarnivoreEndothermAThisCell = 0;
                // Mega carnivores ectotherm
                double MCarnivoreEctothermBMThisCell = 0;
                double MCarnivoreEctothermAThisCell = 0;

                // Mega Herbivores endotherm
                double MHerbivoreEndothermBMThisCell = 0;
                double MHerbivoreEndothermAThisCell = 0;
                // Mega Herbivores ectotherm
                double MHerbivoreEctothermBMThisCell = 0;
                double MHerbivoreEctothermAThisCell = 0;
                
                // Mega Omnivores endotherm
                double MOmnivoreEndothermBMThisCell = 0;
                double MOmnivoreEndothermAThisCell = 0;
                // Mega Omnivores ectotherm
                double MOmnivoreEctothermBMThisCell = 0;
                double MOmnivoreEctothermAThisCell = 0;

                ////////////////////////////////////////
                // Store largest Carnivores
                double MaxCarnivores = 0;
                // Store largest Omnivores
                double MaxOmnivores = 0;
                // Store largest Herbivores
                double MaxHerbivores = 0;
                /////////////////////////////////////

                //######  end additional megafauna specific cohorts statistics

                organicMatterPool += organicMatterThisCell;
                respiratoryPool += respirationThisCell;

                gridCell.ApplyFunctionToAllCohorts( [&]( Cohort* c ) {
                    totalCohorts += 1;
                    totalCohortAbundance += c->mCohortAbundance;

                            double cohortBiomass = ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                            cohortBiomassThisCell += cohortBiomass;
                            totalCohortBiomass += cohortBiomass;
                            cohortAbundanceThisCell += c->mCohortAbundance;

                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" ) {
                        herbivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" ) {
                        omnivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" ) {
                        carnivoreBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreAbundanceThisCell += c->mCohortAbundance;
                    }

                    if( mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        ectothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        ectothermAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        endothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        endothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    if( mCohortReproductiveStrategy[ c->mFunctionalGroupIndex ] == "iteroparity" ) {
                        iteroparousBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        iteroparousAbundanceThisCell += c->mCohortAbundance;
                    } else if( mCohortReproductiveStrategy[ c->mFunctionalGroupIndex ] == "semelparity" ) {
                        semelparousBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        semelparousAbundanceThisCell += c->mCohortAbundance;
                    }

                    //###### custom madingley outputs (turn on/off in OutputControlParameters.csv)
                    // carni + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        carnivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // carni + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        carnivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        carnivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // omni + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        omnivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // omni + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        omnivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        omnivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // herbi + endo
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" ) {
                        herbivoreEndothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreEndothermAbundanceThisCell += c->mCohortAbundance;
                    }

                    // herbi + ecto
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" ) {
                        herbivoreEctothermBiomassThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        herbivoreEctothermAbundanceThisCell += c->mCohortAbundance;
                    }
                    //###### end new code


                    //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                    // Mega carnivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 21000 ) {
                        MCarnivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MCarnivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega carnivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 21000) {
                        MCarnivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MCarnivoreEctothermAThisCell  += c->mCohortAbundance;
                    }

                    // Mega Herbivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MHerbivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MHerbivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega Herbivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MHerbivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MHerbivoreEctothermAThisCell  += c->mCohortAbundance;
                    }
                    
                    // Mega Omnivores endotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "endotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MOmnivoreEndothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MOmnivoreEndothermAThisCell  += c->mCohortAbundance;
                    }
                    // Mega Omnivores ectotherm
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" && 
                        mCohortThermoregulation[ c->mFunctionalGroupIndex ] == "ectotherm" && 
                        c->mIndividualBodyMass > 35000 ) {
                        MOmnivoreEctothermBMThisCell += ( c->mIndividualBodyMass + c->mIndividualReproductivePotentialMass ) * c->mCohortAbundance / 1000.;
                        MOmnivoreEctothermAThisCell  += c->mCohortAbundance;
                    }
                    
                    ////////////////////////////////////////
                    // Store largest Carnivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "carnivore" ) {
                        if( c->mAdultMass > MaxCarnivores ) MaxCarnivores = c->mAdultMass;
                    }
                    // Store largest Omnivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "omnivore" ) {
                        if( c->mAdultMass > MaxOmnivores ) MaxOmnivores = c->mAdultMass;
                    }
                    // Store largest Herbivores
                    if( mCohortNutritionSource[ c->mFunctionalGroupIndex ] == "herbivore" ) {
                        if( c->mAdultMass > MaxHerbivores ) MaxHerbivores = c->mAdultMass;
                    }
                    /////////////////////////////////////


                    //######  end additional megafauna specific cohorts statistics
                } );
                
        gridCell.ApplyFunctionToAllStocks( [&]( Stock & s ) {
            double thisStockBiomass = s.mTotalBiomass / 1000.; //convert from g to kg
            stockBiomassThisCell += thisStockBiomass;
            totalStockBiomass += thisStockBiomass;

            if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "na" ) phytoplanktonBiomassThisCell += thisStockBiomass;
            else if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "deciduous" ) deciduousBiomassThisCell += thisStockBiomass;
            else if( mStockLeafStrategy[ s.mFunctionalGroupIndex ] == "evergreen" ) evergreenBiomassThisCell += thisStockBiomass;
            } );

        double biomassThisCell = cohortBiomassThisCell + stockBiomassThisCell + respirationThisCell + organicMatterThisCell;

                //# outputs included in MonthlyGridOutputs.nc
                DataRecorder::Get( )->SetDataOn( "BiomassDensity", gridCell.GetIndex( ), biomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "AbundanceDensity", gridCell.GetIndex( ), cohortAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "AutotrophBiomassDensity", gridCell.GetIndex( ), stockBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HeterotrophBiomassDensity", gridCell.GetIndex( ), cohortBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "PhytoplanktonBiomassDensity", gridCell.GetIndex( ), phytoplanktonBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "DeciduousBiomassDensity", gridCell.GetIndex( ), deciduousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EvergreenBiomassDensity", gridCell.GetIndex( ), evergreenBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreBiomassDensity", gridCell.GetIndex( ), herbivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreAbundanceDensity", gridCell.GetIndex( ), herbivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreBiomassDensity", gridCell.GetIndex( ), omnivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreAbundanceDensity", gridCell.GetIndex( ), omnivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreBiomassDensity", gridCell.GetIndex( ), carnivoreBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreAbundanceDensity", gridCell.GetIndex( ), carnivoreAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EctothermBiomassDensity", gridCell.GetIndex( ), ectothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EctothermAbundanceDensity", gridCell.GetIndex( ), ectothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EndothermBiomassDensity", gridCell.GetIndex( ), endothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "EndothermAbundanceDensity", gridCell.GetIndex( ), endothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "IteroparityBiomassDensity", gridCell.GetIndex( ), iteroparousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "IteroparityAbundanceDensity", gridCell.GetIndex( ), iteroparousAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "SemelparityBiomassDensity", gridCell.GetIndex( ), semelparousBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "SemelparityAbundanceDensity", gridCell.GetIndex( ), semelparousAbundanceThisCell / gridCell.GetCellArea( ) );


                //###### new madingley grid outputs (turn on/off in OutputControlParameters.csv)
                DataRecorder::Get( )->SetDataOn( "CarnivoreEndothermBiomass", gridCell.GetIndex( ), carnivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEndothermAbundance", gridCell.GetIndex( ), carnivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEctothermBiomass", gridCell.GetIndex( ), carnivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "CarnivoreEctothermAbundance", gridCell.GetIndex( ), carnivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                DataRecorder::Get( )->SetDataOn( "OmnivoreEndothermBiomass", gridCell.GetIndex( ), omnivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEndothermAbundance", gridCell.GetIndex( ), omnivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEctothermBiomass", gridCell.GetIndex( ), omnivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "OmnivoreEctothermAbundance", gridCell.GetIndex( ), omnivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                DataRecorder::Get( )->SetDataOn( "HerbivoreEndothermBiomass", gridCell.GetIndex( ), herbivoreEndothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEndothermAbundance", gridCell.GetIndex( ), herbivoreEndothermAbundanceThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEctothermBiomass", gridCell.GetIndex( ), herbivoreEctothermBiomassThisCell / gridCell.GetCellArea( ) );
                DataRecorder::Get( )->SetDataOn( "HerbivoreEctothermAbundance", gridCell.GetIndex( ), herbivoreEctothermAbundanceThisCell / gridCell.GetCellArea( ) );

                //# output gridcell and environment properties
                DataRecorder::Get( )->SetDataOn( "GridCellIndex", gridCell.GetIndex( ), gridCell.GetIndex( ) ); // gridcell index
                DataRecorder::Get( )->SetDataOn( "GridCellArea", gridCell.GetIndex( ), gridCell.GetCellArea( ) ); // gridcell area km2
                DataRecorder::Get( )->SetDataOn( "GridCellWidth", gridCell.GetIndex( ), gridCell.GetCellWidth( ) ); // gridcell width km
                DataRecorder::Get( )->SetDataOn( "GridCellHeight", gridCell.GetIndex( ), gridCell.GetCellHeight( ) ); // gridcell height km

                DataRecorder::Get( )->SetDataOn( "GridCellTemperature", gridCell.GetIndex( ), Environment::Get( "Temperature", gridCell.GetIndex( ) ) ); // gridcell temperature
                DataRecorder::Get( )->SetDataOn( "GridCelluVel", gridCell.GetIndex( ), Environment::Get( "uVel", gridCell.GetIndex( ) ) ); // gridcell uVel
                DataRecorder::Get( )->SetDataOn( "GridCellvVel", gridCell.GetIndex( ), Environment::Get( "vVel", gridCell.GetIndex( ) ) ); // gridcell vVel
                DataRecorder::Get( )->SetDataOn( "GridCellDiurnalTemperatureRange", gridCell.GetIndex( ), Environment::Get( "DiurnalTemperatureRange", gridCell.GetIndex( ) ) ); // DiurnalTemperatureRange
                DataRecorder::Get( )->SetDataOn( "GridCellTotalPrecip", gridCell.GetIndex( ), Environment::Get( "TotalPrecip", gridCell.GetIndex( ) ) ); // gridcell total precipitation
                DataRecorder::Get( )->SetDataOn( "GridCellPrecipitation", gridCell.GetIndex( ), Environment::Get( "Precipitation", gridCell.GetIndex( ) ) ); // gridcell precipitation
                DataRecorder::Get( )->SetDataOn( "GridCellNPP", gridCell.GetIndex( ), Environment::Get( "NPP", gridCell.GetIndex( ) ) ); // gridcell NPP
                //###### end new code

                //######  additional megafauna specific cohorts statistics (turn on/off in OutputControlParameters.csv)

                // Mega carnivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEndothermBiomassDensity", gridCell.GetIndex( ), MCarnivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEndothermAbundance", gridCell.GetIndex( ), MCarnivoreEndothermAThisCell  / gridCell.GetIndex( ) ); 
                // Mega carnivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEctothermBiomassDensity", gridCell.GetIndex( ), MCarnivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaCarnivoreEctothermAbundance", gridCell.GetIndex( ), MCarnivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 

                // Mega Herbivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEndothermBiomassDensity", gridCell.GetIndex( ), MHerbivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEndothermAbundance", gridCell.GetIndex( ), MHerbivoreEndothermAThisCell  / gridCell.GetIndex( ) ); 
                // Mega Herbivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEctothermBiomassDensity", gridCell.GetIndex( ), MHerbivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaHerbivoreEctothermAbundance", gridCell.GetIndex( ), MHerbivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 
                    
                // Mega Omnivores endotherm
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEndothermBiomassDensity", gridCell.GetIndex( ), MOmnivoreEndothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEndothermAbundance", gridCell.GetIndex( ), MOmnivoreEndothermAThisCell  / gridCell.GetIndex( ) );   
                // Mega Omnivores ectotherm
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEctothermBiomassDensity", gridCell.GetIndex( ), MOmnivoreEctothermBMThisCell / gridCell.GetIndex( ) ); 
                DataRecorder::Get( )->SetDataOn( "MegaOmnivoreEctothermAbundance", gridCell.GetIndex( ), MOmnivoreEctothermAThisCell  / gridCell.GetIndex( ) ); 
                
                ////////////////////////////////////////
                // Store largest Carnivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassCarnivores", gridCell.GetIndex( ), MaxCarnivores );
                // Store largest Omnivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassOmnivores", gridCell.GetIndex( ), MaxOmnivores );
                // Store largest Herbivores
                DataRecorder::Get( )->SetDataOn( "MaxBodyMassHerbivores", gridCell.GetIndex( ), MaxHerbivores );
                /////////////////////////////////////

                //######  end additional megafauna specific cohorts statistics

    } );

    //# outputs included in MonthlyBasicOutputs.nc
    totalLivingBiomass = totalCohortBiomass + totalStockBiomass;
    totalBiomass = totalCohortBiomass + totalStockBiomass + respiratoryPool + organicMatterPool;

    DataRecorder::Get( )->SetDataOn( "InCellTime", mEcologyTimer.GetElapsedTimeSecs( ) );
    DataRecorder::Get( )->SetDataOn( "DispersalTime", mDispersalTimer.GetElapsedTimeSecs( ) );

    DataRecorder::Get( )->SetDataOn( "TotalBiomass", totalBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalLivingBiomass", totalLivingBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalStockBiomass", totalStockBiomass );
    DataRecorder::Get( )->SetDataOn( "TotalCohortBiomass", totalCohortBiomass );
    DataRecorder::Get( )->SetDataOn( "OrganicMatterPool", organicMatterPool );
    DataRecorder::Get( )->SetDataOn( "RespiratoryCO2Pool", respiratoryPool );

    DataRecorder::Get( )->SetDataOn( "NumberOfStocks", mGlobalDiagnosticVariables["NumberOfStocksInModel"] );
    DataRecorder::Get( )->SetDataOn( "NumberOfCohorts", mGlobalDiagnosticVariables["NumberOfCohortsInModel"] );

    DataRecorder::Get( )->SetDataOn( "CohortsProduced", mGlobalDiagnosticVariables["NumberOfCohortsProduced"] );
    DataRecorder::Get( )->SetDataOn( "CohortsExtinct", mGlobalDiagnosticVariables["NumberOfCohortsExtinct"] );
    DataRecorder::Get( )->SetDataOn( "CohortsCombined", mGlobalDiagnosticVariables["NumberOfCohortsCombined"] );
    DataRecorder::Get( )->SetDataOn( "CohortsDispersed", mDispersals );
    DataRecorder::Get( )->SetDataOn( "CohortAbundance", totalCohortAbundance );

}
