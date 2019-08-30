#include "Parameters.h"

#include "Constants.h"
#include "Convertor.h"
#include "Maths.h"
#include "Processor.h"
#include "DataCoords.h"
#include "DataIndices.h"

Types::ParametersPointer Parameters::mThis = NULL;

Types::ParametersPointer Parameters::Get( ) {
    if( mThis == NULL ) {
        mThis = new Parameters( );
    }
    return mThis;
}

Parameters::~Parameters( ) {

    delete[ ] mMonthlyTimeStepArray;
    delete[ ] mAnnualTimeStepArray;
    delete[ ] mDataLongitudeArray;
    delete[ ] mDataLatitudeArray;
    delete[ ] mUserLongitudeArray;
    delete[ ] mUserLatitudeArray;

    if( mThis != NULL ) {
        delete mThis;
    }
}

Parameters::Parameters( ) {
}

bool Parameters::Initialise( const Types::StringMatrix& rawInputParameterData ) {
    bool success = false;

    if( rawInputParameterData.size( ) > 0 ) {
        if( rawInputParameterData[0].size( ) == Constants::eParameterValue + 2 ) {
            for( unsigned rowIndex = 0; rowIndex < rawInputParameterData.size( ); ++rowIndex ) {
                
                std::string parameterName = Convertor::Get( )->RemoveWhiteSpace( Convertor::Get( )->ToLowercase( rawInputParameterData[ rowIndex ][ Constants::eParameterName ] ) );

                if( parameterName == "rootdatadirectory" ) SetRootDataDirectory( Convertor::Get( )->RemoveWhiteSpace( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "timestepunits" ) SetTimeStepUnits( Convertor::Get( )->RemoveWhiteSpace( Convertor::Get( )->ToLowercase( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) ) );
                else if( parameterName == "lengthofsimulationinyears" ) SetLengthOfSimulationInMonths( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "minimumlongitude" ) SetUserMinimumLongitude( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "maximumlongitude" ) SetUserMaximumLongitude( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "minimumlatitude" ) SetUserMinimumLatitude( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "maximumlatitude" ) SetUserMaximumLatitude( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "gridcellsize" ) SetGridCellSize( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "extinctionthreshold" ) SetExtinctionThreshold( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "maximumnumberofcohorts" ) SetMaximumNumberOfCohorts( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "planktonsizethreshold" ) SetPlanktonSizeThreshold( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "drawrandomly" ) SetDrawRandomly( Convertor::Get( )->RemoveWhiteSpace( Convertor::Get( )->ToLowercase( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) ) );
                else if( parameterName == "humannppscenariotype" ) SetHumanNPPScenarioType( Convertor::Get( )->RemoveWhiteSpace( Convertor::Get( )->ToLowercase( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) ) );
                else if( parameterName == "humannppextractionscale" ) SetHumanNPPExtractionScale( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "humannppscenarioduration" ) SetHumanNPPScenarioDuration( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "burninsteps" ) SetBurninSteps( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "impactsteps" ) SetImpactSteps( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "recoverysteps" ) SetRecoverySteps( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                //######################
                else if( parameterName == "runparallel" ) SetRunParallel( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "threadnumber" ) SetThreadNumber( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                
                //# require model ini preferences
                else if( parameterName == "applymodelspinup" ) SetApplyModelSpinup( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "cohortcsvlocation" ) SetCohortCSVLocation( Convertor::Get( )->RemoveWhiteSpace( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "stockcsvlocation" ) SetStockCSVLocation( Convertor::Get( )->RemoveWhiteSpace( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                //######################
                else if( parameterName == "timestepstartextinction" ) SetTimeStepStartExtinction( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "startbodymass" ) SetStartBodyMass( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "endbodymass" ) SetEndBodyMass( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "stepbodymass" ) SetStepBodyMass( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "selectcarnivores" ) SetSelectCarnivores( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "selectomnivores" ) SetSelectOmnivores( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                else if( parameterName == "selectherbivores" ) SetSelectHerbivores( Convertor::Get( )->StringToNumber( rawInputParameterData[ rowIndex ][ Constants::eParameterValue ] ) );
                //######################
            }
            CalculateParameters( );
            
            success = true;
        }
    }
    return success;
}

void Parameters::CalculateParameters( ) {

    // Calculate temporal parameters
    mLengthOfSimulationInMonths = mLengthOfSimulationInYears * 12;

    mMonthlyTimeStepArray = new unsigned[ mLengthOfSimulationInMonths ];
    for( unsigned monthIndex = 0; monthIndex < mLengthOfSimulationInMonths; ++monthIndex ) {
        mMonthlyTimeStepArray[ monthIndex ] = monthIndex;
    }

    mAnnualTimeStepArray = new unsigned[ mLengthOfSimulationInYears ];
    for( unsigned yearIndex = 0; yearIndex < mLengthOfSimulationInYears; ++yearIndex ) {
        mAnnualTimeStepArray[ yearIndex ] = yearIndex;
    }

    // Calculate spatial parameters
    mLengthDataLongitudeArray = 360 / mGridCellSize;
    mDataLongitudeArray = new float[ mLengthDataLongitudeArray ];
    for( unsigned longitudeIndex = 0; longitudeIndex < mLengthDataLongitudeArray; ++longitudeIndex ) {
        mDataLongitudeArray[ longitudeIndex ] = ( -180 + ( ( float )mGridCellSize / 2 ) ) + ( longitudeIndex * ( float )mGridCellSize );
    }

    mLengthDataLatitudeArray = 180 / mGridCellSize;
    mDataLatitudeArray = new float[ mLengthDataLatitudeArray ];
    for( unsigned latitudeIndex = 0; latitudeIndex < mLengthDataLatitudeArray; ++latitudeIndex ) {
        mDataLatitudeArray[ latitudeIndex ] = ( -90 + ( ( float )mGridCellSize / 2 ) ) + ( latitudeIndex * ( float )mGridCellSize );
    }

    mDataIndexOfUserMinimumLongitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLongitudeArray, mLengthDataLongitudeArray, mUserMinimumLongitude );
    mDataIndexOfUserMaximumLongitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLongitudeArray, mLengthDataLongitudeArray, mUserMaximumLongitude );
    mLengthUserLongitudeArray = ( mDataIndexOfUserMaximumLongitude - mDataIndexOfUserMinimumLongitude)+1; //+ 1;
    mLengthUserLongitudeArray2 = mLengthUserLongitudeArray -1;

    mUserLongitudeArray = new float[ mLengthUserLongitudeArray ];
    for( unsigned userLongitudeIndex = 0; userLongitudeIndex < mLengthUserLongitudeArray; ++userLongitudeIndex ) {
        mUserLongitudeArray[ userLongitudeIndex ] = mDataLongitudeArray[ userLongitudeIndex + mDataIndexOfUserMinimumLongitude ];
    }

    mDataIndexOfUserMinimumLatitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLatitudeArray, mLengthDataLatitudeArray, mUserMinimumLatitude );
    mDataIndexOfUserMaximumLatitude = Processor::Get( )->CalculateArrayIndexOfValue( mDataLatitudeArray, mLengthDataLatitudeArray, mUserMaximumLatitude );
    mLengthUserLatitudeArray = ( mDataIndexOfUserMaximumLatitude - mDataIndexOfUserMinimumLatitude )+1; //+ 1;
    mLengthUserLatitudeArray2 = mLengthUserLatitudeArray -1;

    mUserLatitudeArray = new float[ mLengthUserLatitudeArray ];
    for( unsigned userLatitudeIndex = 0; userLatitudeIndex < mLengthUserLatitudeArray; ++userLatitudeIndex ) {
        mUserLatitudeArray[ userLatitudeIndex ] = mDataLatitudeArray[ userLatitudeIndex + mDataIndexOfUserMinimumLatitude ];
    }

    mNumberOfGridCells = mLengthUserLongitudeArray  * mLengthUserLatitudeArray ;
    mSizeOfMonthlyGridDatum = mNumberOfGridCells * mLengthOfSimulationInMonths;
    mSizeOfAnnualGridDatum = mNumberOfGridCells * mLengthOfSimulationInYears;

    unsigned cellIndex = 0;
    mCoordsIndicesLookup.resize( mNumberOfGridCells );
    for( unsigned latitudeIndex = 0; latitudeIndex < mLengthUserLatitudeArray; ++latitudeIndex ) {
        for( unsigned longitudeIndex = 0; longitudeIndex < mLengthUserLongitudeArray; ++longitudeIndex ) {

            float longitude = mUserLongitudeArray[ longitudeIndex ];
            float latitude = mUserLatitudeArray[ latitudeIndex ];

            Types::DataCoordsPointer coords = new DataCoords( longitude, latitude );
            Types::DataIndicesPointer indices = new DataIndices( longitudeIndex, latitudeIndex );

            mCoordsIndicesLookup[ cellIndex ] = std::make_pair( coords, indices );

            cellIndex += 1;
        }
    }
}

std::string Parameters::GetRootDataDirectory( ) const {
    return mRootDataDirectory;
}

std::string Parameters::GetTimeStepUnits( ) const {
    return mTimeStepUnits;
}

unsigned Parameters::GetLengthOfSimulationInYears( ) const {
    return mLengthOfSimulationInYears;
}

int Parameters::GetUserMinimumLongitude( ) const {
    return mUserMinimumLongitude;
}

int Parameters::GetUserMaximumLongitude( ) const {
    return mUserMaximumLongitude;
}

int Parameters::GetUserMinimumLatitude( ) const {
    return mUserMinimumLatitude;
}

int Parameters::GetUserMaximumLatitude( ) const {
    return mUserMaximumLatitude;
}

unsigned Parameters::GetGridCellSize( ) const {
    return mGridCellSize;
}

float Parameters::GetExtinctionThreshold( ) const {
    return mExtinctionThreshold;
}

unsigned Parameters::GetMaximumNumberOfCohorts( ) const {
    return mMaximumNumberOfCohorts;
}

float Parameters::GetPlanktonSizeThreshold( ) const {
    return mPlanktonSizeThreshold;
}

bool Parameters::GetDrawRandomly( ) const {
    return mDrawRandomly;
}


std::string Parameters::GetHumanNPPScenarioType( ) const {
    return mHumanNPPScenarioType;
}
double Parameters::GetHumanNPPExtractionScale( ) const{
    return mHumanNPPExtractionScale;
}
double Parameters::GetHumanNPPScenarioDuration( ) const{
    return mHumanNPPScenarioDuration;
}
unsigned Parameters::GetBurninSteps( ) const{
    return mBurninSteps;
}
unsigned Parameters::GetImpactSteps( ) const{
    return mImpactSteps;
}
unsigned Parameters::GetRecoverySteps( ) const{
    return mRecoverySteps;
}
//############################################
unsigned Parameters::GetRunParallel( ) const{
    return mRunParallel;
}
unsigned Parameters::GetThreadNumber( ) const{
    return mThreadNumber;
}

//############################################
unsigned Parameters::GetApplyModelSpinup( ) const{
    return mApplyModelSpinup;
}
std::string Parameters::GetCohortCSVLocation( ) const{
    return mCohortCSVLocation;
}
std::string Parameters::GetStockCSVLocation( ) const{
    return mStockCSVLocation;
}
//############################################
unsigned Parameters::GetTimeStepStartExtinction( ) const{
    return mTimeStepStartExtinction;
}
unsigned Parameters::GetStartBodyMass( ) const{
    return mStartBodyMass;
}
unsigned Parameters::GetEndBodyMass( ) const{
    return mEndBodyMass;
}
unsigned Parameters::GetStepBodyMass( ) const{
    return mStepBodyMass;
}
unsigned Parameters::GetSelectCarnivores( ) const{
    return mSelectCarnivores;
}
unsigned Parameters::GetSelectOmnivores( ) const{
    return mSelectOmnivores;
}
unsigned Parameters::GetSelectHerbivores( ) const{
    return mSelectHerbivores;
}
//############################################

void Parameters::SetRootDataDirectory( const std::string& rootDataDirectory ) {
    mRootDataDirectory = rootDataDirectory;
}

void Parameters::SetTimeStepUnits( const std::string& timeStepUnits ) {
    mTimeStepUnits = timeStepUnits;
}

void Parameters::SetLengthOfSimulationInMonths( const unsigned& lengthOfSimulationInYears ) {
    mLengthOfSimulationInYears = lengthOfSimulationInYears;
}

void Parameters::SetUserMinimumLongitude( const int& userMinimumLongitude ) {
    mUserMinimumLongitude = userMinimumLongitude;
}

void Parameters::SetUserMaximumLongitude( const int& userMaximumLongitude ) {
    mUserMaximumLongitude = userMaximumLongitude;
}

void Parameters::SetUserMinimumLatitude( const int& userMinimumLatitude ) {
    mUserMinimumLatitude = userMinimumLatitude;
}

void Parameters::SetUserMaximumLatitude( const int& userMaximumLatitude ) {
    mUserMaximumLatitude = userMaximumLatitude;
}

void Parameters::SetGridCellSize( const unsigned& gridCellSize ) {
    mGridCellSize = gridCellSize;
}

void Parameters::SetExtinctionThreshold( const float& extinctionThreshold ) {
    mExtinctionThreshold = extinctionThreshold;
}

void Parameters::SetMaximumNumberOfCohorts( const unsigned& maximumNumberOfCohorts ) {
    mMaximumNumberOfCohorts = maximumNumberOfCohorts;
}

void Parameters::SetPlanktonSizeThreshold( const float& planktonSizeThreshold ) {
    mPlanktonSizeThreshold = planktonSizeThreshold;
}

void Parameters::SetDrawRandomly( const std::string& drawRandomlyString ) {
    if( drawRandomlyString == "yes" )
        mDrawRandomly = true;
    else
        mDrawRandomly = false;
}

void Parameters::SetHumanNPPScenarioType(const std::string& humanNPPScenarioType){
    mHumanNPPScenarioType=humanNPPScenarioType;
}
void Parameters::SetHumanNPPExtractionScale(const double& humanNPPExtractionScale ){
    mHumanNPPExtractionScale=humanNPPExtractionScale;
}
void Parameters::SetHumanNPPScenarioDuration(const double & humanNPPScenarioDuration){
    mHumanNPPScenarioDuration=humanNPPScenarioDuration;
}
void Parameters::SetBurninSteps(const unsigned& burninSteps){
    mBurninSteps=burninSteps;
}
void Parameters::SetImpactSteps(const unsigned& impactSteps){
    mImpactSteps=impactSteps;
}
void Parameters::SetRecoverySteps(const unsigned& recoverySteps){
    mRecoverySteps=recoverySteps;
}

//############################################
void Parameters::SetRunParallel(const unsigned& runParallel){
    mRunParallel=runParallel;
}
void Parameters::SetThreadNumber(const unsigned& threadNumber){
    mThreadNumber=threadNumber;
}

//############################################
void Parameters::SetApplyModelSpinup(const unsigned& applyModelSpinup){
    mApplyModelSpinup=applyModelSpinup;
}
void Parameters::SetCohortCSVLocation(const std::string& cohortCSVLocation){
    mCohortCSVLocation=cohortCSVLocation;
}
void Parameters::SetStockCSVLocation(const std::string& stockCSVLocation){
    mStockCSVLocation=stockCSVLocation;
}

//############################################
void Parameters::SetTimeStepStartExtinction(const unsigned& timeStepStartExtinction){
    mTimeStepStartExtinction=timeStepStartExtinction;
}
void Parameters::SetStartBodyMass(const unsigned& startBodyMass){
    mStartBodyMass=startBodyMass;
}
void Parameters::SetEndBodyMass(const unsigned& endBodyMass){
    mEndBodyMass=endBodyMass;
}
void Parameters::SetStepBodyMass(const unsigned& stepBodyMass){
    mStepBodyMass=stepBodyMass;
}
void Parameters::SetSelectCarnivores(const unsigned& selectCarnivores){
    mSelectCarnivores=selectCarnivores;
}
void Parameters::SetSelectHerbivores(const unsigned& selectHerbivores){
    mSelectHerbivores=selectHerbivores;
}
void Parameters::SetSelectOmnivores(const unsigned& selectOmnivores){
    mSelectOmnivores=selectOmnivores;
}
//############################################

unsigned Parameters::GetNumberOfGridCells( ) const {
    return mNumberOfGridCells;
}

unsigned Parameters::GetLengthOfSimulationInMonths( ) const {
    return mLengthOfSimulationInMonths;
}

unsigned Parameters::GetLengthDataLongitudeArray( ) const {
    return mLengthDataLongitudeArray;
}

unsigned Parameters::GetLengthDataLatitudeArray( ) const {
    return mLengthDataLatitudeArray;
}

unsigned Parameters::GetDataIndexOfUserMinimumLongitude( ) const {
    return mDataIndexOfUserMinimumLongitude;
}

unsigned Parameters::GetDataIndexOfUserMaximumLongitude( ) const {
    return mDataIndexOfUserMaximumLongitude;
}

unsigned Parameters::GetDataIndexOfUserMinimumLatitude( ) const {
    return mDataIndexOfUserMinimumLatitude;
}

unsigned Parameters::GetDataIndexOfUserMaximumLatitude( ) const {
    return mDataIndexOfUserMaximumLatitude;
}

unsigned Parameters::GetLengthUserLongitudeArray( ) const {
    return mLengthUserLongitudeArray;
}

unsigned Parameters::GetLengthUserLongitudeArray2( ) const {
    return mLengthUserLongitudeArray2;
}

unsigned Parameters::GetLengthUserLatitudeArray( ) const {
    return mLengthUserLatitudeArray;
}

unsigned Parameters::GetLengthUserLatitudeArray2( ) const {
    return mLengthUserLatitudeArray2;
}

unsigned Parameters::GetSizeOfMonthlyGridDatum( ) const {
    return mSizeOfMonthlyGridDatum;
}

unsigned Parameters::GetSizeOfAnnualGridDatum( ) const {
    return mSizeOfAnnualGridDatum;
}

float Parameters::GetDataLongitudeAtIndex( const unsigned& index ) const {
    return mDataLongitudeArray[ index ];
}

float Parameters::GetDataLatitudeAtIndex( const unsigned& index ) const {
    return mDataLatitudeArray[ index ];
}

float Parameters::GetUserLongitudeAtIndex( const unsigned& index ) const {
    return mUserLongitudeArray[ index ];
}

float Parameters::GetUserLatitudeAtIndex( const unsigned& index ) const {
    return mUserLatitudeArray[ index ];
}

//float* Parameters::GetDataLongitudeArray( ) const {
//    return mDataLongitudeArray;
//}
//
//float* Parameters::GetDataLatitudeArray( ) const {
//    return mDataLatitudeArray;
//}

unsigned* Parameters::GetMonthlyTimeStepArray( ) const {
    return mMonthlyTimeStepArray;
}

unsigned* Parameters::GetAnnualTimeStepArray( ) const {
    return mAnnualTimeStepArray;
}

float* Parameters::GetUserLongitudeArray( ) const {
    return mUserLongitudeArray;
}

float* Parameters::GetUserLatitudeArray( ) const {
    return mUserLatitudeArray;
}

int Parameters::GetCellIndexFromDataIndices( const unsigned& longitudeIndex, const unsigned& latitudeIndex ) const {

    int cellIndex = Constants::cMissingValue;
    for( unsigned index = 0; index < mNumberOfGridCells; ++index ) {
        Types::DataIndicesPointer indices = mCoordsIndicesLookup[ index ].second;

        if( indices->GetX( ) == longitudeIndex && indices->GetY( ) == latitudeIndex ) {
            cellIndex = index;
            break;
        }
    }
    return cellIndex;
}

Types::DataCoordsPointer Parameters::GetDataCoordsFromCellIndex( const unsigned& cellIndex ) const {
    return mCoordsIndicesLookup[ cellIndex ].first;
}

Types::DataIndicesPointer Parameters::GetDataIndicesFromCellIndex( const unsigned& cellIndex ) const {
    return mCoordsIndicesLookup[ cellIndex ].second;
}
