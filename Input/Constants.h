#ifndef CONSTANTS
#define CONSTANTS

#include <string>

namespace Constants {
    
    enum eParametersMetadata {
        eParameterName,
        eParameterValue
    };

    enum eOutputControlParametersMetadata {
        eDatumName,
        eDatumType,
        eTimeUnit,
        eDataUnit
    };

    enum eEnvironmentalDataLayersMetadata {
        eInternalName,
        eFilePath,
        eDefaultVariableName
    };

    enum eDataLayerTypes {
        eDataLayer2D, // Spatial: land-sea mask
        eDataLayer2DwithTime, // Two-dimensional with time, e.g. SST.
        eDataLayer3D, // Three-dimensional, e.g. ocean temperature.
        eDataLayer3DwithTime // Three-dimensional with time.
    };

    enum eVariableTypes {
        eLongitude,
        eLatitude,
        eTime,
        eDepth,
        eOther
    };

    const std::string cLongitudeVariableNames[ ] = { "lon", "long", "longitude", "x" };
    const std::string cLatitudeVariableNames[ ] = { "lat", "lats", "latitude", "y" };
    const std::string cDepthVariableNames[ ] = { "dep", "depth", "z" };
    const std::string cTimeVariableNames[ ] = { "time", "month", "t" };

    const std::string cBasicDatumTypeName = "basic";
    const std::string cGridDatumTypeName = "grid";

    const std::string cConfigurationDirectory = "C://Users/science.intern2/Desktop/Madingley-cp/src/Model_setup/input/Model-setup/";
    const std::string cInputParametersFileName = "SimulationControlParameters.csv";
    const std::string cInputDataFileName = "EnvironmentalDataLayers.csv";
    const std::string cOutputVariablesFileName = "OutputControlParameters.csv";

    const std::string cCohortDefinitionsFileName = "CohortFunctionalGroupDefinitions.csv";
    const std::string cStockDefinitionsFileName = "StockFunctionalGroupDefinitions.csv";
    const std::string cMassBinDefinitionsFileName = "MassBinDefinitions.csv";

    const std::string cOutputBaseDirectory = "C://Users/science.intern2/Desktop/output/";
    const std::string cDataSetNameFormat = "%Y-%m-%d_%H-%M-%S";
    const std::string cCompleteDateFormat = "%c";
    const std::string cUnitsString = "units";
    const std::string cTimeVariableUnits = "months";

    ///////////////////////////////////////////////////////////////////////////
    const std::string cAnnualBasicOutputsFileName = "AnnualBasicOutputs.nc4";
    const std::string cMonthlyBasicOutputsFileName = "MonthlyBasicOutputs.nc4";
    const std::string cAnnualGridOutputsFileName = "AnnualGridOutputs.nc4";
    const std::string cMonthlyGridOutputsFileName = "MonthlyGridOutputs.nc4";
    const std::string cMonthlyTimeUnitName = "month";
    const std::string cAnnualTimeUnitName = "year";
    const std::string cLongitudeVariableUnit = "degrees east";
    const std::string cLatitudeVariableUnit = "degrees north";
    ///////////////////////////////////////////////////////////////////////////

    const int cMissingValue = -9999;
    const int cOutputFolderPermissions = 0775;
    const int cDateTimeBufferSize = 25;

    const char cDataDelimiterValue = ',';
    const char cTextFileCommentCharacter = '#';
    const char cFolderDelimiter = '/';
    const char cWhiteSpaceCharacter = ' ';
    
    const bool cHibernationEnabled = false;
    const float cminLatforHib = 37;
    const float cHibernationSlowdown = 0.5; //Dispersal speed reduced to 50% of normal
    const float cWinterSlowDown = 0.05; // time active reduced to 5% of normal
    const float cHibernationMetabolism = 0.05; // metabolism slowed to 5% of normal
    const float cHibernatingFoodSlowDown = 0.05; // cant search for food all the time if you're hibernating
    const float cHibernatingProtectionInADen = 0.2; // harder to be preyed on if you're in a den
    const float cStarvationThreshold = 0.4; // fraction of maximum body weight at which starvation deaths kick in
}

#endif
