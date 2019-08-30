#include "UtilityFunctions.h"

void UtilityFunctions::ConvertToM180To180( std::vector<double>& lons ) {
    // Loop over longitudinal coordinates of the model grid cells
    for( unsigned jj = 0; jj < lons.size( ); jj++ ) {
        // If longitudinal coordindates exceed 180, then subtract 360 to correct the coorindates
        if( lons[jj] >= 180.0 ) {
            lons[jj] -= 360.0;
        }
    }
    // Re-sort the longitudinal coordinates
    sort( lons.begin( ), lons.end( ) );
}

std::vector<unsigned> UtilityFunctions::RandomlyOrderedCohorts( unsigned cohortNumber ) {
    //A vector to hold indices of cohorts in order
    std::vector<unsigned> RandomOrderCohorts( cohortNumber );
    for( unsigned i = 0; i < cohortNumber; i++ ) RandomOrderCohorts[i] = i;
    // Return the randomly ordered vector of cohort indices - randomly off system clock...c++11 style
    unsigned seed = std::chrono::system_clock::now( ).time_since_epoch( ).count( );
    std::shuffle( RandomOrderCohorts.begin( ), RandomOrderCohorts.end( ), std::default_random_engine( seed ) );
    return RandomOrderCohorts;
}

std::vector<unsigned> UtilityFunctions::NonRandomlyOrderedCohorts( unsigned cohortNumber, unsigned currentTimeStep ) {
    //A vector to hold indices of cohorts in order
    std::vector<unsigned> RandomOrderCohorts( cohortNumber ), OrderedCohorts( cohortNumber );
    for( unsigned i = 0; i < cohortNumber; i++ ) {
        RandomOrderCohorts[i] = i;
        OrderedCohorts[i] = i;
    }
    NonStaticSimpleRNG randomizer;
    randomizer.SetSeed( currentTimeStep );
    // Loop over cohorts
    for( int ii = 0; ii < OrderedCohorts.size( ); ii++ ) {
        // Generate a pseudo-random integer to swap this cohort index with
        int SwapIndex = ( int )( randomizer.GetUniform( )* ( OrderedCohorts.size( ) - 1 ) + 0.5 ); // random.Next(ii, OrderedCohorts.Length);
        // If the cohort index to swap is not the same as the active cohort index, then swap the values
        if( SwapIndex != ii ) {
            unsigned Temp = RandomOrderCohorts[ii];
            RandomOrderCohorts[ii] = RandomOrderCohorts[SwapIndex];
            RandomOrderCohorts[SwapIndex] = Temp;
        }
    }
    return RandomOrderCohorts;
}

unsigned UtilityFunctions::GetCurrentMonth( unsigned currentTimestep ) {
    unsigned Month;

    double DaysInYear = 360.0;
    double MonthsInYear = 12.0;
    double DaysInWeek = 7.0;

    //C++ can only use integers or enums in a switch, so make a map to convert.
    std::map<std::string, int> units;

    units["year" ] = 0;
    units["month"] = 1;
    units["week" ] = 2;
    units["day" ] = 3;

    //lowercase the string - a bit clunky...but then C++ strings are a bit
    switch( units[ Parameters::Get( )->GetTimeStepUnits( ) ] ) {
        case 0://year
            Month = 0;
            break;
        case 1://month
            Month = currentTimestep % 12;
            break;
        case 2://week
            Month = ( unsigned )floor( currentTimestep / ( ( DaysInYear / MonthsInYear ) / DaysInWeek ) ) % 12;
            break;
        case 3://day
            Month = ( unsigned )floor( currentTimestep / ( DaysInYear / MonthsInYear ) ) % 12;
            break;
        default://should the program bomb out at this point?
            std::cout << "Requested model time units not currently supported" << std::endl;
            Month = 100;
            break;

    }
    return Month;
}

double UtilityFunctions::ConvertTimeUnits( std::string fromUnit, std::string toUnit ) {
    //transform( fromUnit.begin( ), fromUnit.end( ), fromUnit.begin( ), ::tolower );
    toUnit = Convertor::Get( )->ToLowercase( toUnit );
    // Variable to hold the conversion factor
    double ConversionValue;
    double DaysInYear = 360.0;
    double MonthsInYear = 12.0;
    double DaysInWeek = 7.0;

    //C++ can only use intergers or enums in a switch, so make a map to convert.
    std::map<std::string, int> units;

    units["year"] = 0;
    units["month"] = 1;
    units["bimonth"] = 2;
    units["week"] = 3;
    units["day"] = 4;
    units["second"] = 5;

    // Determine which combination of time units is being requested and return the appropriate scaling factor
    switch( units[fromUnit] ) {
        case 0:// "year":
            switch( units[toUnit] ) {
                case 0:// "year":
                    ConversionValue = 1.0;
                    break;
                case 1://"month":
                    ConversionValue = MonthsInYear;
                    break;
                case 2://"bimonth":
                    ConversionValue = MonthsInYear * 2;
                    break;
                case 3://"week":
                    ConversionValue = DaysInYear / DaysInWeek;
                    break;
                case 4://"day":
                    ConversionValue = DaysInYear;
                    break;
                default:
                    std::cout << "Requested combination of time units not currently supported" << std::endl;
                    ConversionValue = 0;
                    break;
            }
            break;
        case 1://"month":
            switch( units[toUnit] ) {
                case 0:// "year":
                    ConversionValue = 1.0 / MonthsInYear;
                    break;
                case 1://"month":
                    ConversionValue = 1.0;
                    break;
                case 2://"bimonth":
                    ConversionValue = 2.0;
                    break;
                case 3://"week":
                    ConversionValue = ( DaysInYear / MonthsInYear ) / DaysInWeek;
                    break;
                case 4://"day":
                    ConversionValue = ( DaysInYear / MonthsInYear );
                    break;
                case 5://"second":
                    ConversionValue = ( DaysInYear / MonthsInYear ) * 24.0 * 60.0 * 60.0;
                    break;
                default:
                    std::cout << "Requested combination of time units not currently supported" << std::endl;
                    ConversionValue = 0;
                    break;
            }
            break;
        case 2://"bimonth":
            switch( units[toUnit] ) {
                case 0:// "year":
                    ConversionValue = 1.0 / ( MonthsInYear * 2 );
                    break;
                case 1://"month":
                    ConversionValue = 1 / 2.0;
                    break;
                case 2://"bimonth":
                    ConversionValue = 1.0;
                    break;
                case 3://"week":
                    ConversionValue = ( DaysInYear / ( MonthsInYear * 2 ) ) / DaysInWeek;
                    break;
                case 4://"day":
                    ConversionValue = ( DaysInYear / ( MonthsInYear * 2 ) );
                    break;
                case 5://"second":
                    ConversionValue = ( DaysInYear / ( MonthsInYear * 2 ) ) * 24.0 * 60.0 * 60.0;
                    break;
                default:
                    std::cout << "Requested combination of time units not currently supported" << std::endl;
                    ConversionValue = 0;
                    break;
            }
            break;

        case 3: //"week":
            switch( units[toUnit] ) {
                case 0:// "year":
                    ConversionValue = DaysInWeek / DaysInYear;
                    break;
                case 1://"month":
                    ConversionValue = DaysInWeek / ( DaysInYear / MonthsInYear );
                    break;
                case 2://"bimonth":
                    ConversionValue = DaysInWeek / ( DaysInYear / ( MonthsInYear * 2 ) );
                    break;
                case 3://"week":
                    ConversionValue = 1.0;
                    break;
                case 4://"day":
                    ConversionValue = DaysInWeek;
                    break;
                case 5://"second":
                    ConversionValue = DaysInWeek * 24.0 * 60.0 * 60.0;
                    break;
                default:
                    std::cout << "Requested combination of time units not currently supported" << std::endl;
                    ConversionValue = 0;
                    break;
            }
            break;
        case 4: //"day":
            switch( units[toUnit] ) {
                case 0:// "year":
                    ConversionValue = 1.0 / DaysInYear;
                    break;
                case 1://"month":
                    ConversionValue = 1.0 / ( DaysInYear / MonthsInYear );
                    break;
                case 2://"bimonth":
                    ConversionValue = 1.0 / ( DaysInYear / ( MonthsInYear * 2 ) );
                    break;
                case 3://"week":
                    ConversionValue = 1.0 / DaysInWeek;
                    break;
                case 4://"day":
                    ConversionValue = 1.0;
                    break;
                default:
                    std::cout << "Requested combination of time units not currently supported" << std::endl;
                    ConversionValue = 0;
                    break;
            }
            break;
        default:
            std::cout << "Requested combination of time units not currently supported" << std::endl;
            ConversionValue = 0;
            break;
    }

    // Return the conversion factor
    return ConversionValue;
}

std::vector<int> UtilityFunctions::FindJaggedArrayIndex( unsigned valueToFind, std::vector< std::vector < unsigned> > arrayToSearch, unsigned totalNumberOfCohorts ) {
    // Create a vector to hold the location of the cohort in the jagged array
    std::vector<int> ValueLocation( 2 );

    // Check to make sure that specified cohort index is not greater than the total number of cohorts
    assert( valueToFind < totalNumberOfCohorts && "Value searched for in jagged array is bigger than the biggest value in the jagged array" );

    // Variables to hold the row and colum indices of the cohort in the jaggged array
    int RowIndex = 0;
    int ColumnIndex = 0;

    // Loop over rows (functional groups) and locate the one in which the specified cohort is located
    while( arrayToSearch[RowIndex].size( ) == 0 || valueToFind > arrayToSearch[RowIndex][arrayToSearch[RowIndex].size( ) - 1] ) {
        RowIndex++;
    }

    // Add the located row to the vector of values to return
    ValueLocation[0] = RowIndex;

    // Loop over columns (cohorts within the functional group) and locate the one in which the specified cohort is located
    while( valueToFind != arrayToSearch[RowIndex][ColumnIndex] ) {
        ColumnIndex++;
    }

    // Add the located column to the vector of values to return
    ValueLocation[1] = ColumnIndex;

    // Return the vector of two values correpsonding to the located position in the jagged array of grid cell cohorts
    return ValueLocation;
}

double UtilityFunctions::ConvertSqMToSqDegrees( double valueToConvert, double latitude ) {
    // Convert the value to per square degree using the cosine of latitude and assuming cell dimensions of 110km by 110km at the Equator
    return valueToConvert * 110000.0 * 110000.0 * cos( DegreesToRadians( latitude ) );
}

double UtilityFunctions::LogNormalPDF( double xValue, double meanIdentity, double standardDeviation ) {
    const double PI = acos( -1. );
    // Calculate the mean of the log-normal distribution in log space
    double meanLog = log( meanIdentity );
    // Calculate and return the probability of the specified value under the specified log-normal distribution
    return (1 / sqrt( 2 * PI * pow( standardDeviation, 2 ) ) )*exp( -( pow( log( xValue ) - meanLog, 2 ) / ( 2 * pow( standardDeviation, 2 ) ) ) );
}

double UtilityFunctions::NormalPDF( double xValue, double meanValue, double standardDeviation ) {
    const double PI = acos( -1. );
    // Calculate and return the probability of the specified value under the specified normal distribution
    return (1 / sqrt( 2 * PI * pow( standardDeviation, 2 ) ) ) * exp( -( pow( xValue - meanValue, 2 ) / ( 2 * pow( standardDeviation, 2 ) ) ) );
}

double UtilityFunctions::CalculateGridCellArea( double latitude, double cellSize ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Angular eccentricity
    double AngularEccentricity = acos( DegreesToRadians( PolarRadius / EquatorialRadius ) );

    // First eccentricity squared
    //double ESquared = pow(sin(DegreesToRadians(AngularEccentricity)), 2);

    // Flattening
    //double Flattening = 1 - cos(DegreesToRadians(AngularEccentricity));

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Meridional radius of curvature
    double MPhi = pow( EquatorialRadius * PolarRadius, 2 ) / pow( TempVal, 1.5 );

    // Normal radius of curvature
    double NPhi = pow( EquatorialRadius, 2 ) / sqrt( TempVal );

    // Length of latitude (km)
    double LatitudeLength = PI / 180 * MPhi / 1000;

    // Length of longitude (km)
    double LongitudeLength = PI / 180 * cos( latitudeRad ) * NPhi / 1000;

    // Return the cell area in km^2
    return LatitudeLength * cellSize * LongitudeLength * cellSize;
}

double UtilityFunctions::CalculateLengthOfDegreeLatitude( float latitude ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Angular eccentricity
    //double AngularEccentricity = acos(DegreesToRadians(PolarRadius / EquatorialRadius));

    // First eccentricity squared
    //double ESquared = pow(sin(DegreesToRadians(AngularEccentricity)), 2);

    // Flattening
    //double Flattening = 1 - cos(DegreesToRadians(AngularEccentricity));

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Meridional radius of curvature
    double MPhi = pow( EquatorialRadius * PolarRadius, 2 ) / pow( TempVal, 1.5 );

    // Length of latitude (km)
    return PI / 180 * MPhi / 1000;
}

double UtilityFunctions::CalculateLengthOfDegreeLongitude( float latitude ) {
    const double PI = acos( -1. );
    // Convert from degrees to radians
    double latitudeRad = DegreesToRadians( latitude );

    // Equatorial radius in metres
    double EquatorialRadius = 6378137;

    // Polar radius in metres
    double PolarRadius = 6356752.3142;

    // Temporary value to save computations
    double TempVal = pow( ( EquatorialRadius * cos( latitudeRad ) ), 2 ) + pow( ( PolarRadius * sin( latitudeRad ) ), 2 );

    // Normal radius of curvature
    double NPhi = pow( EquatorialRadius, 2 ) / sqrt( TempVal );

    // Length of longitude (km)
    return PI / 180 * cos( latitudeRad ) * NPhi / 1000;
}

double UtilityFunctions::DegreesToRadians( double degrees ) {
    const double PI = acos( -1. );
    return (degrees * PI / 180.0 );
}

std::pair<double, double> UtilityFunctions::NewLocationGivenDistance( Location location, double distance, double bearing){
    const double R = 6371; //radius of earth in km
    const double pi = acos(-1);
    const double lat1 = DegreesToRadians(location.mLatitude);
    const double long1 = DegreesToRadians(location.mLongitude);
    
    
    double lat2 = asin( sin(lat1)*cos(distance/R) + cos(lat1)*sin(distance/R)*cos(bearing));
    double long2 = long1 + atan2(sin(bearing)*sin(distance/R)*cos(lat1), cos(distance/R) - sin(lat1)*sin(lat2));
    lat2 = lat2 * 180/pi;
    long2 = long2 * 180/pi;
    
    std::pair <double, double> coords = std::make_pair(lat2, long2);
    return coords;
}

unsigned UtilityFunctions::FindGridCellGivenLongandLat( double lat, double lon ){
    float maxLat = Parameters::Get()->GetUserMaximumLatitude();
    float minLat = Parameters::Get()->GetUserMinimumLatitude();
    float maxLon = Parameters::Get()->GetUserMaximumLongitude();
    float minLon = Parameters::Get()->GetUserMinimumLongitude();
    
    if (lat < maxLat && lat > minLat) ;
    else {
        if (lat > maxLat) {lat = maxLat;}
        if (lat < minLat) {lat = minLat;}
    }
    
    if (lon < maxLon && lon > minLon) ;
    else { 
        if (lon > maxLon) {lon = maxLon;}
        if (lon < minLon) {lon = minLon;}
    }
    float* longs = Parameters::Get()->GetUserLongitudeArray();
    float* lats = Parameters::Get()->GetUserLatitudeArray();
    float halfGridCell = (float)Parameters::Get()->GetGridCellSize()/2;
    unsigned longIndex, latIndex;  
    
    for (int a = 0; a <= Parameters::Get()->GetLengthUserLongitudeArray(); a++){
        if (longs[a] - halfGridCell <= lon && lon < longs[a] + halfGridCell) { longIndex = a; break; }
    }
    for (int b = 0; b <= Parameters::Get()->GetLengthUserLatitudeArray(); b++){
        if (lats[b] - halfGridCell <= lat && lat < lats[b] + halfGridCell) { latIndex = b; break; } 
    }
    return (Parameters::Get()->GetLengthUserLongitudeArray() * latIndex) + longIndex;
}

double UtilityFunctions::DistanceBetweenLocations( Location loc1, Location loc2 ){
    const double R = 6371; //radius of earth in km
    const double lat1 = DegreesToRadians(loc1.mLatitude);
    const double lat2 = DegreesToRadians(loc2.mLatitude);
    const double deltaLat = DegreesToRadians(loc1.mLatitude) - DegreesToRadians(loc2.mLatitude);
    const double deltaLong = DegreesToRadians(loc1.mLongitude) - DegreesToRadians(loc2.mLongitude);
   
    double a = pow(sin(deltaLat/2),2) + cos(lat1)*cos(lat2)*pow(sin(deltaLong/2), 2);
    double c = 2*atan2(sqrt(a), sqrt(1-a));   
    return R*c;
}


double UtilityFunctions::InitialBearingBetweenLocations( Location loc1, Location loc2 ){
    const double lat1 = DegreesToRadians(loc1.mLatitude);
    const double lat2 = DegreesToRadians(loc2.mLatitude);
    const double deltaLong = DegreesToRadians(loc1.mLongitude) - DegreesToRadians(loc2.mLongitude);
   
    double y = sin(deltaLong) * cos(lat2);
    double x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(deltaLong);
    double bearing = atan2(y, x);
    
    return bearing;
}

double UtilityFunctions::FindCircleOverlapGivenCentres( double r1, double r2, Location loc1, Location loc2 ) {
    const double pi = acos(-1);
    double d = DistanceBetweenLocations(loc1, loc2);
    double segmentArea =  0;
        
    // checking if there's no overlap or total overlap
    if (d >= r1 + r2 || r1 == 0 || r2 == 0) { return 0; }
    else {
        if ((d == 0 && r1>0 && r2>0) || d + r1 < r2 || d + r2 < r1) { return 1; }
        else {
            double x1 = abs((pow(d, 2) + pow(r1, 2) - pow(r2, 2)) / (2 * d)); //distance centre to line of intersection for centre 1
            double x2 = abs((pow(d, 2) + pow(r2, 2) - pow(r1, 2)) / (2 * d)); //distance centre to line of intersection for centre 2
            double a = sqrt((r1 - r2 - d)*(r2 - r1 - d)*(r1 + r2 - d)*(d + r1 + r2)) / d; //length of line of intersection
                
            segmentArea = atan(a / (2 * x1)) * pow(r1, 2) - 0.5*x1*a; //overlap area 1
            segmentArea += atan(a / (2 * x2)) * pow(r2, 2) - 0.5*x2*a; //overlap area 2
                
            return (segmentArea / (pi*pow(r1, 2)));
        }
    }
};


