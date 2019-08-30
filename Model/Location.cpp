#include "Location.h"

Location::Location( ) {
    mLongitudeIndex = 0;
    mLatitudeIndex = 0;
}

Location::Location( const Location& location ) {
    mLongitudeIndex = location.mLongitudeIndex;
    mLatitudeIndex = location.mLatitudeIndex;
    mLatitude = location.mLatitude;
    mLongitude = location.mLongitude;
}

void Location::SetIndices( unsigned latitudeIndex, unsigned longitudeIndex ) {
    mLatitudeIndex = latitudeIndex;
    mLongitudeIndex = longitudeIndex;
}

void Location::SetCoordinates ( double latitude, double longitude) {
    float maxLat = Parameters::Get()->GetUserMaximumLatitude();
    float minLat = Parameters::Get()->GetUserMinimumLatitude();
    float maxLon = Parameters::Get()->GetUserMaximumLongitude();
    float minLon = Parameters::Get()->GetUserMinimumLongitude();
    
    if (latitude < maxLat && latitude > minLat) mLatitude = latitude;
    else {
        if (latitude > maxLat) mLatitude = maxLat;
        if (latitude < minLat) mLatitude = minLat;
    }
    
    if (longitude < maxLon && longitude > minLon) mLongitude = longitude;
    else { 
        if (longitude > maxLon) mLongitude = maxLon;
        if (longitude < minLon) mLongitude = minLon;
    }
}

