
############################
### Identity conversions ###
############################

ECEF(ecef::ECEF, datum) = ecef
LLA(lla::LLA, datum) = lla
ENU(enu::ENU, datum) = enu
UTMZ(utm::UTMZ, datum) = utm


################################
### LLA <-> ECEF coordinates ###
################################

ECEF(lla::LLA, datum) = transform(ECEFfromLLA(datum), lla)
LLA(ecef::ECEF, datum) = transform(LLAfromECEF(datum), ecef)


################################
### ECEF <-> ENU coordinates ###
################################

ENU(ecef::ECEF, origin_lla::LLA, datum) = transform(ENUfromECEF(origin_lla, datum), ecef)
ENU(ecef::ECEF, origin_ecef::ECEF, datum) = transform(ENUfromECEF(origin_ecef, datum), ecef)

ECEF(enu::ENU, origin_lla::LLA, datum) = transform(ECEFfromENU(origin_lla, datum), enu)
ECEF(enu::ENU, origin_ecef::ECEF, datum) = transform(ECEFfromENU(origin_ecef, datum), enu)


################################
### LLA <-> ENU coordinates ###
################################

ENU(lla::LLA, origin, datum) =  transform(ENUfromLLA(origin, datum), lla)
LLA(enu::ENU, origin, datum) = transform(LLAfromENU(origin, datum), enu)


################################
### LLA <-> UTMZ coordinates ###
################################

LLA(utm::UTMZ, datum) = transform(LLAfromUTMZ(datum), utm)
UTMZ(lla::LLA, datum) = transform(UTMZfromLLA(datum), lla)


#################################
### ECEF <-> UTMZ coordinates ###
#################################

ECEF(utm::UTMZ, datum) = transform(ECEFfromUTMZ(datum), utm)
UTMZ(ecef::ECEF, datum) = transform(UTMZfromECEF(datum), ecef)


################################
### ENU <-> UTMZ coordinates ###
################################

ENU(utm::UTMZ, origin, datum) = transform(ENUfromUTMZ(origin, datum), utm)
UTMZ(enu::ENU, origin, datum) = transform(UTMZfromENU(origin, datum), enu)


###############################
### LLA <-> UTM coordinates ###
###############################

LLA(utm::UTM, zone::Integer, hemisphere::Bool, datum) = transform(LLAfromUTM(zone, hemisphere, datum), utm)
UTM(lla::LLA, zone::Integer, hemisphere::Bool, datum) = transform(UTMfromLLA(zone, hemisphere, datum), lla)


################################
### ECEF <-> UTM coordinates ###
################################

ECEF(utm::UTM, zone::Integer, hemisphere::Bool, datum) = transform(ECEFfromUTM(zone, hemisphere, datum), utm)
UTM(ecef::ECEF, zone::Integer, hemisphere::Bool, datum) = transform(UTMfromECEF(zone, hemisphere, datum), lla)


###############################
### ENU <-> UTM coordinates ###
###############################

ENU(utm::UTM, zone::Integer, hemisphere::Bool, origin, datum) = transform(ENUfromUTM(origin, zone, hemisphere, datum), utm)
UTM(enu::ENU, zone::Integer, hemisphere::Bool, origin, datum) = transform(UTMfromENU(origin, zone, hemisphere, datum), enu)
