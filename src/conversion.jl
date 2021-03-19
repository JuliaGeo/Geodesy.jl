
############################
### Identity conversions ###
############################

ECEF(ecef::ECEF, datum) = ecef
LLA(lla::LLA, datum) = lla
ENU(enu::ENU, datum) = enu
NED(ned::NED, datum) = ned
UTM(utm::UTM, datum) = utm
UTMZ(utmz::UTMZ, datum) = utmz


################################
### ENU <-> NED coordinates  ###
################################

ENU(ned::NED) = ENUfromNED()(ned)
NED(enu::ENU) = NEDfromENU()(enu)


################################
### LLA <-> ECEF coordinates ###
################################

ECEF(lla::LLA, datum) = ECEFfromLLA(datum)(lla)
LLA(ecef::ECEF, datum) = LLAfromECEF(datum)(ecef)


################################
### ECEF <-> ENU coordinates ###
################################

ENU(ecef::ECEF, origin, datum) = ENUfromECEF(origin, datum)(ecef)
ECEF(enu::ENU, origin, datum) = ECEFfromENU(origin, datum)(enu)


################################
### ECEF <-> NED coordinates ###
################################

NED(ecef::ECEF, origin, datum) = NEDfromECEF(origin, datum)(ecef)
ECEF(ned::NED, origin, datum) = ECEFfromNED(origin, datum)(ned)


################################
### LLA <-> ENU coordinates ###
################################

ENU(lla::LLA, origin, datum) =  ENUfromLLA(origin, datum)(lla)
LLA(enu::ENU, origin, datum) = LLAfromENU(origin, datum)(enu)


################################
### LLA <-> NED coordinates ###
################################

NED(lla::LLA, origin, datum) =  NEDfromLLA(origin, datum)(lla)
LLA(ned::NED, origin, datum) = LLAfromNED(origin, datum)(ned)


################################
### LLA <-> UTMZ coordinates ###
################################

LLA(utm::UTMZ, datum) = LLAfromUTMZ(datum)(utm)
UTMZ(lla::LLA, datum) = UTMZfromLLA(datum)(lla)


#################################
### ECEF <-> UTMZ coordinates ###
#################################

ECEF(utm::UTMZ, datum) = ECEFfromUTMZ(datum)(utm)
UTMZ(ecef::ECEF, datum) = UTMZfromECEF(datum)(ecef)


################################
### ENU <-> UTMZ coordinates ###
################################

ENU(utm::UTMZ, origin, datum) = ENUfromUTMZ(origin, datum)(utm)
UTMZ(enu::ENU, origin, datum) = UTMZfromENU(origin, datum)(enu)


################################
### NED <-> UTMZ coordinates ###
################################

NED(utm::UTMZ, origin, datum) = NEDfromUTMZ(origin, datum)(utm)
UTMZ(ned::NED, origin, datum) = UTMZfromNED(origin, datum)(ned)


###############################
### LLA <-> UTM coordinates ###
###############################

LLA(utm::UTM, zone::Integer, hemisphere::Bool, datum) = LLAfromUTM(zone, hemisphere, datum)(utm)
UTM(lla::LLA, zone::Integer, hemisphere::Bool, datum) = UTMfromLLA(zone, hemisphere, datum)(lla)


################################
### ECEF <-> UTM coordinates ###
################################

ECEF(utm::UTM, zone::Integer, hemisphere::Bool, datum) = ECEFfromUTM(zone, hemisphere, datum)(utm)
UTM(ecef::ECEF, zone::Integer, hemisphere::Bool, datum) = UTMfromECEF(zone, hemisphere, datum)(ecef)


###############################
### ENU <-> UTM coordinates ###
###############################

ENU(utm::UTM, zone::Integer, hemisphere::Bool, origin, datum) = ENUfromUTM(origin, zone, hemisphere, datum)(utm)
UTM(enu::ENU, zone::Integer, hemisphere::Bool, origin, datum) = UTMfromENU(origin, zone, hemisphere, datum)(enu)


###############################
### NED <-> UTM coordinates ###
###############################

NED(utm::UTM, zone::Integer, hemisphere::Bool, origin, datum) = NEDfromUTM(origin, zone, hemisphere, datum)(utm)
UTM(ned::NED, zone::Integer, hemisphere::Bool, origin, datum) = UTMfromNED(origin, zone, hemisphere, datum)(ned)


###############################
### UTMZ <-> UTM coordinates ###
###############################

UTMZ(utm::UTM, zone::Integer, hemisphere::Bool, datum) = UTMZfromUTM(zone, hemisphere, datum)(utm)
UTM(utmz::UTMZ, zone::Integer, hemisphere::Bool, datum) = UTMfromUTMZ(zone, hemisphere, datum)(utmz)
