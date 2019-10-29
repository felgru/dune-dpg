#
# Module that checks whether boost::hana is available and usable.
#
# Sets the follwing variable:
#
# HAVE_BOOST_HANA True if boost::hana is available.
#
include(DuneBoost)

if(HAVE_DUNE_BOOST)
  if(${BoostHana_FIND_VERSION_MAJOR})
    if(${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} VERSION_LESS
        ${BoostHana_FIND_VERSION_MAJOR}.${BoostHana_FIND_VERSION_MINOR})
      set(BoostHana_FOUND false)
      if(${BoostHana_FIND_REQUIRED})
        message(SEND_ERROR "Boost::Hana ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION} found, but ${BoostHana_FIND_VERSION_MAJOR}.${BoostHana_FIND_VERSION_MINOR} required")
      endif()
      return()
    endif()
  endif()
  message(STATUS "Checking whether the Boost::Hana library is available.")
  set(OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_INCLUDES ${BOOST_INCLUDEDIR})
  check_cxx_source_compiles("
\#include <boost/hana.hpp>
int main(){
  boost::hana::tuple<int,char,double> v;
  return 0;
}" HAVE_BOOST_HANA)
  set(CMAKE_REQUIRED_INCLUDES OLD_CMAKE_REQUIRED_INCLUDES)
  if(HAVE_BOOST_HANA)
    message(STATUS "Boost::Hana is available")
    set(BoostHana_FOUND true)
  endif(HAVE_BOOST_HANA)
else()
  message(STATUS "Skipping check for Boost::Hana as Boost is not available.")
endif()
