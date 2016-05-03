/* begin dune-dpg
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* Boost::Fusion definitions */

/* 7 would be enough, but better take some more, in case we
 * change our procedures later. */
#define BOOST_FUSION_INVOKE_PROCEDURE_MAX_ARITY 10

/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of dune-dpg */
#define DUNE_DPG_VERSION "@DUNE_DPG_VERSION@"

/* Define to the major version of dune-dpg */
#define DUNE_DPG_VERSION_MAJOR @DUNE_DPG_VERSION_MAJOR@

/* Define to the minor version of dune-dpg */
#define DUNE_DPG_VERSION_MINOR @DUNE_DPG_VERSION_MINOR@

/* Define to the revision of dune-dpg */
#define DUNE_DPG_VERSION_REVISION @DUNE_DPG_VERSION_REVISION@

/* Define if the Eigen3 headers were found */
#cmakedefine HAVE_EIGEN3 1

/* end dune-dpg
   Everything below here will be overwritten
*/
