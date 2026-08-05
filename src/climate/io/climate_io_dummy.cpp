/*
 * Placeholder translation unit for the climate NetCDF convenience library.
 *
 * The climate grid readers/writers in this directory require NetCDF or
 * PNetCDF.  When neither is configured, libmoabclimateio still has to exist
 * so that the link in src/Makefile.am resolves; this file gives it a single
 * (empty) object so automake can infer a linker language and so that the
 * resulting archive is not empty, which some archivers reject.
 */

namespace moab
{
namespace climate
{
// Intentionally empty; see the file comment above.
}  // namespace climate
}  // namespace moab
