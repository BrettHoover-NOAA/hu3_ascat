Repositories:

./production
This repository contains the code-changes necessary to read the HU3-ASCAT BUFR file, select
an optimal wind from the available ambiguities (if present), and send appropriate data to
cdata_all for further processing in setupw with HU3-ASCAT winds defined as Type=290 and
Subtype=6 as an assimilated wind:

read_obs --> read_hu3ascat --> setupw
* Requires hard-wiring the global_convinfo.txt file found in production/fix/in config.anal
* Requires extra qc in setupw.f90 to reject HU3ASCAT winds with pre < 950. hPa

./lkcs_diag
This repository achieves the monitoring of HU3-ASCAT winds that is achieved in ./production,
but in addition:
 - HU3-ASCAT optimal wind likelihood (LKCS) is retained for the diag_conv_uv*.nc4 file by
   being passed to setupw through the slot reserved for station elevation. The station
   elevation MUST be set to zero for proper handling in setupw, so an if-trap is introduced
   in setupw to identify when cdata(ielev,i) comes from an HU3-ASCAT wind, at which point
   the local station elevation variable (dstn) is reassigned to zero. This will allow for
   proper handling of the HU3-ASCAT winds in setupw while letting the LKCS value still pass
   through as cdata(ielev,i) to the diag file under the Station_Elevation variable.

read_obs --> read_hu3ascat --> setupw
* Requires hard-wiring the global_convinfo.txt file found in production/fix/in config.anal
* Diag file contains HU3-ASCAT wind likelihood in Station_Elevation variable for HU3-ASCAT
  winds only, with the understanding that the *real* station elevation is zero.
