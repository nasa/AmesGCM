# NASA Ames Mars Global Climate Model
## Description of the default simulation
The Ames Mars GCM is released with **fms\_mars_default**, the default run script.  
*    The default simulation has a horizontal simulation of 24x24 cells per cube face, or approximately 240x240 km grid cell size.  
*    The vertical resolution is 30 vertical layers with a top at 0.02 Pa.  
*    The surface albedo are derived from TES observations, with thermal inertia values derived to produce observed surface temperature diurnal cycles.  
*    CO<sub>2</sub> is condensed at the surface when predicted temperatures fall below the condensation temperature. The condensed mass is removed from the atmosphere. Soil properties are tuned to produce observed annual pressure cycles from the Viking Landers. Atmospheric CO<sub>2</sub> is condensed when the temperature falls below the condensation temperature. Mass is removed from the layer where condensation occurred, then deposited on the surface.   
*    For radiative heating purposes, the default simulation assumes a fixed dust distribution that matches a background dust scenario representing a minimum dust climatology over 10 Mars Years. The dust field assumes a modified Conrath profile in the vertical, scaled to a horizontally and temporally varying prescribed top level.
*    Dust is lifted and advected as passive tracers of moment mass and number with a log-normal size distribution. Lifting is performed as an assimilation, reading the background dust scenario and lifting according to the column opacity deficit compared to the map. The dust scenario is a series of 2D lat-lon snapshots of column dust opacity, with a temporal resolution of 6 L<sub>s</sub>. After being lifted into the atmosphere, the dust is spread throughout the PBL and advected by the dynamical core.  
*    Turbulence is calculated with a Mellor-Yamada level 2.0 turbulence closure scheme. Heat, momentum, and tracers are mixed in the boundary layer, and the top of the boundary layer is saved for calculating injection of dust and water vapor.  
*    Radiatively inert clouds are formed assuming supersaturated water vapor is immediately condensed into solid ice. Ice is assumed to have a fixed number scaled by pressure.  
*    The radiative heating uses a 12 band correlated-k radiative transfer scheme with a non-LTE correction to the visible heating.  
*    A convective adjustment is iteratively applied from the surface to the top of the atmosphere to remove atmospheric instabilities.  

## License
This software is released under the [NASA Open Source Agreement Version 1.3](NOSA.pdf).

## Notices:

Copyright © 2023 United States Government as represented by the Administrator of the National Aeronautics and Space Administration.  All Rights Reserved.

## Disclaimers

No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER, CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.  FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE, AND DISTRIBUTES IT "AS IS."

Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES, EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL TERMINATION OF THIS AGREEMENT. 
  
