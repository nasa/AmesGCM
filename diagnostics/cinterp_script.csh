#!/bin/csh -f

#PBS -l select=1:ncpus=4:mem=56GB
#PBS -l walltime=1:00:00
#PBS -q ldan


# Set up defaults 
#alias diagnos "/u/mkahre/MCMC/analysis/working/variable_diagnostics.py"

set histDir = ( `pwd` )

set RES = c24
set npe = 0
set INFILE = ( )
set NFILES = 0
set date = 0
set VARS = ( ) 
set NVARS = 0
set conserv = 0
set slow = 0
set schmidt = 0
set stretch = 1.
set tlon = 0.
set tlat = -90.
set lonb = 0.
set lone = 360.
set latb = -90.
set late = 90.

set argv = (`getopt o:h:i:n:d:v:c:p:s:q:w:e:r:t:y:u:j: $*`)

while ("$argv[1]" != "--")
  switch ($argv[1])
    case -h:
      set histDir = $argv[2]; shift argv; breaksw
    case -o:
      set outDir = $argv[2]; shift argv; breaksw 
    case -d:
      set date = $argv[2]; shift argv; breaksw 
    case -n:
      set RES = $argv[2]; shift argv; breaksw 
    case -i:
      set INFILE = ( $INFILE $argv[2] ); @ NFILES ++; shift argv; breaksw 
    case -v:
      set VARS = ( $VARS $argv[2] ); @ NVARS ++; shift argv; breaksw 
    case -c:
      set conserv = $argv[2]; shift argv; breaksw 
    case -p:
      set npe = $argv[2]; shift argv; breaksw
    case -s:
      set slow = $argv[2]; shift argv; breaksw
    case -q:
      set schmidt = $argv[2]; shift argv; breaksw
    case -w:
      set stretch = $argv[2]; shift argv; breaksw
    case -e:
      set tlon = $argv[2]; shift argv; breaksw
    case -r:
      set tlat = $argv[2]; shift argv; breaksw
    case -t:
      set lonb = $argv[2]; shift argv; breaksw
    case -y:
      set lone = $argv[2]; shift argv; breaksw
    case -u:
      set latb = $argv[2]; shift argv; breaksw
    case -j:
      set late = $argv[2]; shift argv; breaksw
  endsw
  shift argv
end
shift argv

if ($schmidt > 0) then
  echo 'schmidt transfer: stretch = ' $stretch
  echo 'tlon tlat = ' $tlon ' ' $tlat
  echo 'lonb lone = ' $lonb ' ' $lone
  echo 'latb late = ' $latb ' ' $late
endif
set workdir = $histDir/goof
if (! -d $workdir) then
 mkdir $workdir
endif

#     If no INFILEs are specified; default to all (ie 0)
if ( ! $NFILES > 0  ) then
  echo 'INFILE not specified:  defaulting to all files in each archive'
  set INFILE = 0
endif
echo ' Processing  INFILE = ' $INFILE

#     If  VARS is not specified; default to all ( ie 0)
if ( ! $NVARS > 0  ) then
  echo 'VARS not specified:  defaulting to all variables  in each file'
  set VARS = 0
endif
echo ' Processing  VARS = ' $VARS

#     List all available history files in directory 

set FILES = ( *.nc.tar     )

#      Get only the lastest file...... unless a specific date is requested  
set high = $#FILES
set HistFiles = ( $FILES[$high] )

#       ... unless a specific date is requested           
if ( $date != 0 ) then 
  set HistFiles =  ( $date.nc.* )
endif 


#SLES11
#set regridpath = /u/mkahre/MCMC/analysis/lib/SLES11
#SLES12
set regridpath = /u/mkahre/MCMC/analysis/lib/fregrid_150422
#set regridpath = /u/rurata/fregrid_150422
set fregrid = $regridpath/fregrid
set fregrid_parallel = $regridpath/fregrid_parallel
set HGRID = $regridpath/make_hgrid
set MOSAIC = $regridpath/make_solo_mosaic

source /usr/share/modules/init/csh
module purge

module load comp-intel/2020.4.304
module load mpi-hpe/mpt
module load python3
#SLES12
module load hdf4/4.2.12
module load hdf5/1.8.18_serial
module load netcdf/4.4.1.1_serial 
module load pkgsrc
module load nco/4.6.7

limit stacksize unlimited

if( $conserv == 0 ) then
      set interp_meth = 'bilinear  --center_y'
else if ( $conserv == 1 ) then
      set interp_meth = 'conserve_order1'
else
      set interp_meth = 'conserve_order2'
endif

if ( $RES == c6 ) then
 set NLON = 24
 set NLAT = 12
endif

if ( $RES == c12 ) then
 set NLON = 48
 set NLAT = 24
endif

if ( $RES == c18 ) then
 set NLON = 72
 set NLAT = 36
endif

if ( $RES == c22 ) then
 set NLON = 90
 set NLAT = 45
endif

if ( $RES == c24 ) then
 set NLON = 96
 set NLAT = 48
endif

if ( $RES == c32 ) then
 set NLON = 120
 set NLAT = 60
endif

if ( $RES == c36 ) then
 set NLON = 144
 set NLAT = 72
endif

if ( $RES == c45 ) then
 set NLON = 144
 set NLAT = 90
endif

if ( $RES == c44 ) then
 set NLON = 144
 set NLAT = 90
endif

if ( $RES == c48 ) then
 set NLON = 180
 set NLAT = 90
endif

if ( $RES == c96 ) then
 set NLON = 360
 set NLAT = 180
endif

if ( $RES == c90 ) then
 set NLON = 288
 set NLAT = 180
endif

if ( $RES == c180 ) then
 set NLON = 576
 set NLAT = 360
endif

if ( $RES == c192 ) then
 set NLON = 720
 set NLAT = 360
endif

if ( $RES == c360 ) then
 set NLON = 1152
 set NLAT = 720
endif

if ( $RES == c384 ) then
 set NLON = 1440
 set NLAT = 720
endif

if ( $RES == c720 ) then
 set NLON = 2304
 set NLAT = 1440
endif

if ( $RES == c768 ) then
 set NLON = 2800
 set NLAT = 1440
endif

echo 'Interpolating to resolution:  ' $RES

set icube = ( `echo $RES | cut -c 2- ` )

 @ nlon2 = 2 * $icube

#echo ${nlon2}
cd $workdir

if ( $icube > 99 ) then
#     @ npe = 4 
endif 


if ( $schmidt > 0 ) then
   echo $NLON $NLAT $stretch
#   @ NLON = $NLON * $stretch
#   @ NLAT = $NLAT * $stretch
   echo 'nlon nlat = ' $NLON $NLAT
   $HGRID  --grid_type gnomonic_ed --nlon $nlon2 --do_schmidt --stretch_factor $stretch  --target_lon $tlon  --target_lat $tlat
else
if ( $schmidt < 0 ) then
   echo 'doing conformal cubic grid'
   $HGRID  --grid_type conformal_cubic_grid --nlon $nlon2 --nratio 2
else
   $HGRID  --grid_type gnomonic_ed --nlon $nlon2   
endif
endif

#   Either copy  or create a mosaic file 

### cp $fregrid_parallel:h/{$RES}_cube_mosaic.nc  cube_mosaic.nc

### cp {$RES}_cube_mosaic.nc  cube_mosaic.nc

$MOSAIC --num_tiles 6 --dir ./  --mosaic_name cube_mosaic

#    Otherwise copy  necessary mosaic and horizontal_grid.tile? files 
#               (using appropriate resolution) 

set input_mosaic = cube_mosaic.nc


#     Begin loop over the history files to be regridded 

foreach infile  ( $HistFiles  )

     set histfile = $histDir/$infile


      set name = `echo $infile | sed -e "s/\.tile1\.nc//"`
      echo $name

   set DATE = $name:r:r
   
   if( $INFILE == 0 ) then
    echo " Getting all files in the archive "
          tar -xf $histfile

   else
      echo "Getting only grid_spec and requested file(s):  "  $INFILE
      	  set ext_args  =  ""
    foreach nfdex  (  $INFILE  )
	    echo $nfdex
	                  set ext_args =   " $ext_args   -e $nfdex "
     end
	  set tar_list = ( `tar -tf $histfile  | grep -e grid_spec  $ext_args ` )
          tar -xf $histfile  $tar_list
   endif 


#     List file types, while omitting  grid_spec and horizontal_grid files 

   set diagFiles = (`/bin/ls -1 *.tile1.nc | grep -v grid_spec | grep -v horizontal_grid.tile`)
   set latlonfiles = ()

   set first_fregrid = 1


foreach File ( $diagFiles )
#             Either extract all available variables in the history file ..........
  if( ! $NVARS > 0 ) then
     set variables = (`ncdump -h $File | grep 'grid_yt, grid_xt' | awk '{print $2}' | cut -d\( -f1`)
  else
#         ........... or extract the requested variables
     set variables = ( $VARS )
  endif

#    determine the number of 1-d variables----- in this case, ones that are a function of phalf  or areo
  set onedvars = `ncdump -h $File | egrep  '(float).*((phalf)|(areo))' | awk '{print $2} ' | cut -d\( -f1`

  set variables = `echo $variables |sed 's/ /,/g'`

  set basename = $File:r:r

#  set data_out = ${basename}_lat_lon      
  set data_out = ${basename}      

  echo "Encountered  ' $#onedvars  ' 1-D variables " in  " $File   $basename  ': " $onedvars 

  echo $File
  echo $variables

  if ( $slow != 0 ) then
    if ( $npe > 1 ) then
         mpiexec -np $npe $fregrid_parallel  --input_mosaic $input_mosaic --input_file $basename --interp_method $interp_meth \
                                 --nlon $NLON --nlat $NLAT --scalar_field $variables  \
                                 --lonBegin $lonb --lonEnd $lone --latBegin $latb --latEnd $late  \
                                 --output_file  $data_out
    else
         $fregrid  --input_mosaic $input_mosaic --input_file $basename --interp_method $interp_meth \
                                 --nlon $NLON --nlat $NLAT --scalar_field $variables  \
                                 --lonBegin $lonb --lonEnd $lone --latBegin $latb --latEnd $late  \
                                 --output_file  $data_out
    endif
  else
    if ( $npe > 1 ) then 
         mpiexec -np $npe $fregrid_parallel  --input_mosaic $input_mosaic --input_file $basename --interp_method $interp_meth \
                                 --remap_file fregrid_remap_file --nlon $NLON --nlat $NLAT --scalar_field $variables  \
                                 --lonBegin $lonb --lonEnd $lone --latBegin $latb --latEnd $late  \
                                 --output_file  $data_out
    else
         $fregrid  --input_mosaic $input_mosaic --input_file $basename --interp_method $interp_meth \
                                 --remap_file fregrid_remap_file --nlon $NLON --nlat $NLAT --scalar_field $variables  \
                                 --lonBegin $lonb --lonEnd $lone --latBegin $latb --latEnd $late  \
                                 --output_file  $data_out
    endif
  endif



    if ($status != 0) then
        echo "Error in execution of $fregrid while working on $basename"
        exit 1
    endif
    set latlonfiles = ($latlonfiles $basename  )


   ncrename -h -d grid_xt,lon -d grid_yt,lat -v grid_xt,lon -v grid_yt,lat $data_out.nc

   ncatted  -a long_name,lon,m,c,"longitude" -a long_name,lat,m,c,"latitude" $data_out.nc

#        Add the 1-d fields to the gridded output
   if(  $#onedvars > 0 ) then
          ncks -O -v `echo $onedvars |sed 's/ /,/g'`    $File   akpk.nc
          ncks -A akpk.nc $data_out.nc
          /bin/rm -f akpk.nc
   endif

      # remove files from working directory
###    rm -f $basename.tile*.nc

#      move resulting gridded file back to archive 
    mv $data_out.nc  $histDir/$DATE.$data_out.nc 
end


end     #    Loop over different history files                  

rm -f *.grid_spec.tile?.nc

rm -rf $workdir

#echo $histDir, $DATE,$data_out
#echo 'Printing out diagnostics'
#diagnos -f $DATE.atmos_average.nc $DATE.atmos_diurn.nc $DATE.atmos_daily.nc $DATE.fixed.nc -dir $histDir

exit 0


