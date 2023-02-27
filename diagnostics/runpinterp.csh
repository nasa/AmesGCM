#!/bin/csh -f

source /usr/share/modules/init/csh
module purge

module load comp-intel
module load mpi-sgi/mpt
#SLES12
module load hdf4/4.2.12
module load hdf5/1.8.18_mpt
module load netcdf/4.4.1.1_mpt
module load pkgsrc/2018Q3
module load nco/4.6.7

setenv PATH /u/$USER/FRE-NCtools/bin:$PATH
set PINTERP = /u/$USER/FRE-NCtools/bin/plevel.sh    #may require absolute path

set date = 0
set ftype = atmos_average
set all = 0
set levels = 0
set sfile = 0

set argv = (`getopt d:f:a:g:s: $*`)

while ("$argv[1]" != "--")
  switch ($argv[1])
    case -d:
      set date = $argv[2]; shift argv; breaksw
    case -f:
      set ftype = $argv[2]; shift argv; breaksw
    case -a:
      set all = $argv[2]; shift argv; breaksw
    case -g:
      set levels = $argv[2]; shift argv; breaksw
    case -s:
      set sfile = $argv[2]; shift argv; breaksw
  endsw
  shift argv
end
shift argv

if ( $levels == 0 ) then
set levels = (  1.e-5 3.e-5 5.e-5 1.e-4 3.e-4 5.e-4 \
                0.003 0.005 0.01 0.03 0.05 0.1 0.2  \
                0.3 0.5 0.01e2 0.02e2 0.03e2 0.05e2 \
                0.07e2 0.1e2 0.2e2 0.3e2 0.5e2 0.7e2 \
                1.0e2 1.5e2 2.0e2 2.5e2 3.0e2 3.5e2 \
                4.0e2 4.5e2 5.0e2 5.3e2 5.5e2 5.9e2 \
                6.0e2 6.3e2 6.5e2 6.9e2 7.e2 7.5e2 \
                8.e2 8.5e2 9.0e2 9.5e2 10.e2 )
else
set levtmp = `cat $levels`
set levels = ( ${levtmp} )
endif

if ( $ftype != atmos_average ) then
   set FILES = ( *.${ftype}.nc )
else
   set FILES = ( *atmos_average.nc )
endif

if ( $date != 0 ) then
  set FILES =  ( $date.$ftype.nc )
endif

if ( $sfile != 0 ) then
  set FILES = ( ${sfile} )
endif


set index =  $FILES[1]

touch xxx

foreach index ( $FILES )

   set INFILE = $index;

   set base = $index:r
   set OUTFILE =    {$base}_pstd.nc
   ncdump -h $INFILE > dump.out
   if ( { grep -q "float ak" dump.out } ) then
       set pk_name = "ak"
   else
       set pk_name = "pk"
   endif
   echo $pk_name
   rm -f dump.out

   echo $OUTFILE
  $PINTERP  -a  -p "$levels" -k $pk_name -i $INFILE -o $OUTFILE  >> xxx
   ncks  -A  -v areo $INFILE   $OUTFILE
   ncks --overwrite --mk_rec_dmn time $OUTFILE $OUTFILE

   

end



