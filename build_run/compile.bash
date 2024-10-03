#!/bin/bash
shopt -s expand_aliases
alias unalias='unalias -a'
set -e
cwd=$(pwd)
#--------------------------------------------------------------------------------------------------------
platform="intel"                                     # A unique identifier for your platform
template="$cwd/mkmf.template.$platform"     # path to template for your platform
mkmf="$cwd/../src/mkmf/bin/mkmf"                        # path to executable mkmf
sourcedir="$cwd/../src"                             # path to directory containing model source code
src_share="$cwd/../src/FMS"
src_fv3="$cwd/../src/GFDL_atmos_cubed_sphere"
src_mars="$cwd/../src/AmesGCM"
listpaths="$cwd/../src/mkmf/bin/list_paths"

export NETCDF="/nasa/netcdf/4.4.1.1_mpt"    
export NETCDFPATH="/nasa/netcdf/4.4.1.1_mpt"    
#export GCC_PATH="/usr/bin"
export MPI=""     
#export INTEL_LICENSE_FILE="/nasa/intel/licenses/*.lic"

##set the pre-compiler defs for the model.
##include -DMARS_GDIAGS for global diagnostics output. compiling will take longer
##include RELEASE for release version compile
##use optim = bridge for Sandy Bridge or Ivy Bridge processors
##use optim = well for Haswell or Broadwell processors
##use optim = lake for Cascade Lake or Sky Lake processors
##use optim = rome for Rome processors
##leave optim blank for generic processors
optim=""
cppDefs="-DMARS_GCM -DMARS_SURFACE -DRELEASE"

#--------------------------------------------------------------------------------------------------------
execdir="$cwd/exec.$platform.mars3.1"  # where code is compiled and executable is created
execname="FMS_MARS.x"
executable="$execdir/$execname"

source "$cwd/loadmods.bash"

export NC_BLKSZ="64K"

#--------------------------------------------------------------------------------------------------------
# setup directory structure
if [ ! -d "$execdir" ]; then
  mkdir -p "$execdir"
fi
cd "$execdir"
#--------------------------------------------------------------------------------------------------------

#rjw -----  added new code below and  added $SLIST to the mkmf line.

SLIST=( $sourcedir/atmos_phys/atmos_param/cg_drag/cg_drag.F90 )

echo '----------------*********--------------' 
echo ' Additional Source:  '  "$SLIST"
echo '----------------*********--------------' 


#---------------- write main Makefile

sed -e 's/<TAB>/\t/' >"$execdir/Makefile" <<END

SRCROOT = $sourcedir/
BUILDROOT = $execdir/

MK_TEMPLATE = $template
include \$(MK_TEMPLATE)

${execname}: mars_atmos_phys/libmars_atmos_phys.a atmos_cubed/lib_atmos_cubed.a fms/libfms.a
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@

fms/libfms.a:  FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=fms \$(@F)

atmos_cubed/lib_atmos_cubed.a: fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos_cubed \$(@F)

mars_atmos_phys/libmars_atmos_phys.a: atmos_cubed/lib_atmos_cubed.a fms/libfms.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=mars_atmos_phys \$(@F)


FORCE:

.DEFAULT:
<TAB>-echo \$@ does not exist.
all: a.out
SRC =
OBJ =
clean: neat
<TAB>-rm -f .a.out.cppdefs \$(OBJ) a.out
<TAB>-rm -r fms
<TAB>-rm -r mars_atmos_phys
neat:
<TAB>-rm -f \$(TMPFILES)
TAGS: \$(SRC)
<TAB>etags \$(SRC)
tags: \$(SRC)
<TAB>ctags \$(SRC)
a.out: \$(OBJ)
<TAB>\$(LD) \$(OBJ) -o a.out  \$(LDFLAGS)

END

echo "made Makefile"

if [[ "$cppDefs" == *MARS_GDIAGS* ]]; then

echo "found gdiags"
${listpaths} -o "$execdir/path_names" "$src_fv3/driver/mars/atmosphere.F90" "$sourcedir/atmos_drivers/mars/atmos_model.F90" "${SLIST[@]}" "$src_share" "$src_fv3" "$src_mars"

/bin/mv path_names pathnames.all
egrep -v "*cmip*" pathnames.all > path_names
/bin/mv path_names pathnames.all
egrep -v "*fv_cmp*" pathnames.all > path_names
/bin/mv path_names pathnames.all
egrep -v "atmos_cubed_sphere/driver/GFDL/*"  pathnames.all  > path_names
/bin/mv path_names pathnames.all
egrep -v "atmos_cubed_sphere/driver/SHiELD/*"  pathnames.all  > path_names
mv path_names pathnames.all
egrep -v "*/monin_obukhov/monin_obukhov*"  pathnames.all  > path_names

${mkmf} -a "$sourcedir" -t "$template" -p "${executable##*/}" -c "$cppDefs" "$execdir/path_names" "$src_share/include" "$src_share/mpp/include" /usr/local/include

else

echo "no gdiags"

mkdir -p "$execdir/fms"
cd "$execdir"

/bin/rm -f "$execdir/fms/*.html"
"$listpaths" -o "$execdir/fms/pathnames_fms" "$src_share"

pushd fms

${mkmf} -m Makefile -a "$sourcedir" -p libfms.a -t "$template" -c "$cppDefs" "$execdir/fms/pathnames_fms" "$src_share/include" "$src_share/mpp/include"
popd


mkdir -p "$execdir/atmos_cubed"

/bin/rm -f "$execdir/atmos_cubed/*.html"
${listpaths} -o "$execdir/atmos_cubed/pathnames_atmos_cubed" "$src_fv3"

pushd atmos_cubed


/bin/mv pathnames_atmos_cubed pathnames.all
egrep -v "atmos_cubed_sphere/driver/GFDL/*"  pathnames.all  > pathnames_atmos_cubed
/bin/mv pathnames_atmos_cubed pathnames.all
egrep -v "atmos_cubed_sphere/driver/SHiELD/*"  pathnames.all  > pathnames_atmos_cubed
/bin/mv pathnames_atmos_cubed pathnames.all
egrep -v "*cmip*" pathnames.all > pathnames_atmos_cubed
/bin/mv pathnames_atmos_cubed pathnames.all
egrep -v "*fv_cmp*" pathnames.all > pathnames_atmos_cubed

${mkmf} -m Makefile -a "$sourcedir" -p lib_atmos_cubed.a -t "$template" -c "$cppDefs" -o "-I$execdir/fms" "$execdir/atmos_cubed/pathnames_atmos_cubed" "$src_share/include" "$src_share/mpp/include"
popd


mkdir -p "$execdir/mars_atmos_phys"

/bin/rm -f "$execdir/mars_atmos_phys/*.html"
${listpaths} -o "$execdir/mars_atmos_phys/pathnames_mars_atmos_phys" "$src_fv3/driver/mars/atmosphere.F90" "$sourcedir/atmos_drivers/mars/atmos_model.F90" "$sourcedir/AmesGCM/atmos_cubed_sphere_mars/mars_physics_update.F90" "${SLIST[@]}" "$sourcedir/AmesGCM/atmos_param_mars"


pushd mars_atmos_phys

/bin/mv pathnames_mars_atmos_phys pathnames.all
egrep -v "atmos_cubed_sphere/driver/solo/hswf*"  pathnames.all  > pathnames_mars_atmos_phys

${mkmf} -m Makefile -a "$sourcedir" -p libmars_atmos_phys.a -t "$template" -c "$cppDefs" -o "-I$execdir/fms -I$execdir/atmos_cubed" "$execdir/mars_atmos_phys/pathnames_mars_atmos_phys" "$src_share/include" "$src_share/mpp/include"
popd




# ---------------- adjust the main Makefile

cat Makefile | sed -e 's/<TAB>/\t/' > Makefile.$$ && mv -f Makefile.$$ Makefile


fi

status=$?
if [[ ${status} != 0 ]]; then
   unset echo
   echo "ERROR: mkmf failed for FMS_MARS"
   exit 1
fi

# --- execute make ---
make "${executable##*/}"
status=$?
if [[ ${status} != 0 ]]; then
   unset echo
   echo "ERROR: make failed for FMS_MARS"
   exit 1
fi

unset echo
echo "NOTE: make successful for FMS_MARS"

