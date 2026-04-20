#!/bin/sh
#Copyright © 2023, United States Government, as represented by the Administrator of the National Aeronautics and Space Administration. All rights reserved. 
#The NASA Ames Mars Global Climate Model Patches is licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE‐2.0.
#Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
#script to apply AM4 patches for Ames Mars model

cd ../../FMS && git reset --hard 71957e9d
git checkout -b mars_branch && git apply --reject --whitespace=fix ../AmesGCM/patches/srcFMS.patch
cd ../GFDL_atmos_cubed_sphere && git reset --hard 29654ea
git checkout -b mars_branch && git apply --reject --whitespace=fix ../AmesGCM/patches/srcGFDL_acs.patch
cd ../atmos_drivers && git reset --hard 03d6c84
git checkout -b mars_branch && git apply --reject --whitespace=fix ../AmesGCM/patches/srcatmos_drivers.patch
cd ../FMScoupler && git checkout tags/2021.02 -b mars_branch
cd ../MOM6 && git reset --hard 2e39f174b
cd ../atmos_phys && git reset --hard 46817a1
git checkout -b mars_branch && git apply --reject --whitespace=fix ../AmesGCM/patches/srcatmos_phys.patch
cd ../land_lad2 && git reset --hard c03c4f6
cd ../ocean_shared && git reset --hard e113821
cd ../AmesGCM
cp -r build_run/bin ../../
cp build_run/compile.archives ../../exec
cp build_run/diag_* ../../exec
cp build_run/fms_* ../../exec

exit 0
