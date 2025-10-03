file_in=$1
file_out=$2

pls40=$(mktemp)
cdo -L fldmean -zonmean -sellonlatbox,0.0,360.0,-41,-39 $file_in $pls40
cdo -L  ymonsub $pls40  -ymonmean -select,year=1981/2010/1 $pls40 $pls40.anom   
cdo -L  ymondiv $pls40.anom  -ymonstd -select,year=1981/2010/1 $pls40 $pls40.stdanom   


pls65=$(mktemp)
cdo -L fldmean -zonmean -sellonlatbox,0.0,360.0,-66,-64 $file_in $pls65
cdo -L  ymonsub $pls65  -ymonmean -select,year=1981/2010/1 $pls65 $pls65.anom   
cdo -L  ymondiv $pls65.anom  -ymonstd -select,year=1981/2010/1 $pls65 $pls65.stdanom   


cdo -L chname,psl,sam -sub $pls40.stdanom $pls65.stdanom $file_out 

rm $pls40 $pls40.anom $pls40.stdanom $pls65 $pls65.anom $pls65.stdanom