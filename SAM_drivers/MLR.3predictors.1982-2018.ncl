
begin
 
 samf = addfile("/scratch/v45/SAMworkshop2024/data/SAM_GW_1m_1979-2023.nc","r")
 msam = doubletofloat(samf->SAM(:527))
 msam@_FillValue = -999
 printVarSummary(msam)

 ssam = runave_Wrap(msam,3,0)
 ssam@_FillValue = -999
 ;print(ssam(:11))
 ssam_sub = ssam(36:479)

 epENSOf = addfile("/scratch/v45/SAMworkshop2024/data/msked.hurrell-reynolds.monthly.nino3.upto2022.clim1981-2010.nc","r")
 nino = epENSOf->nino3
 printVarSummary(nino)

 snino = runave_Wrap(nino,3,0)
; print(snino(:11))
 snino@_FillValue = -999
; xtime = ut_calendar(snino&time,-2)
; print(xtime+"  "+snino)
 snino_sub = snino(36:479)

 cpENSOf = addfile("/scratch/v45/SAMworkshop2024/data/msked.hurrell-reynolds.monthly.emi.upto2022.clim1981-2010.nc","r")
 emi = cpENSOf->emi
 printVarSummary(emi)
 semi = runave_Wrap(emi,3,0)
; print(semi(:11))
 semi@_FillValue = -999
 semi_sub = semi(36:479)

 spvf = addfile("/scratch/v45/SAMworkshop2024/data/zmUat10hPa_60S.1979-2023.nc","r")
 spv = spvf->u10anom(:527)
 printVarSummary(spv)

 sspv = runave_Wrap(spv,3,0)
; print(sspv(:11))
 sspv@_FillValue = -999
 sspv_sub = sspv(36:479)

;;-----------------------------------------
 
yr = ispan(1982,2018,1)
nyr = dimsizes(yr)

cc = new(12,float)

ptor = new((/3,nyr,12/),float)

ptor(0,:,:) = onedtond(sspv_sub,(/nyr,12/))
ptor(1,:,:) = onedtond(snino_sub,(/nyr,12/))
ptor(2,:,:) = onedtond(semi_sub,(/nyr,12/))

ptand = onedtond(ssam_sub,(/nyr,12/))

ptor!0 = "predictors"
ptor!1 = "year"
ptor&year = ispan(1982,2018,1)
ptor!2 = "month"
ptor&month = ispan(1,12,1)

ptand!0 = "year"
ptand&year = ispan(1982,2018,1)
ptand!1 = "month"
ptand&month = ispan(1,12,1)

sim_cc = new((/3,12/),float) ; simultaneous correlation between drivers and sam

do m=0,11                    ; m=0 DJF m=11 NDJ 
 do pp = 0,2
   sim_cc(pp,m) = escorc(ptor(pp,:,m),ptand(:,m))
 end do
end do

print(sim_cc(0,:)+"  "+sim_cc(1,:)+"  "+sim_cc(2,:))
exit

;-------------------------------------------------------

y_hat = ptand   
y_hat = -999
y_hat@_FillValue = -999

do i=0,nyr-1   ; Predicting SAM with cross-validation

 idx = ind(yr.ne.yr(i))

 do m=0,11
;
   p1 = ptor(0,idx,m)
   p2 = ptor(1,idx,m)
   p3 = ptor(2,idx,m)
; 
   q = ptand(idx,m)
;
   y = q

   x = new((/4,nyr-1/),float)
   x(0,:) = 1.0
   x(1,:) = p1 
   x(2,:) = p2  
   x(3,:) = p3
;
   beta = reg_multlin(y,x,False)
;
   y_hat(i,m)  = beta(0) + beta(1)*ptor(0,i,m)+ beta(2)*ptor(1,i,m) + beta(3)*ptor(2,i,m) 
;
   delete([/p1,p2,p3,q,x,y,beta/])
 end do
end do

ptand@_FillValue = -999
y_hat@_FillValue = -999

do m=0,11
 cc(m) = escorc(ptand(:,m),y_hat(:,m))  
end do
 print(cc)

end
