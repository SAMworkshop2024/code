begin

 Corr = new((/4,12/),float)

 samf = addfile("/scratch/v45/SAMworkshop2024/data/SAM_GW_1m_1979-2023.nc","r")

 msam = doubletofloat(samf->SAM(:527))
 msam@_FillValue = -999
 printVarSummary(msam)

 ssam = runave_Wrap(msam,3,0)
 ssam@_FillValue = -999
 ptand = onedtond(ssam,(/44,12/))

 p24abcf = asciiread("/scratch/v45/SAMworkshop2024/data/access-s2.9memM.SAManom.1982-2018.txt",(/37*12,6/),"float")

 yr = ispan(1982,2018,1)
 nyr = dimsizes(yr)

 obs = ptand(3:39,0:)    ; reading from DJF

 odata = new((/nyr,12/),float)
 odata(:,0:10) = obs(:,1:11)
 odata(0:nyr-2,11) = obs(1:nyr-1,0)
 odata(nyr-1,11) = -999

 p24abc = p24abcf(:,2:)

 pdata = new((/12,nyr,4/),float)
 pdata@_FillValue = -999

 do i=0,11
  poama = p24abc(i::12,:)
  printVarSummary(poama)
  pdata(i,:,:) = poama
 end do

 odata!0 = "year"
 odata!1 = "month"
 odata@_FillValue = -999

 pdata!0 = "month"
 pdata!1 = "year"
 pdata!2 = "leadT"

;  asciiwrite("pdata.txt",pdata(11,:,0))
;  print(pdata(leadT|2,month|:,year|22))

  do i=0,3                                    ; lead time
   do j=0,11
       vm = j + i
       sy = 0
       ey = nyr-1
       if(vm.ge.12) then
          vm = vm - 12
          sy = 1
          ey = nyr-2 
       end if
       
       Corr(i,j) = escorc(pdata(leadT|i,month|j,year|:ey),odata(month|vm,year|sy:))
       print(pdata(leadT|i,month|j,year|:ey)+"  "+odata(month|vm,year|sy:))
   end do
  end do

  Corr!0 = "leadT"
  Corr!1 = "FStartM"

  Corr&leadT = fspan(0,3.0,4)
  Corr&FStartM = fspan(1,12.0,12)
  printVarSummary(Corr)

  print(Corr(0,:))  ; Forecast skill at lead time 0

end
