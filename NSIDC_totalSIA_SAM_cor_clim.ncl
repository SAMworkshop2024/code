;load "~/scripts/NCL/ncl_funcs.ncl"

;correlat total SIA with SAm by month, and plot

begin

  yrst = 1979      ;analysis start time
  yren = 2023      ;analysis end time
  clst = 1981     ;climatology start date
  clen   = 2010   ;climatiology end date

  diri = "~/data/NSIDC/"
  fili = "NSIDC_SH_totalSIA_monthly_v04r00.nc"
  iNam = "SIE_cdr"

  dira = "~/data/indices/"
  fila = "SAM_GW_1m_1979-2023.nc"
  aNam = "SAM"

  diro = "~/plots/ANTARCTIC/NSIDC/"
  filo = "NSIDC_"+iNam+"_SAM_cor_clim"
  oType = "png"


;*******************
;read data
;********************


;sea ice
  fi   = addfiles(diri+fili, "r")
  YYYYMM = fi[:]->YYYYMM
  st     = ind(YYYYMM.eq.yrst*100+1)
  en     = ind(YYYYMM.eq.yren*100+12)
  SIA    = fi[:]->$iNam$(st:en)
  SIA&time := YYYYMM(st:en)

 
;SAM
  fa  = addfile(dira+fila,"r")
  YYYYMM := cd_calendar(fa->time, -1)
  st  = ind(YYYYMM.eq.(yrst*100+1))
  en  = ind(YYYYMM.eq.(yren*100+12))
  SAM = fa->$aNam$(st:en)
  SAM&time := YYYYMM(st:en)

  SAM@_FillValue := -999.

  replace_ieeenan(SAM, SAM@_FillValue, 0)

;*********************
;detrend for each month
;**********************

  procedure dtren_by_month(x[*]:numeric)
  local tdum
  begin

    cdum = ispan(0, clen-clst, 1)
    tdum = ispan(0,yren-yrst,1)


    do m = 1, 12
      tren     = regCoef(cdum, x({clst*100+m:clen*100+m:12}))
      tren    := tren * tdum
      x((m-1)::12) = x((m-1)::12) - tren
    end do

  end

  dtren_by_month(SAM)
  dtren_by_month(SIA)


;**********************
;correlate
;**********************

  r = new(12,typeof(SAM))
  
  do m = 0, 11
    r(m) = escorc(SAM(m::12), SIA(m::12))
  end do


  N     = yren-yrst
  rdum  = fspan(0.1,0.6,11)
  p     = 1.-rtest(rdum,N, 0)
  rcrit = rdum(min(ind(p.ge.0.9)))
  
  delete([/rdum,p/])


;*********************
;plot correlations
;*********************

  wks = gsn_open_wks(oType, diro+filo)

  res                   = True
  res@gsnFrame          = True
  res@gsnDraw           = True
  res@gsnMaximize       = True
  res@vpHeightF         = 0.5
  res@vpWidthF          = 1.

  res@tiYAxisString     = "r"
  res@gsnYRefLine       = (/0.,-rcrit, rcrit/)
  res@gsnYRefLineDashPatterns = (/0, 1, 1/)
  res@xyLineThicknessF  = 3.5
  res@trYMaxF           = 0.5
  res@trYMinF           = - res@trYMaxF

  res@tmXBMode          = "Explicit"
  res@trXMaxF           = 12
  res@trXMinF           = 1
  res@tmXBValues        = ispan(2, res@trXMaxF,2)
  res@tmXBLabels        = (/"FEB","APR","JUN","AUG","OCT","DEC"/)
  res@tmXBMinorValues   = ispan(res@trXMinF, res@trXMaxF,1)



  plot = gsn_csm_xy(wks,res@tmXBMinorValues, r, res)



end
