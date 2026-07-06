;Wrote by JuanLi, IPRC,UH
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"    
;load "./vertical_thete.ncl"

  
;******************** MAIN **********************************

begin

 nameRegion = "EIO"        ; Indian Ocean base region
  latS_IO    = -10.
  latN_IO    =  10.
  lonL_IO    =  80.
  lonR_IO    = 100.

  latn       =  5.        ; lat band for (lag,lon) Fig 5
  lats       = -5.

  pltName    = "mjo_theta_850_400"       ; output plot name
  pltType    = "png"       ; x11, ps, eps, pdf, png
  pltDir     = getenv("FIG_DIR")        ; output plot directory
  ncDir      = getenv("NC_DIR")
  modout     = getenv("modout")

  dir   = "./OBS/"
  pName  = "PR"
  filp   = "GPCP_1DD_v1.2_199610-201312.nc"
  rName  = "Q"
  filr   = "ERA_q_1997-2013.nc"
  aName  = "T"
  fila   = "ERA_t_1997-2013.nc"

  twStrt = 19790101        ; Period for OBS
  twLast = 19951231        ;

  cwStrtm = stringtoint(getenv("ymdStrt"))        ; Period for MODEL
  cwLastm = stringtoint(getenv("ymdLast"))

  dirm   = "../modeldata/NUIST/test/"
  film   = "day.nc"
  pNamem  = "PRECT"
  aNamem  = "T"
  rNamem  = "Q"

  nameSeason = (/ "Winter"/)
;************************************************
; create Lanczos BandPass Filter
;************************************************
  ihp    = 2                          ; bpf=>band pass filter
  nWgt   = 141
  sigma  = 1.0                        ; Lanczos sigma
  fca    = 1./70.                    ; MJO clivar
  fcb    = 1./20.
  wgt    = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
;***************PRE************************************
;OBS*****
  f      = addfile(getenv("obs_data"), "r")
  date_p = cd_calendar(f->time, -2)   ; entire file
  iStrt  = ind(date_p.eq.twStrt)      ; desired dates
  iLast  = ind(date_p.eq.twLast)
  delete(date_p)
                                      ; P(time,lat,lon)
  P  = f->$pName$(iStrt:iLast,{-30:30},:)
  printMinMax(P, True)

  wyIO   = f->lat({latS_IO:latN_IO})
  wyIO   = cos(0.017459*wyIO)         ; spatial wgt
;MODEL****
  f1      = addfile(getenv("mod_data_daily"),"r")
  date_p = cd_calendar(f1->time, -2)   ; entire file
  iStrt1  = ind(date_p.eq.cwStrtm)      ; desired dates
  iLast1  = ind(date_p.eq.cwLastm)
  iStrt1 = 1
  iLast1 = 4379
  delete(date_p)
                                      ; P(time,lat,lon)
  Pm  = f1->$pNamem$(iStrt1:iLast1,:,:)
  Pm  =Pm*86400000
  printMinMax(Pm, True)

  wyIOm   = f1->lat({latS_IO:latN_IO})
  wyIOm   = cos(0.017459*wyIOm)
;************************************************
; Create wgted area average of the base IO precip series  (time)
; Really, no need to area weight here .... area is very small.
;************************************************
  PIO    = wgt_areaave_Wrap(P(:,{latS_IO:latN_IO},{lonL_IO:lonR_IO}), wyIO, 1., 0)
 ; PIO    = dtrend (PIO, False)                                   ; rmv overall trend
  PIO    = wgt_runave_leftdim( PIO, wgt, 0 )                     ; apply filter

  PIOm    = wgt_areaave_Wrap(Pm(:,{latS_IO:latN_IO},{lonL_IO:lonR_IO}), wyIOm, 1., 0)
 ; PIOm    = dtrend (PIOm, False)                                   ; rmv overall trend
  PIOm    = wgt_runave_leftdim( PIOm, wgt, 0 )
delete([/P,Pm/])
;*********************************************************************
;***********************theta_eqv*************************************
;OBS*******
  f      =addfile(getenv("obs_data_t"), "r")
  date_a = cd_calendar(f->time, -2)     ; entire file
  iStrt  = ind(date_a.eq.twStrt)        ; desired dates
  iLast  = ind(date_a.eq.twLast)
  
  print(f)
  print(iStrt)
  print(iLast)
  delete(date_a)

;  p=(/1000,850,700,500,250/)
  p=850

  f      = addfile(getenv("obs_data_q"), "r")
  print(f)
  SH  = short2flt( f->$rName$(iStrt:iLast,{p},{-30:30},:))
  f      = addfile(getenv("obs_data_t"), "r")
  AIR  = short2flt( f->$aName$(iStrt:iLast,{p},{-30:30},:))

  p2=400
  f      = addfile(getenv("obs_data_q"), "r")
  print(f)
  SH2  = short2flt( f->$rName$(iStrt:iLast,{p2},{-30:30},:))
  f      = addfile(getenv("obs_data_t"), "r")
  AIR2  = short2flt( f->$aName$(iStrt:iLast,{p2},{-30:30},:))

  nlev=dimsizes(p)
 ; X = mixhum_ptrh (conform(AIR,p,1), AIR, RH, 1); convert rel humidity to mixing ratio
  X = SH/(1-SH) ; convert spc humidity to mixing ratio
  p0 = 1000
  P = 850
  p0_at_units = "hPa" ; must match units of p
;*****cal theta************
   L_v = 2510000 ; latent heat of evaporation [approx 2400 kJ/kg at 25C]
   c_pd = 1004. ; specific heat at constant pressure for air [approx 1004 J/(kg-K)]
   R = 287. ; specific gas constant for air [J/(kg-K)]
   Rcpd = R/c_pd
   T_e = AIR + (L_v/c_pd)*X ; common approximation
  THETAE  = T_e*(p0/P)^Rcpd
   copy_VarCoords(AIR,THETAE)
   printMinMax(THETAE,True)

;*****cal theta************
 X2 = SH2/(1-SH2) ; convert spc humidity to mixing ratio
  p0 = 1000
  P2 = 400
   T_e2 = AIR2 + (L_v/c_pd)*X2 ; common approximation
  THETAE2  = T_e2*(p0/P2)^Rcpd
   copy_VarCoords(AIR2,THETAE2)
   printMinMax(THETAE2,True)
  THETAE = THETAE - THETAE2


;*********MODEL*******
 p = 850
 AIRm  = f1->$aNamem$(iStrt1:iLast1,{p},{-30:30},:)
 printMinMax(AIRm,True)
 SHm  = f1->$rNamem$(iStrt1:iLast1,{p},{-30:30},:)
 printMinMax(SHm,True)
 Xm = SHm/(1-SHm)
 delete(SHm)

 p2 = 400
 AIRm2  = f1->$aNamem$(iStrt1:iLast1,{p2},{-30:30},:)
 printMinMax(AIRm,True)
 SHm2  = f1->$rNamem$(iStrt1:iLast1,{p2},{-30:30},:)
 printMinMax(SHm2,True)
 Xm2 = SHm2/(1-SHm2)
 delete(SHm2)

 P = 850
 p0 = 1000
 p0_at_units = "Pa" ; must match units of p
;*****cal theta************
   L_v = 2510000 ; latent heat of evaporation [approx 2400 kJ/kg at 25C]
   c_pd = 1004. ; specific heat at constant pressure for air [approx 1004 J/(kg-K)]
   R = 287. ; specific gas constant for air [J/(kg-K)]
   Rcpd = R/c_pd
   T_em = AIRm + (L_v/c_pd)*Xm ; common approximation
  THETAEm  = T_em*(p0/P)^Rcpd
     copy_VarCoords(AIRm,THETAEm) 
   printMinMax(THETAEm,True)

   P2 = 400
   T_em2 = AIRm2 + (L_v/c_pd)*Xm2 ; common approximation
  THETAEm2  = T_em2*(p0/P2)^Rcpd
     copy_VarCoords(AIRm2,THETAEm2)
   printMinMax(THETAEm2,True)

    THETAEm = THETAEm - THETAEm2

  date=cd_calendar(THETAE&time,-2)
  datem=cd_calendar(THETAEm&time,-2)
;**********************OBS************************
  nDay = 180
  mmddStrt = 1101
  yyyy  = date/10000
  mmdd  = date - yyyy*10000
  nSea  = num(mmdd.eq.mmddStrt)-1           ; number of seasons
  iSea  = ind(mmdd.eq.mmddStrt)

  pio=new((/nDay*nSea/),"float")
  q_time=new((/nDay*nSea,dimsizes(THETAE&lat),dimsizes(THETAE&lon)/),"float")
  u_time=q_time
  v_time=q_time
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSea(ns)
    iLast  = iSea(ns)+nDay-1
    pio(ns*180:180*(ns+1)-1)=tofloat(PIO(iStrt:iLast))
    q_time(ns*180:180*(ns+1)-1,:,:)=THETAE(iStrt:iLast,:,:)
  end do

pio=dim_standardize(pio,1)
regq=regCoef_n(pio,q_time, 0, 0)*3
copy_VarMeta(THETAE(0,:,:),regq)
;**********************MODEL************************

  yyyym  = datem/10000
  mmddm  = datem - yyyym*10000
  nSea  = num(mmddm.eq.mmddStrt)-1           ; number of seasons
  iSeam  = ind(mmddm.eq.mmddStrt)

  piom=new((/nDay*nSea/),"float")
  q_timem=new((/nDay*nSea,dimsizes(THETAEm&lat),dimsizes(THETAEm&lon)/),"float")
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSeam(ns)
    iLast  = iSeam(ns)+nDay-1
    piom(ns*180:180*(ns+1)-1)=tofloat(PIOm(iStrt:iLast))
    q_timem(ns*180:180*(ns+1)-1,:,:)=THETAEm(iStrt:iLast,:,:)
  end do
piom=dim_standardize(piom,1)
regqm=regCoef_n(piom,q_timem, 0, 0)*3
copy_VarMeta(THETAEm(0,:,:),regqm)
 maxf = max (abs(regqm))
    regqm = regqm /maxf*2.

 maxf = max (abs(regq))
    regq = regq /maxf*2.

;*********************************************************************
  wks = gsn_open_wks(pltType,pltDir+"/"+modout+"."+pltName)
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@gsnAddCyclic = False
res@tmXBLabelFontHeightF = 0.012
res@tmYLLabelFontHeightF = 0.012
res@tmYLLabelFontThicknessF = 0.6
res@tmXBLabelFontThicknessF = 0.6
res@tmYLLabelStride = 0
res@tmXBLabelFont = 25
res@tmYLLabelFont = 25
res@tmYLMode = "Explicit"
res@tmXBMode = "Explicit"
res@tmXBValues = fspan(40,180,8)
res@tmYLValues = fspan(-30,30,5)
res@tmYLLabels = (/"30~S~o~N~S","15~S~o~N~S","EQ","15~S~o~N~N","30~S~o~N~N"/)
res@tmXBLabels = (/"40~S~o~N~S","60~S~o~N~E","80~S~o~N~E","100~S~o~N~E","120~S~o~N~E","140~S~o~N~E","160~S~o~N~E","180~S~o~N"/)
res@tmYROn     = False
res@tmXTOn     = False
res@cnFillOn   = True
res@mpMinLatF            = -25.         ; zoom in on map
res@mpMaxLatF            = 25.
res@mpMinLonF            = 40.
res@mpMaxLonF            = 190.
res@mpCenterLonF         = 115.
res@cnInfoLabelOn        = False
res@cnLinesOn  = False
res@cnLineLabelsOn = False
res@gsnStringFontHeightF = 0.014
res@gsnLeftString = "OBS "
res@gsnRightString = ""
res@gsnSpreadColors  = False
res@cnLevelSelectionMode ="ExplicitLevels"          ;"AutomaticLevels"     ; "ManualLevels"
res@cnLevels = (/-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6/)
res@cnLevels = (/-1.2,-0.9,-0.6,-0.3,0.3,0.6,0.9,1.2/)
res@cnLevels = (/-1.0,-0.7,-0.4,-0.1,0.1,0.4,0.7,1.0/)

res@cnFillColors=(/3,18,38,58,0,170,190,210,230/)
res@lbLabelBarOn= False
res@pmLabelBarOrthogonalPosF= 0.12
res@lbLabelFontHeightF = 0.011
res@pmLabelBarWidthF = 0.5
res@pmLabelBarHeightF = 0.03
res@mpGeophysicalLineColor = "gray30"
res@lbOrientation        = "Vertical"

  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02
; change label spacing to avoid overlap
  res@lbLabelStride    = 0.1                        ; every other label
  res@tmXBTickSpacingF = 10.                      ; every 60 degree
  res@vpHeightF= 0.4                    ; change aspect ratio of plot
  res@vpWidthF = 0.8
  res@lbLabelStride = 2
  res@pmLabelBarWidthF = 0.7
  res@pmLabelBarHeightF = 0.07
  res@lbLabelFontHeightF = 0.018
  res@lbLabelFont = 25
  res@pmLabelBarOrthogonalPosF= 0.05
  res@tmYROn               = False    ; Turn off right tickmarks
  res@tmXTOn               = False   ; Turn off top tickmarks
  res@tiYAxisString        = "lag (days)"





res2=True
res2@gsnDraw = False
res2@gsnFrame = False
res2@gsnAddCyclic = False
res2@vcRefAnnoOrthogonalPosF = -0.285
res2@vcRefAnnoParallelPosF   = 1.0
res2@vcRefMagnitudeF         = 1.0
res2@vcRefLengthF            = 0.03
res2@vcGlyphStyle            = "CurlyVector"
res2@vcLineArrowThicknessF   = 1.2
res2@vcMinDistanceF          = 0.018
res2@vcRefAnnoFontHeightF    = 0.007
res2@vcRefAnnoString1        = "1.0"
res2@vcRefAnnoString2On      = False
res2@vcRefAnnoFont           = 25
res2@gsnLeftString           = ""
res2@gsnRightString           = ""

res1=True
res1@gsnSpreadColors  = True
res1@gsnSpreadColorStart = 1
res1@gsnSpreadColorEnd   = 1
res1@gsnDraw = False
res1@gsnFrame = False
res1@cnFillOn   = False
res1@cnLinesOn  = True
res1@cnLineLabelsOn = False
res1@cnLineLabelPlacementMode="Computed"   ;"Constant"
res1@cnLineLabelDensityF=1.
res1@cnLineLabelInterval=1
res1@cnLineLabelFontHeightF=0.007
res1@cnLineLabelFontThicknessF=2.
res1@cnLabelScalingMode="TrimZeros"
res1@cnLineLabelFormat = "0@;*.1f"
res1@cnLineLabelFont=22;"times-roman"
res1@cnMonoLineLabelFontColor= True
res1@cnLineLabelFontColors = 1
res1@cnInfoLabelOn  = False
res1@cnMonoLineColor = False
res1@cnLineColors = (/"blue","blue","blue","blue","blue","blue","blue","blue"/)
res1@cnLevelSelectionMode = "ExplicitLevels"
res1@cnLevels = (/-0.2,-0.15,-0.1,-0.05,0.05,0.1,0.15,0.2/)
res1@gsnContourNegLineDashPattern = 2 
res1@cnLineThicknesses = True
res1@cnLineThicknessF = 2.
res1@cnLineLabelBackgroundColor=0
res1@gsnLeftString     = ""
res1@gsnRightString    =""

   ;a1=ndtooned(regqm({-10:10},{80:150}))
   ;a2=ndtooned(regq({-10:10},{80:150}))
   ;cor=escorc(a1,a2)
   ;print(cor)
   ;rmse=dim_rmsd(a1,a2)
   ;y2=dim_stddev(a2)
   ;nrmse=rmse/y2
   ;print(nrmse)

;   a1=ndtooned(regdvm({-10:10},{100:140}))
;   a2=ndtooned(regdv({-10:10},{100:140}))
   ;cor1=escorc(a1,a2)
   ;print(cor1)
   ;rmse=dim_rmsd(a1,a2)
   ;y2=dim_stddev(a2)
   ;nrmse1=rmse/y2
   ;print(nrmse1)

plot = new(2,graphic)
plot(0) = gsn_csm_contour_map(wks,regq,res)
;plot1=gsn_csm_vector(wks,regu,regv,res2)
;overlay(plot(0),plot1)
;plot2= gsn_csm_contour(wks,regdv,res1)
;overlay(plot(0),plot2)
res@gsnLeftString           = "MODEL"
;res@gsnRightString    = "Q250:PCC="+sprintf("%5.2f",cor)+" NRMSE="+sprintf("%5.2f",nrmse)+"   Div250:PCC="+sprintf("%5.2f",cor1)+" NRMSE="+sprintf("%5.2f",nrmse1)
plot(1) = gsn_csm_contour_map(wks,regqm,res)
;plot1=gsn_csm_vector(wks,regum,regvm,res2)
;overlay(plot(1),plot1)
;plot2= gsn_csm_contour(wks,regdvm,res1)
;overlay(plot(1),plot2)

dum1=new(4,graphic)
dum2=new(4,graphic)
ypts1 = (/ -10.0,  -10.0,  10.0, 10.0, -10.0/) ;mega naf,saf
xpts1 = (/ 150.0, 80., 80., 150.0, 150.0/)
gres                  = True                ; polyline mods desired
gres@gsLineThicknessF = 3.5
gres@gsLineColor      = "black"

do i = 0 , 3
  dum1(i)=gsn_add_polyline(wks,plot(0),xpts1(i:i+1),ypts1(i:i+1),gres)
  dum2(i)=gsn_add_polyline(wks,plot(1),xpts1(i:i+1),ypts1(i:i+1),gres)
end do

  resP=True
  resP@gsnPanelLabelBar     = False
  resP@gsnPanelYWhiteSpacePercent = 3
  resP@txFont =25
  resP@txFontHeightF =0.018
  resP@txString            =  "Regressed convective instability (EPT850-EPT400,K)"
  gsn_panel(wks,plot,(/2,1/),resP)
;----------------------------------

;------------------output------------------------

dirName = "./NC_OUT/"
; Check if the directory exists and create it if it doesn't
dirCheck = systemfunc("test -d " + dirName + " ; echo $?")

if (dirCheck .eq. "0") then
    print("Directory " + dirName + " already exists.")
else
    system("mkdir -p " + dirName)
    print("Directory " + dirName + " created.")
end if

filename = "./NC_OUT/" + modout + ".mjo_conindex_obs.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->regq=regq
fout->pio=pio
fout->q_time=q_time

filename = "./NC_OUT/" + modout + ".mjo_conindex_model.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->regqm=regqm
fout->piom=piom
fout->q_timem=q_timem

;-----------------------------------------------






end 
