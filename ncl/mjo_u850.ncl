;Wrote by JuanLi, IPRC,UH
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"    
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"    

  
;******************** MAIN **********************************
begin
  cwStrt = 19790101        ; Period for OBS 
  cwLast = 19951230 
 
  dir   = "./OBS/"
  pName  = "PR"
  filp   = "GPCP_1DD_v1.2_199610-201312.nc"
  uName  = "U"
  filu   = "ERA_u_1997-2013.nc"
  vName  = "V"
  filv   = "ERA_v_1997-2013.nc"
  rName  = "Q"
  filr   = "ERA_q_1997-2013.nc"


  cwStrtm = stringtoint(getenv("ymdStrt"))        ; Period for MODEL
  cwLastm = stringtoint(getenv("ymdLast"))

  dirm   = "../modeldata/NUIST/test/"
  film   = "day.nc"
  pNamem  = "PRECT"
  uNamem  = "U"
  vNamem  = "V"
  rNamem  = "Q" 

  nameRegion = "IO"        ; Indian Ocean base region
  latS_IO    =  -10.
  latN_IO    =   10. 
  lonL_IO    =  80.
  lonR_IO    =  100. 

  if (isvar("pm")) then
    print("The parameter value is: " + pm)
  else
    print("No parameter value passed. Using default.")
    pm = 850   ; Set a default value if no parameter is passed
  end if
  p=850  ;850hpa

  nameSeason = (/ "Winter"/)

  latn       =  5.        ; lat band for (lag,lon) Fig 5
  lats       = -5.



  pltName    = "mjo_u850"       ; output plot name         
  pltType    = "png"       ; x11, ps, eps, pdf, png 
  pltDir     = getenv("FIG_DIR")        ; output plot directory
  ncDir      = getenv("NC_DIR")
  modout     = getenv("modout")
;************************************************
; create Lanczos BandPass Filter
;************************************************
  ihp    = 2                          ; bpf=>band pass filter
  nWgt   = 141
  sigma  = 1.0                      ; Lanczos sigma
  fca    = 1./70.                    ; MJO clivar
  fcb    = 1./20.
  wgt    = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )
;************************************************************
; PRECIPITATION:
;     time indices corresponding to the desired time window
;     Read user specified period
;**********************OBS*************************************

  f      = addfile(getenv("obs_data_p"), "r")  
  date_p = cd_calendar(f->time, -2)   ; entire file
  iStrt  = ind(date_p.eq.cwStrt)      ; desired dates
  iLast  = ind(date_p.eq.cwLast)
  delete(date_p)
                                      ; P(time,lat,lon)
  P  = f->$pName$(iStrt:iLast,:,:) 
 ; printVarSummary( P )
;  printMinMax(P, True)
 
  wyIO   = f->lat({latS_IO:latN_IO})  
  wyIO   = cos(0.017459*wyIO)         ; spatial wgt
;*****************MODEL*******************************
  f1      = addfile(getenv("mod_data_daily"),"r")
  date_p = cd_calendar(f1->time, -2)   ; entire file
  iStrt1  = ind(date_p.eq.cwStrtm)      ; desired dates
  iLast1  = ind(date_p.eq.cwLastm)
  iStrt1 = 1
  iLast1 = 4379
  delete(date_p)
                                      ; P(time,lat,lon)
  print(iStrt1)
  print(iLast1)

  Pm  = f1->$pNamem$(iStrt1:iLast1,:,:)
  Pm  =Pm*86400000
;  printMinMax(Pm, True)

  wyIOm   = f1->lat({latS_IO:latN_IO})
  wyIOm   = cos(0.017459*wyIOm)         ; spatial wgt
;************************************************************
; U ANOMALIES:
;     time indices corresponding to the desired time window
;     Read user specified period
;*********************OBS***********************************
  f      = addfile(getenv("obs_data_u"), "r")  
  date_p = cd_calendar(f->time, -2)   ; entire file
  iStrt  = ind(date_p.eq.cwStrt)      ; desired dates
  iLast  = ind(date_p.eq.cwLast)
  delete(date_p)

  U  = short2flt( f->$uName$(iStrt:iLast,{p},:,:))
  date_u=cd_calendar(U&time, -2)
  Um  = f1->$uNamem$(iStrt1:iLast1,{pm},:,:)
  date_um=cd_calendar(Um&time, -2)



  f = addfile(getenv("obs_data_v"), "r")
  V  = short2flt( f->$vName$(iStrt:iLast,{p},:,:))
  Vm  =  f1->$vNamem$(iStrt1:iLast1,{pm},:,:)


;*******rhsm->shsm**********************
 f      = addfile(getenv("obs_data_q"), "r")
 SH  = short2flt( f->$rName$(iStrt:iLast,{p},:,:))
; printMinMax(SH,True)
;
 SHm  = f1->$rNamem$(iStrt1:iLast1,{pm},:,:)
; printMinMax(SHm,True)

;*******U850************************
  u=U
  v=V
 
  um=Um
  vm=Vm
  

  date=cd_calendar(U&time,-2)
  datem=cd_calendar(Um&time,-2)
  ;print("date m = ")
  ;print(datem)

;************************************************
; Create wgted area average of the base IO precip series  (time)
; Really, no need to area weight here .... area is very small.
;************************************************
  PIO    = wgt_areaave_Wrap(P(:,{latS_IO:latN_IO},{lonL_IO:lonR_IO}), wyIO, 1., 0)
;  PIO    = dtrend (PIO, False)                                   ; rmv overall trend
  PIO    = wgt_runave_leftdim( PIO, wgt, 0 )                     ; apply filter
 
  PIOm    = wgt_areaave_Wrap(Pm(:,{latS_IO:latN_IO},{lonL_IO:lonR_IO}), wyIOm, 1., 0)
;  PIOm    = dtrend (PIOm, False)                                   ; rmv overall trend
  PIOm    = wgt_runave_leftdim( PIOm, wgt, 0 ) 
;************************************************
; Create LAT average of the global Precip and U series  (time,lon)
; Really, no need to area weight here
;**********************OBS************************
  nDay = 180
  mmddStrt = 1101
  yyyy  = date/10000
  mmdd  = date - yyyy*10000
  nSea  = num(mmdd.eq.mmddStrt)-1           ; number of seasons
  iSea  = ind(mmdd.eq.mmddStrt)

  pio=new((/nDay*nSea/),"float")
  u_time=new((/nDay*nSea,dimsizes(u&lat),dimsizes(u&lon)/),"float")
  v_time=u_time
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSea(ns)
    iLast  = iSea(ns)+nDay-1
    pio(ns*180:180*(ns+1)-1)=tofloat(PIO(iStrt:iLast))
    u_time(ns*180:180*(ns+1)-1,:,:)=U(iStrt:iLast,:,:)
    v_time(ns*180:180*(ns+1)-1,:,:)=V(iStrt:iLast,:,:)
  end do

pio=dim_standardize(pio,1)
regu=regCoef_n(pio,u_time, 0, 0)*3
regv=regCoef_n(pio,v_time, 0, 0)*3
printMinMax(regu,True)
copy_VarMeta(U(0,:,:),regu)
copy_VarMeta(V(0,:,:),regv)
;**********************MODEL************************

  yyyym  = datem/10000
  mmddm  = datem - yyyym*10000
  nSea  = num(mmddm.eq.mmddStrt)-1           ; number of seasons
  iSeam  = ind(mmddm.eq.mmddStrt)

  piom=new((/nDay*nSea/),"float")
  u_timem=new((/nDay*nSea,dimsizes(um&lat),dimsizes(um&lon)/),"float")
  v_timem=u_timem

  print(dimsizes(u_timem))
  print(dimsizes(v_timem))
    

  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSeam(ns)
    iLast  = iSeam(ns)+nDay-1
    ;print(iStrt)
    ;print(iLast)
    ;print(ns*180)
    ;print(180*(ns+1)-1)
    piom(ns*180:180*(ns+1)-1)=tofloat(PIOm(iStrt:iLast))
    u_timem(ns*180:180*(ns+1)-1,:,:)=Um(iStrt:iLast,:,:)
    v_timem(ns*180:180*(ns+1)-1,:,:)=Vm(iStrt:iLast,:,:)
  end do

piom=dim_standardize(piom,1)
regum=regCoef_n(piom,u_timem, 0, 0)*3
regvm=regCoef_n(piom,v_timem, 0, 0)*3
printMinMax(regum,True)
copy_VarMeta(Um(0,:,:),regum)
copy_VarMeta(Vm(0,:,:),regvm)
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
res@tmXBValues = fspan(40,180,7)
res@tmYLValues = fspan(-20,20,5)
res@tmYLLabels = (/"20~S~o~N~S","10~S~o~N~S","EQ","10~S~o~N~N","20~S~o~N~N"/)
res@tmXBLabels = (/"40~S~o~N","60~S~o~N~E","80~S~o~N~E","100~S~o~N~E","120~S~o~N~E","160~S~o~N~E","180~S~o~N"/)
res@tmYROn     = False
res@tmXTOn     = False
res@cnFillOn   = True
res@mpFillOn   = False
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
;res@cnLevels = (/-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7/)
;res@cnLevels = (/-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4/)
res@cnLevels = (/-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6/)
res@cnFillColors=(/15,30,45,60,0,170,190,210,230/)
res@lbLabelBarOn= False
res@pmLabelBarOrthogonalPosF= 0.15
res@lbLabelFontHeightF = 0.011
res@pmLabelBarWidthF = 0.5
res@pmLabelBarHeightF = 0.03

  res@tmXBLabelFontHeightF = 0.02                 ; resize tick labels
  res@tmYLLabelFontHeightF = 0.02
; change label spacing to avoid overlap
  res@lbLabelStride    = 0.1                        ; every other label
  res@tmXBTickSpacingF = 10.                      ; every 60 degree


res2=True
res2@gsnDraw = False
res2@gsnFrame = False
res2@gsnAddCyclic = False
res2@vcRefAnnoOrthogonalPosF = -0.4
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
res@cnLevels = (/-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6/)
res1@gsnContourNegLineDashPattern = 2
res1@cnLineThicknesses = True
res1@cnLineThicknessF = 2.
res1@cnLineLabelBackgroundColor=0
res1@gsnLeftString     = ""
res1@gsnRightString    =""


   ;a1=ndtooned(regum({-10:10},{80:150}))
   ;a2=ndtooned(regu({-10:10},{80:150}))
   ;cor=escorc(a1,a2)
   ;print(cor)
   ;rmse=dim_rmsd(a1,a2)
   ;y2=dim_stddev(a2)
   ;nrmse=rmse/y2
   ;print(nrmse)

plot = new(2,graphic)
plot(0) = gsn_csm_contour_map(wks,regu,res)
;plot0 = gsn_csm_contour_map(wks,regv,res1)
;overlay(plot(0),plot0)

;res@gsnLeftString           = "MODEL"
;res@gsnRightString    = "U850: PCC="+sprintf("%5.2f",cor)+" NRMSE="+sprintf("%5.2f",nrmse)
plot(1) = gsn_csm_contour_map(wks,regum,res)
;plot0 = gsn_csm_contour(wks,regvm,res1)
;overlay(plot(1),plot0)

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
  resP@txString            =  "Regressed U850 winds (m/s)"
  gsn_panel(wks,plot,(/2,1/),resP)
;---------------------------------

;------------------output------------------------


dirName = "./NC_OUT"
; Check if the directory exists and create it if it doesn't
dirCheck = systemfunc("test -d " + dirName + " ; echo $?")

if (dirCheck .eq. "0") then
    print("Directory " + dirName + " already exists.")
else
    system("mkdir -p " + dirName)
    print("Directory " + dirName + " created.")
end if

filename = ncDir + modout + ".mjo_U850_model." + pm + "mb.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->regum=regum

filename = ncDir + modout + ".mjo_U850_obs." + pm + "mb.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->regu=regu

filename = ncDir + modout + ".mjo_U850_model_all." + pm + "mb." + lonL_IO + "." + lonR_IO + ".nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->piom=piom
fout->u_timem=u_timem
fout->v_timem=v_timem

filename = ncDir + modout + ".mjo_U850_obs_all." + pm + "mb.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->pio=pio
fout->u_time=u_time
fout->v_time=v_time


end
;-----------------------------------------------
