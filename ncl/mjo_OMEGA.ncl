;********************************************
; These files are loaded by default in NCL V6.2.0 and newer
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./vertical_thete.ncl"
;********************************************
  nameRegion = "EIO"        ; Indian Ocean base region
  latS_IO    = -10.
  latN_IO    =  10.
  lonL_IO    =  80.
  lonR_IO    = 100.

  latn       =  5.        ; lat band for (lag,lon) Fig 5
  lats       = -5.

  pltName    = "mjo_Ws"       ; output plot name
  pltType    = "png"       ; x11, ps, eps, pdf, png
  pltDir     = getenv("FIG_DIR")        ; output plot directory
  ncDir      = getenv("NC_DIR")
  modout     = getenv("modout")

  dir   = "./OBS/"
  pName  = "PR"
  filp   = "GPCP_1DD_v1.2_199610-201312.nc"
  rName  = "W"
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
  rNamem  = "OMEGA"

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
  P  = f->$pName$(iStrt:iLast,{-10:10},:)
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
   p=(/1000,925,850,700,600,500,400,300,250,200,150,100/)
  f      = addfile(getenv("obs_data_w"), "r")
  print(f)
  SH  = short2flt( f->$rName$(iStrt:iLast,{p},{-10:10},:))
  f      = addfile(getenv("obs_data_t"), "r")
  AIR  = short2flt( f->$aName$(iStrt:iLast,{p},{-10:10},:))
 
  nlev=dimsizes(p)
 ; X = mixhum_ptrh (conform(AIR,p,1), AIR, RH, 1); convert rel humidity to mixing ratio
  p0 = 1000
  p0_at_units = "hPa" ; must match units of p
;MODEL******* 
 AIRm  = f1->$aNamem$(iStrt1:iLast1,{p},{-10:10},:)
 printMinMax(AIRm,True)
 SHm  = f1->$rNamem$(iStrt1:iLast1,{p},{-10:10},:)
 printMinMax(SHm,True)
 p0 = 1000
 p0_at_units = "Pa" ; must match units of p
;************************************************
; Create LAT average of the global Precip and U series  (time,lon)
; Really, no need to area weight here
;************************************************
  the_timeLon   = dim_avg_Wrap(SH( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
 ; the_timeLon   = dtrend_leftdim (the_timeLon, False)                ; rmv overall trend
  the_timeLon   = wgt_runave_leftdim( the_timeLon, wgt, 0 )          ; apply filter

  the_timeLonm   = dim_avg_Wrap(SHm( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
 ; the_timeLonm   = dtrend_leftdim (the_timeLonm, False)                ; rmv overall trend
  the_timeLonm   = wgt_runave_leftdim(the_timeLonm, wgt, 0 )
;**************winter**********************************
;OBS*************
  date=cd_calendar(SH&time,-2)
  nDay = 180
  mmddStrt = 1101
  yyyy  = date/10000
  mmdd  = date - yyyy*10000
  nSea  = num(mmdd.eq.mmddStrt)-1           ; number of seasons
  iSea  = ind(mmdd.eq.mmddStrt)

  pio=new((/nDay*nSea/),"float")
  u_timelon=new((/nDay*nSea,nlev,dimsizes(the_timeLon&lon)/),"float")
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSea(ns)
    iLast  = iSea(ns)+nDay-1
    pio(ns*180:180*(ns+1)-1)=tofloat(PIO(iStrt:iLast))
    u_timelon(ns*180:180*(ns+1)-1,:,:)=the_timeLon(iStrt:iLast,:,:)
  end do
pio=dim_standardize(pio,1)
u_lon=regCoef_n(pio,u_timelon, 0, 0)*100
copy_VarMeta(SH(0,:,0,:),u_lon)
printMinMax(u_lon,True)
;MODEL**********
  datem=cd_calendar(SHm&time,-2)
  yyyym  = datem/10000
  mmddm  = datem - yyyym*10000
  nSea  = num(mmddm.eq.mmddStrt)-1           ; number of seasons
  iSeam  = ind(mmddm.eq.mmddStrt)

  piom=new((/nDay*nSea/),"float")
  u_timelonm=new((/nDay*nSea,nlev,dimsizes(the_timeLonm&lon)/),"float")
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSeam(ns)
    iLast  = iSeam(ns)+nDay-1
    piom(ns*180:180*(ns+1)-1)=tofloat(PIOm(iStrt:iLast))
    u_timelonm(ns*180:180*(ns+1)-1,:,:)=the_timeLonm(iStrt:iLast,:,:)
  end do
piom=dim_standardize(piom,1)
u_lonm=regCoef_n(piom,u_timelonm, 0, 0)*100
printMinMax(u_lonm,True)
copy_VarMeta(SHm(0,:,0,:),u_lonm)
; ===========================
; color plot
; ===========================
  plot = new( 2, "graphic")
  pltNameSeason = modout+"."+pltName+"."+nameSeason(0)
  wks = gsn_open_wks(pltType,pltDir+"/"+pltNameSeason) 
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")

  res                      = True                  ; plot mods desired
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
  res@cnLevelSpacingF      = 0.1                   ; contour interval
  res@cnMinLevelValF       = -1.5                ; min level
  res@cnMaxLevelValF       =  1.5                  ; max level
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False                  ; turn on line labels
  res@cnFillOn             = True                  ; turn on color fill
  res@lbLabelStride = 2
  res@pmLabelBarWidthF = 0.6
  res@pmLabelBarHeightF = 0.06
  res@lbLabelFontHeightF = 0.02
  res@lbLabelFont = 25
  res@pmLabelBarOrthogonalPosF= 0.02   
  res@gsnLeftString        ="OBS"
  res@gsnRightString=""

  res3                      = True
  res3@gsnDraw              = False              ; don't draw
  res3@gsnFrame             = False              ; don't advance frame
  res3@gsnAddCyclic         = True
  res3@cnFillOn             = False
  res3@cnLineColor          ="blue"
  res3@cnLineThicknessF     = 1.5
  res3@gsnContourNegLineDashPattern = 1
  res3@cnLevelSelectionMode = "ExplicitLevels"
  res3@cnLevels = (/-0.3,0.3/)
  res3@cnLineLabelsOn       = False
  res3@cnInfoLabelOn        = False
  res3@tiMainOn             = False
  res3@gsnLeftString        =""

;   orp=where(abs(u_lon({1000:250},{40:180})).ge.0.3,u_lon({1000:250},{40:180}),u_lon@_FillValue)
   orp=u_lon({1000:100},{40:180})
   orp({1000:100},{85:95})=u_lon@_FillValue
   a1=ndtooned(u_lonm({1000:100},{40:180}))  
   a2=ndtooned(orp)
   ;cor=escorc(a1,a2)
   ;print(cor)
   ;rmse=dim_rmsd(a1,a2)
   ;y2=dim_stddev(a2)
   ;nrmse=rmse/y2
   ;print(nrmse)

  plot(0) = gsn_csm_pres_hgt(wks,u_lon({1000:100},{40:180}),res)
  plot3 = gsn_csm_contour(wks,u_lon,res3)
  overlay(plot(0),plot3)
  res@gsnLeftString    = "MODEL"
  ;res@gsnRightString    = "PCC="+sprintf("%5.2f",cor)+" NRMSE="+sprintf("%5.2f",nrmse)
  plot(1) = gsn_csm_pres_hgt(wks,u_lonm({1000:100},{40:180}),res)

  resP=True
  resP@gsnPanelLabelBar     = False
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@txFont =25
  resP@txFontHeightF =0.018
  resP@txString            =  "Regressed OMEGA"
  gsn_panel(wks,plot,(/2,1/),resP)
;-------------------------------

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

filename = "./NC_OUT/" + modout + ".mjo_W_model.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->THETAEm=u_timelonm
fout->piom=piom

filename = "./NC_OUT/" + modout + ".mjo_W_obs.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->THETAE=u_timelon
fout->pio=pio

filename = "./NC_OUT/" + modout + ".mjo_reg_W.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->u_lonm=u_lonm

filename = "./NC_OUT/" + modout + ".mjo_reg_W_obs.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->u_lon=u_lon
;-----------------------------------------------

 
