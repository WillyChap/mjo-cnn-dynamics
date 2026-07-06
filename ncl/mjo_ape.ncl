;********************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;********************************************
  nameRegion = "EIO"        ; Indian Ocean base region
  latS_IO    = -10.
  latN_IO    =  10.
  lonL_IO    =  80.
  lonR_IO    = 100.

  latn       =  5.        ; lat band 
  lats       = -5.

  pltName    = "mjo_ape"       ; output plot name
  pltType    = "png"       ; x11, ps, eps, pdf, png
  pltDir     = getenv("FIG_DIR")        ; output plot directory
  ncDir      = getenv("NC_DIR")
  modout     = getenv("modout")

  dir   = "./obs/"
  uName  = "U"
  filu   = "ERA_u_1997-2013.nc"
  vName  = "V"
  filv   = "ERA_v_1997-2013.nc"
  wName  = "W"
  filw   = "ERA_w_1997-2013.nc"
  aName  = "T"
  fila   = "ERA_t_1997-2013.nc"
  pName  = "PR"
  filp   = "GPCP_1DD_v1.2_199610-201312.nc"

  twStrt = 19790101        ; Period for OBS
  twLast = 19951231        ; 

  cwStrtm = stringtoint(getenv("ymdStrt"))        ; Period for MODEL
  cwLastm = stringtoint(getenv("ymdLast"))

  dirm   = "../modeldata/NUIST/test/"
  film   = "day.nc"
  uNamem  = "U"
  aNamem  = "T"
  vNamem  = "V"
  wNamem  = "OMEGA"
  pNamem  = "PRECT"

  nameSeason = (/ "Winter"/)
  cp=1004.7
  R=287.
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
  f      = addfile(getenv("obs_data_p"), "r")
  date_p = cd_calendar(f->time, -2)   ; entire file
  iStrt  = ind(date_p.eq.twStrt)      ; desired dates
  iLast  = ind(date_p.eq.twLast)
  delete(date_p)
                                      ; P(time,lat,lon)
  P  = f->$pName$(iStrt:iLast,{-20:20},:)
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
  Pm  = f1->$pNamem$(iStrt1:iLast1,{-20:20},:)
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
;***********************U,V,W,T->ape*************************************
;OBS*******
  f      = addfile(getenv("obs_data_t"), "r")
  date_a = cd_calendar(f->time, -2)     ; entire file
  iStrt  = ind(date_a.eq.twStrt)        ; desired dates
  iLast  = ind(date_a.eq.twLast)
  delete(date_a)

;  lev=(/100,150,200,250,300,400,500,600,700,850,925,1000/)
 lev=(/1000,925,850,700,600,500,400,300,250,200,150,100/)
  AIR  = short2flt( f->$aName$(iStrt:iLast,{lev},{-10:10},:))
  lon =f->lon
  lat =f->lat({-10:10})
  time =AIR&time
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  nlev=dimsizes(lev)
; dTdt
   dTdt = center_finite_diff_n(AIR,1,False,0,0) ;units: K/day
   dTdt = dTdt /86400.                            ; units: K/s
   copy_VarMeta(AIR,dTdt)
; utx
  dtdx =  new(dimsizes(AIR),typeof(AIR),AIR@_FillValue)
  dlon = (lon(2)-lon(1))*0.0174533
  do nl=0,nlat-1                      ; loop over each latitude
      dX = 6378388.*cos(0.0174533*lat(nl))*dlon  ; constant at this latitude
      dtdx(:,:,nl:nl,:) = center_finite_diff_n (AIR(:,:,nl:nl,:), tofloat(dX) , False,0,3)
  end do
  f      = addfile(getenv("obs_data_u"), "r")
  U  = short2flt( f->$uName$(iStrt:iLast,{lev},{-10:10},:))
  utx = U * dtdx  
  copy_VarMeta(U,utx)
  delete([/U,dtdx,dX/]) 
; vty
  dtdy =  new(dimsizes(AIR),typeof(AIR),AIR@_FillValue)
  dlat = 0.0174533 * lat
  dY = 6378388.*dlat
  dtdy = center_finite_diff_n(AIR,tofloat(dY),False,0,2)
  f      = addfile(getenv("obs_data_v"), "r")
  V  = short2flt( f->$vName$(iStrt:iLast,{lev},{-10:10},:))
  vty = V * dtdy    ; units: K/s
  copy_VarMeta(V,vty)
  delete([/V,dtdy,dY/])
; wtp
  dtdp =  new(dimsizes(AIR),typeof(AIR),AIR@_FillValue)
  dtdp = center_finite_diff_n(AIR,tofloat(lev*100),False,0,1)
  f      = addfile(getenv("obs_data_w"), "r")
  W  = short2flt( f->$wName$(iStrt:iLast,{lev},{-10:10},:))
  wtp = W * dtdp 
  copy_VarMeta(W,wtp)
  delete([/dtdp/])
; wt
   wt  = W * AIR * R / cp 
   delete([/W/])
   do i = 0,nlev-1
     wt(:,i,:,:) = wt(:,i,:,:) / tofloat(lev(i))/100.0    ;K/s 
   end do  
   q1=dTdt+utx+vty+wtp-wt
    q1 = q1*AIR

   copy_VarMeta(dTdt,q1)
   q1@long_name = "diabatic heating rate"
   q1@units     = "K/s"
   q1&time = time
   printVarSummary(q1)
   q1=q1*86400 ;K/day
 
delete([/dTdt,utx,vty,wtp,wt/])
;MODEL******U,V,W,T->ape 
  f      = addfile(getenv("mod_data_daily"),"r")
  AIRm  = short2flt( f->$aNamem$(iStrt1:iLast1,{lev},{-10:10},:))
  lonm =f->lon
  latm =f->lat({-10:10})
  time =AIR&time
  nlatm=dimsizes(latm)
  nlonm=dimsizes(lonm)
  tim =AIRm&time
; dTdt
   dTdt = center_finite_diff_n(AIRm,1,False,0,0) ;units: K/day
   dTdt = dTdt /86400.                            ; units: K/s
   copy_VarMeta(AIRm,dTdt)
; utx
  dtdx =  new(dimsizes(AIRm),typeof(AIRm),AIRm@_FillValue)
  dlonm = (lonm(2)-lonm(1))*0.0174533
  do nl=0,nlatm-1                      ; loop over each latitude
      dX = 6378388.*cos(0.0174533*latm(nl))*dlonm  ; constant at this latitude
      dtdx(:,:,nl:nl,:) = center_finite_diff_n (AIRm(:,:,nl:nl,:), tofloat(dX) , False,0,3)
  end do
  U  = short2flt( f->$uNamem$(iStrt1:iLast1,{lev},{-10:10},:))
  utx = U * dtdx  
  copy_VarMeta(U,utx)
  delete([/U,dtdx,dX/]) 
; vty
  dtdy =  new(dimsizes(AIRm),typeof(AIRm),AIRm@_FillValue)
  dlatm = 0.0174533 * latm
  dY = 6378388.*dlatm
  dtdy = center_finite_diff_n(AIRm,tofloat(dY),False,0,2)
  V  = short2flt( f->$vNamem$(iStrt1:iLast1,{lev},{-10:10},:))
  vty = V * dtdy    ; units: K/s
  copy_VarMeta(V,vty)
  delete([/V,dtdy,dY/])
; wtp
  dtdp =  new(dimsizes(AIRm),typeof(AIRm),AIRm@_FillValue)
  dtdp = center_finite_diff_n(AIRm,tofloat(lev*100),False,0,1)
  W  = short2flt( f->$wNamem$(iStrt1:iLast1,{lev},{-10:10},:))
  wtp = W * dtdp 
  copy_VarMeta(W,wtp)
  delete([/dtdp/])
; wt
   wt  = W * AIRm * R / cp 
   delete([/W/])
   do i = 0,nlev-1
     wt(:,i,:,:) = wt(:,i,:,:) / tofloat(lev(i))/100.0    ;K/s 
   end do  
   q1m=dTdt+utx+vty+wtp-wt

   q1m =q1m *AIRm

   copy_VarMeta(dTdt,q1m)
   q1m@long_name = "diabatic heating rate"
   q1m@units     = "K/s"
   q1m&time = tim
   printVarSummary(q1m)
   q1m=q1m*86400 ;K/day

delete([/dTdt,utx,vty,wtp,wt/])
;************************************************
; Create LAT average of the global Precip and U series  (time,lon)
; Really, no need to area weight here
;************************************************
  the_timeLon   = dim_avg_Wrap(q1( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
;  the_timeLon   = dtrend_leftdim (the_timeLon, False)                ; rmv overall trend
  the_timeLon   = wgt_runave_leftdim( the_timeLon, wgt, 0 )          ; apply filter
  the_timeLon1   = dim_avg_Wrap(AIR( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
;  the_timeLon1   = dtrend_leftdim (the_timeLon1, False)                ; rmv overall trend
  the_timeLon1   = wgt_runave_leftdim( the_timeLon1, wgt, 0 )

  the_timeLonm   = dim_avg_Wrap(q1m( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
;  the_timeLonm   = dtrend_leftdim (the_timeLonm, False)                ; rmv overall trend
  the_timeLonm   = wgt_runave_leftdim(the_timeLonm, wgt, 0 )
  the_timeLonm1   = dim_avg_Wrap(AIRm( time|:,level|:,lon|:,{lat|lats:latn}) )  ; (time,lon)
;  the_timeLonm1   = dtrend_leftdim (the_timeLonm1, False)                ; rmv overall trend
  the_timeLonm1   = wgt_runave_leftdim(the_timeLonm1, wgt, 0 )
;**************winter**********************************
;OBS*************
  date=cd_calendar(time,-2)
  nDay = 180
  mmddStrt = 1101
  yyyy  = date/10000
  mmdd  = date - yyyy*10000
  nSea  = num(mmdd.eq.mmddStrt)-1           ; number of seasons
  iSea  = ind(mmdd.eq.mmddStrt)

  pio=new((/nDay*nSea/),"float")
  u_timelon=new((/nDay*nSea,nlev,dimsizes(the_timeLon&lon)/),"float")
  u_timelon1=new((/nDay*nSea,nlev,dimsizes(the_timeLon&lon)/),"float")
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSea(ns)
    iLast  = iSea(ns)+nDay-1
    pio(ns*180:180*(ns+1)-1)=tofloat(PIO(iStrt:iLast))
    u_timelon(ns*180:180*(ns+1)-1,:,:)=the_timeLon(iStrt:iLast,:,:)
    u_timelon1(ns*180:180*(ns+1)-1,:,:)=the_timeLon1(iStrt:iLast,:,:)
  end do
pio=dim_standardize(pio,1)
;u_timelon=dim_rmvmean_n(u_timelon,0)
;u_timelon1=dim_rmvmean_n(u_timelon1,0)
u_lon=regCoef_n(pio,u_timelon, 0, 0)*3
u_lon1=regCoef_n(pio,u_timelon1, 0, 0)*3
ape=u_lon/1000
copy_VarMeta(AIR(0,:,0,:),ape)
printMinMax(ape,True)
delete([/AIR,u_lon,u_lon1/])
;MODEL**********
  datem=cd_calendar(tim,-2)
  yyyym  = datem/10000
  mmddm  = datem - yyyym*10000
  nSea  = num(mmddm.eq.mmddStrt)-1           ; number of seasons
  iSeam  = ind(mmddm.eq.mmddStrt)

  piom=new((/nDay*nSea/),"float")
  u_timelonm=new((/nDay*nSea,nlev,dimsizes(the_timeLonm&lon)/),"float")
  u_timelonm1=new((/nDay*nSea,nlev,dimsizes(the_timeLonm&lon)/),"float")
  do ns=0,nSea-1                           ; seasonal averaging
    iStrt  = iSeam(ns)
    iLast  = iSeam(ns)+nDay-1
    piom(ns*180:180*(ns+1)-1)=tofloat(PIOm(iStrt:iLast))
    u_timelonm(ns*180:180*(ns+1)-1,:,:)=the_timeLonm(iStrt:iLast,:,:)
    u_timelonm1(ns*180:180*(ns+1)-1,:,:)=the_timeLonm1(iStrt:iLast,:,:)
  end do
piom=dim_standardize(piom,1)
;u_timelonm=dim_rmvmean_n(u_timelonm,0)
;u_timelonm1=dim_rmvmean_n(u_timelonm1,0)
u_lonm=regCoef_n(piom,u_timelonm, 0, 0)*3
u_lonm1=regCoef_n(piom,u_timelonm1, 0, 0)*3
apem=u_lonm;*u_lonm1
apem=apem/1000
copy_VarMeta(AIRm(0,:,0,:),apem)
printMinMax(apem,True)
delete([/AIRm,u_lonm,u_lonm1/])
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
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False                  ; turn on line labels
  res@cnFillOn             = True                  ; turn on color fill
  res@cnLevelSelectionMode ="ExplicitLevels"
  res@cnLevels = (/-0.02,0.02,0.3,0.5,0.7,0.9,1.1/)
  res@cnFillColors=(/80,0,165,180,195,210,225,240/)
  res@lbLabelStride = 1
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
  res3@cnLevels = (/-0.1,0.1/)
  res3@cnLineLabelsOn       = False
  res3@cnInfoLabelOn        = False
  res3@tiMainOn             = False
  res3@gsnLeftString        =""
  res3@gsnRightString    =""

   orp=ape({1000:100},{40:180})
   orp({1000:100},{85:95})=ape@_FillValue
   ;a1=ndtooned(apem({100000:10000},{40:180}))  
   ;a2=ndtooned(orp)
   ;cor=escorc(a1,a2)
   ;print(cor)
   ;rmse=dim_rmsd(a1,a2)
   ;y2=dim_stddev(a2)
   ;nrmse=rmse/y2
   ;print(nrmse)

  plot(0) = gsn_csm_pres_hgt(wks,ape({1000:100},{40:180}),res)
  plot3 = gsn_csm_contour(wks,ape,res3)
 ; overlay(plot(0),plot3)
  res@gsnLeftString    = "MODEL"
  ;res@gsnRightString    = "PCC="+sprintf("%5.2f",cor)+" NRMSE="+sprintf("%5.2f",nrmse)
  plot(1) = gsn_csm_pres_hgt(wks,apem({1000:100},{40:180}),res)

  resP=True
  resP@gsnPanelLabelBar     = False
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@txFont =25
  resP@txFontHeightF =0.018
  resP@txString            =  "Regressed APE generation rate (K~S~2~N~/day)"
  gsn_panel(wks,plot,(/2,1/),resP)
;-----------------------------

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

filename = "./NC_OUT/" + modout + ".mjo_ape.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->apem=apem


filename = "./NC_OUT/" + modout + ".mjo_ape_obs.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->ape=ape

filename = "./NC_OUT/" + modout + ".mjo_ape_model_all.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->q1m=u_timelonm
fout->piom=piom

filename = "./NC_OUT/" + modout + ".mjo_ape_obs_all.nc"
system("rm -f "+filename)
fout=addfile(filename,"c")
fout->pio=pio
fout->q1=u_timelon

;----------------------------------------------- 
