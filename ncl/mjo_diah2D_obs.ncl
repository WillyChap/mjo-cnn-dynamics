;********************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;********************************************
  nameRegion = "EIO"        ; Indian Ocean base region
  latS_IO    = -18
  latN_IO    =  18.
  lonL_IO    =  80.
  lonR_IO    =  100. 

  latn       =  5.        ; lat band 
  lats       = -5.

  pltName    = "mjo_diah2"       ; output plot name
  pltType    = "png"       ; x11, ps, eps, pdf, png
  pltDir     = getenv("FIG_DIR")        ; output plot directory
  ncDir      = getenv("NC_DIR")
  modout     = getenv("modout")

  dir   = "./OBS/"
  uName  = "U"
  vName  = "V"
  wName  = "W"
  aName  = "T"
  pName  = "PR"

  twStrt = 19790101        ; Period for OBS
  twLast = 19951231        ; 

  cwStrtm = stringtoint(getenv("ymdStrt"))        ; Period for MODEL
  cwLastm = stringtoint(getenv("ymdLast"))

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
  P  = f->$pName$(iStrt:iLast,{-28:28},{30:200})
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
  Pm  = f1->$pNamem$(iStrt1:iLast1,{-28:28},{30:200})
  Pm  =Pm*86400000
  printMinMax(Pm, True)

  wyIOm   = f1->lat({latS_IO:latN_IO})
  wyIOm   = cos(0.017459*wyIOm)
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
delete([/P,Pm/])
;***********************U,V,W,T->Q1*************************************
;OBS*******
  f      = addfile(getenv("obs_data_t"), "r")
  date_a = cd_calendar(f->time, -2)     ; entire file
  iStrt  = ind(date_a.eq.twStrt)        ; desired dates
  iLast  = ind(date_a.eq.twLast)
  delete(date_a)

;  lev=(/100,250,500,700,850,1000/)
 lev=(/1000,925,850,700,600,500,400,300,250,200,150,100/)

  AIR  = short2flt( f->$aName$(iStrt:iLast,{lev},{-28:28},{30:200}))
  printVarSummary(AIR)
  lon =f->lon
  lat =f->lat({-28:28})
  time =AIR&time
  nlat=dimsizes(lat)
  nlon=dimsizes(lon)
  nlev=dimsizes(lev)
  print("!!!!!!!!!")
  ;print(lat)
  print("!!!!!!!!!")
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
  U  = short2flt( f->$uName$(iStrt:iLast,{lev},{-28:28},{30:200}))
  utx = U * dtdx  
  copy_VarMeta(U,utx)
  delete([/U,dtdx,dX/]) 
; vty
  
  dtdy =  new(dimsizes(AIR),typeof(AIR),AIR@_FillValue)
  dlat = 0.0174533 * lat
  dY = 6378388.*dlat
  print(dimsizes(AIR))
  dtdy = center_finite_diff_n(AIR,tofloat(dY),False,0,2)
  f      = addfile(getenv("obs_data_v"), "r")
  V  = short2flt( f->$vName$(iStrt:iLast,{lev},{-28:28},{30:200}))
  vty = V * dtdy    ; units: K/s
  copy_VarMeta(V,vty)
  delete([/V,dtdy,dY/])
; wtp
  dtdp =  new(dimsizes(AIR),typeof(AIR),AIR@_FillValue)
  dtdp = center_finite_diff_n(AIR,tofloat(lev*100),False,0,1)
  f      = addfile(getenv("obs_data_w"), "r")
  W  = short2flt( f->$wName$(iStrt:iLast,{lev},{-28:28},{30:200}))
  wtp = W * dtdp 
  copy_VarMeta(W,wtp)
  delete([/dtdp/])
; wt
   wt  = W * AIR * R / cp 
   delete([/W,AIR/])
   do i = 0,nlev-1
     wt(:,i,:,:) = wt(:,i,:,:) / tofloat(lev(i))/100.0    ;K/s 
   end do  
   ;q1=dTdt+utx+vty+wtp-wt
   q1=dTdt+utx+vty+wtp-wt
   q1=smth9(q1, 0.50,  0.25, True)

   copy_VarMeta(dTdt,q1)
   copy_VarMeta(dTdt,wt)

   q1@long_name = "diabatic heating rate"
   q1@units     = "K/s"
   q1&time = time
   printVarSummary(q1)
   q1=q1 * 86400;K/day
   printMinMax(q1,True)

   filename = "./NC_OUT/" + modout + ".mjo_diah2d_ll_obs.nc"
   system("rm -f "+filename)
   fout=addfile(filename,"c")
   levout=(/925,700,300/)
   fout->q1_ll=q1(:,{levout},:,:)
   fout->wt_ll=wt(:,{levout},:,:)
delete([/q1,dTdt,utx,vty,wtp,wt/])

;MODEL******U,V,W,T->Q1 
  f      = addfile(getenv("mod_data_daily"),"r")
  AIRm  = short2flt( f->$aNamem$(iStrt1:iLast1,{lev},{-28:28},{30:200}))
  printVarSummary(AIRm)
  lonm =f->lon
  latm =f->lat({-28:28})
  print("!!!!!!!!!")
  ;print(latm)
  print("!!!!!!!!!")
  nlatm=dimsizes(latm)
  nlonm=dimsizes(lonm)
  nlevm=dimsizes(lev)
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
  U  = short2flt( f->$uNamem$(iStrt1:iLast1,{lev},{-28:28},{30:200}))
  utx = U * dtdx  
  copy_VarMeta(U,utx)
  delete([/U,dtdx,dX/]) 
; vty
  dtdy =  new(dimsizes(AIRm),typeof(AIRm),AIRm@_FillValue)
  dlatm = 0.0174533 * latm
  dY = 6378388.*dlatm
  print(dimsizes(AIRm))
  print("!!!!!!!!")
  print(dimsizes(tofloat(dY)))
  dtdy = center_finite_diff_n(AIRm,tofloat(dY),False,0,2)
  V  = short2flt( f->$vNamem$(iStrt1:iLast1,{lev},{-28:28},{30:200}))
  vty = V * dtdy    ; units: K/s
  copy_VarMeta(V,vty)
  delete([/V,dtdy,dY/])
; wtp
  dtdp =  new(dimsizes(AIRm),typeof(AIRm),AIRm@_FillValue)
  dtdp = center_finite_diff_n(AIRm,tofloat(lev*100),False,0,1)
  W  = short2flt( f->$wNamem$(iStrt1:iLast1,{lev},{-28:28},{30:200}))
  wtp = W * dtdp 
  copy_VarMeta(W,wtp)
  delete([/dtdp/])
; wt
   wt  = W * AIRm * R / cp 
   delete([/W,AIRm/])
   do i = 0,nlev-1
      wt(:,i,:,:) = wt(:,i,:,:) / tofloat(lev(i))/100.0    ;K/s 
      ;wt(:,i,:,:) = wt(:,i,:,:) / tofloat(lev(i))    ;K/s 
   end do  
   
  ;q1m=dTdt+utx+vty+wtp-wt   
   q1m=dTdt+utx+vty+wtp-wt
   copy_VarMeta(dTdt,q1m)
   q1m@long_name = "diabatic heating rate"
   q1m@units     = "K/s"
   q1m&time = tim
   printVarSummary(q1m)
   q1m=q1m * 86400 ;K/day
   printMinMax(q1m,True)
   filename = "./NC_OUT/" + modout + ".mjo_diah2d_ll_model.nc"
   system("rm -f "+filename)
   fout=addfile(filename,"c")

   ; Specifying NetCDF4 format
   setfileoption(fout, "Format", "NetCDF4")


   levout=(/1000,700,300/)
   fout->q1m_ll=q1m(:,{levout},:,:)

   copy_VarMeta(dTdt,wt)
   fout->wtm_ll=wt(:,{levout},:,:)

   ;copy_VarMeta(dTdt,wtp)
   ;fout->wtpm_ll=wtp(:,{levout},:,:)

   ;copy_VarMeta(dTdt,vty)
   ;fout->vtym_ll=vty(:,{levout},:,:)

   ;copy_VarMeta(dTdt,utx)
   ;fout->utxm_ll=utx(:,{levout},:,:)

   ;fout->dTdtm_ll=dTdt(:,{levout},:,:)


   delete([/q1m, dTdt,utx,vty,wtp,wt/])

;----------------------------------------------- 
