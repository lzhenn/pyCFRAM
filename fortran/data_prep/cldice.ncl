load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  yrStrt1 = 1986
  yrLast1 = 1995
  yrStrt2 = 2004
  yrLast2 = 2013

;************************************************
; read data
;************************************************
  f1 = addfile("/home/ys17-23/Extension/ERA5/ERA5-monthly/pressure/cloud_ice_water_content_mon_1x1_1979_2019.nc","r")
  url = "/home/ys17-23/Extension2/zyh/AEROSOL/python_data/"
  TIME   = f1->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY1  = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
  iYYYY2  = ind(YYYY.ge.yrStrt2 .and. YYYY.le.yrLast2)

  var1  = short2flt(f1->ciwc(iYYYY1,:,:,:))  ;Pa
  var2  = short2flt(f1->ciwc(iYYYY2,:,:,:))  ;Pa
  
  lat  = f1->latitude
  lon  = f1->longitude
  lev  = f1->level
  Time = f1->time(iYYYY1)
  printVarSummary(var1)
  delete(f1)
  
  var1_annual = month_to_annual(var1,1)
  var2_annual = month_to_annual(var2,1)
  year1    = ispan(yrStrt1,yrLast1,1)
  year2    = ispan(yrStrt2,yrLast2,1)
  var1_annual&year = year1
  var2_annual&year = year2

  printVarSummary(var1_annual)
  
  n  = 2                        ;Returns the dimension sizes of the input variable
  nx = dimsizes(lon)
  ny = dimsizes(lat)
  nz = dimsizes(lev)

  var_b = new((/2,10,nz,ny,nx/),typeof(var1_annual))
  
  var_b(0,:,:,:,:) = var1_annual
  var_b(1,:,:,:,:) = var2_annual
  
  var_avg  = dim_avg_n_Wrap(var_b,1) ;[case,level,lat,lon]

  print(var_avg(0,:,11,11))

;---Remove file just in case
  file1 = "ciwc_base.dat"
  system("rm -f " + file1)
;---Write first variable to file
  fbindirwrite(url+file1,flt2dble(var_avg(0,:,0,0)))
  system("ls -l " + file1)

;---Remove file just in case
  file2 = "ciwc_warm.dat"
  system("rm -f " + file2)
;---Write first variable to file
  fbindirwrite(url+file2,flt2dble(var_avg(1,:,0,0)))
  system("ls -l " + file2)


;   var_zonal = dim_avg_Wrap(var_avg)  ; [case,level,lat]
; ;********************************************************
; ; Create  Plot page1
; ;********************************************************
;   wks   = gsn_open_wks ("eps","cldice")             ; open workstation
;   gsn_define_colormap(wks,"BlueWhiteOrangeRed")             ; choose colormap
;   plot= new(n,graphic)

;   res                      = True
;   res@gsnFrame             = False
;   res@gsnDraw              = False                ; plot mods desired
;   res@cnLineLabelsOn       = True                 ; turn on line labels
;   res@cnFillOn             = True                 ; turn on color fill
;   res@cnLinesOn        = False
;   res@cnLineLabelsOn   = False                    ; no contour line labels
;   res@cnInfoLabelOn    = False
; ;---This resource not needed in V6.1.0
;   res@gsnSpreadColors      = True                 ; use full range of map

; ;---This resource defaults to True in NCL V6.1.0
;   res@lbLabelAutoStride    = True                 ; optimal labels
;   res@tmXBLabelFontHeightF = 0.027                ; resize tick labels
;   res@tmYLLabelFontHeightF = 0.027
;   res@tmYRMode             = "Automatic"
;   res@tmYLMode   = "Explicit"  		          ; explicitly set Y-axis labels
;   res@tmYLValues = (/"1000","850","750","650","500","300","100","70","50"/)
;   res@tmYLLabels = (/"1000","850","750","650","500","300","100","70","50"/)
;   res@tmYUseLeft = False		; use different resources for the right Y-axis and the left Y-axis.
;   res@tmYRLabelFontColor = 18

; ;  res@tmXBTickSpacingF   = 60.
;   res@lbLabelBarOn        = False            ; turn off individual cb's
; ;  res@lbLabelStride       = 2                        ; every other label
; ;  res@lbLabelFontHeightF  = 0.027
; ;  res@pmLabelBarHeightF   = 0.1
; ;  res@pmLabelBarWidthF    = 0.6
; ;  res@lbLabelStrings      =(/"-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7"/)

; ;  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; ;  res@cnMinLevelValF       = -7e-4                ; set min contour level
; ;  res@cnMaxLevelValF       = 7e-4               ; set max contour level
; ;  res@cnLevelSpacingF      =  1e-4                ; set contour spacing

; ;  res@txFontHeightF   = 0.025	; set font height of gsn*String
;   res@gsnLeftStringFontHeightF    = 0.03
;   res@gsnRightStringFontHeightF   = 0.025
;   res@gsnRightString              = "1"   ;"*20~S~-4~N~ kgkg~S~**-1~N~"

;   do i = 0,n-1
;      res@gsnLeftString  = "t case("+tostring(i)+")"
;      plot(i)  = gsn_csm_pres_hgt(wks, var_zonal(i,{1000:50},:), res )   ; plaace holder
;   end do

; ;************************************************
; ; create panel
; ;************************************************
;   resP                            = True
;   resP@gsnFrame                   = True
;   resP@gsnMaximize                = True
;   resP@gsnPanelLabelBar           = True                ; add common colorbar
; ;  resP@gsnPanelYWhiteSpacePercent = 4
; ;  resP@gsnPanelXWhiteSpacePercent = 2
; ;  resP@gsnPanelScalePlotIndex     = 5
;   gsn_panel(wks,plot,(/3,3/),resP)

; ;************************************************
; ; second page
; ;************************************************
;   do i = 0,n-1
; ;     res@gsnLeftString     = "o3 ("+warm_s(i)+"-"+warm_e(i)+")"
; ;     plot(i)  = gsn_csm_pres_hgt(wks, cc_w_zonal(i,{1000:50},:), res )
;   end do
; ;  gsn_panel(wks,plot,(/3,3/),resP)

; ;************************************************
; ; third page
; ;************************************************
; ;  gsn_define_colormap(wks,"MPL_BrBG")
;   do i = 0,n-1
; ;     res@gsnLeftString     = "Delta o3 ("+warm_s(i)+"/"+base_s(i)+")"
; ;     plot(i)  = gsn_csm_pres_hgt(wks, dc_zonal(i,{1000:50},:), res )
;   end do
; ;  gsn_panel(wks,plot,(/3,3/),resP)

; ;************************************************
; ; forth page(x-y plot)
; ;************************************************
;   lev_x = (/1000,925,850,800,750,500,300,200,100/)
;   n = dimsizes(lev_x)
;   res@mpCenterLonF      = 180
;   do i = 0,n-1
;      res@gsnLeftString     = lev_x(i)+ " hPa 1986-1981" ; Specific Humidity("+warm_s(i)+"-"+warm_e(i)+")"
; ;     plot(i) = gsn_csm_contour_map_ce(wks,d_cc(0,{lev_x(i)},:,:),res)
;   end do
;   gsn_panel(wks,plot,(/3,3/),resP)

end
