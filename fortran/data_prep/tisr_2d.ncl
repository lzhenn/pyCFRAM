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
; EP & CP .dat
;************************************************
  f1     = addfile("/home/ys17-23/Extension/ERA5/ERA5-monthly/surface/toa_solarin_mon_1x1_1979_2019.nc","r")
  url = "/home/ys17-23/Extension2/zyh/AEROSOL/python_data/"
  TIME   = f1->time
  YYYY   = cd_calendar(TIME,-1)/100                 ; entire file
  iYYYY1  = ind(YYYY.ge.yrStrt1 .and. YYYY.le.yrLast1)
  iYYYY2  = ind(YYYY.ge.yrStrt2 .and. YYYY.le.yrLast2)

  x1_raw  = short2flt(f1->tisr(iYYYY1,:,:))  ;Pa
  x1      = x1_raw/86400.
  x2_raw  = short2flt(f1->tisr(iYYYY2,:,:))  ;Pa
  x2      = x2_raw/86400.
 
  lat = f1->latitude
  lon = f1->longitude
  Time1= f1->time(iYYYY1)
  time1= ut_calendar(Time1,-1) 
  Time2= f1->time(iYYYY2)
  time2= ut_calendar(Time2,-1) 

  x1!0 = "time"
  x1!1 = "lat"
  x1!2 = "lon"
  x1&time = time1
  x1&lat  = lat
  x1&lon  = lon
  x1&lat@units = "degrees_north"
  x1&lon@units = "degrees_east"
  copy_VarMeta(x1,x2)
  x2&time = time2

  X1      = month_to_annual (x1, 1)     ;for ssrd
  X2      = month_to_annual (x2, 1)     ;for ssru
  print(X1(:,0,0))

  year1    = ispan(yrStrt1,yrLast1,1)
  year2    = ispan(yrStrt2,yrLast2,1)
  X1&year = year1
  X2&year = year2

  n      = 2
  nx     = dimsizes(lat)
  ny     = dimsizes(lon)

  X_b   = new((/n,10,nx,ny/),typeof(X1))

  X_b(0,:,:,:) = X1
  X_b(1,:,:,:) = X2

  X_avg   = dim_avg_n_Wrap(X_b,1) ; [case,lat,lon]

  printVarSummary(X_avg)
  print(X_avg(0,:,11))

;---Remove file just in case
  file1 = "solarin_base.dat"
  system("rm -f " + file1)
;---Write first variable to file
  fbindirwrite(url+file1,flt2dble(X_avg(0,:,:)))
  system("ls -l " + file1)

;---Remove file just in case
  file2 = "solarin_warm.dat"
  system("rm -f " + file2)
;---Write first variable to file
  fbindirwrite(url+file2,flt2dble(X_avg(1,:,:)))
  system("ls -l " + file2)

;************************************************
; create plots
;************************************************
;   wks = gsn_open_wks("pdf","solin")               ; open a ps file
;   gsn_define_colormap(wks,"MPL_YlGnBu")
;   plot = new(n,graphic)                          ; create a plot array

;   res                  = True                      ; plot mods desired
;   res@gsnFrame         = False
;   res@gsnDraw          = False
;   res@gsnSpreadColors  = True                      ; use full colormap

;   res@mpOutlineOn      = True
;   res@mpFillOn         = True
;   res@mpMinLatF        = -90.	; zoom in on a subregion
;   res@mpMaxLatF        = 90.
;   res@mpMinLonF        = 0.
;   res@mpCenterLonF     = 180.
;   res@mpMaxLonF        = 360.
;   res@mpGeophysicalLineThicknessF=1.3

;   res@cnFillOn         = True                      ; turn on color
;   res@cnLinesOn        = False
;   res@cnLineLabelsOn   = False                     ; no contour line labels
;   res@cnInfoLabelOn    = False                     ; no info label
;   res@lbLabelBarOn     = False

; ;  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; ;  res@cnMinLevelValF       = 0.
; ;  res@cnMaxLevelValF       = 1.
; ;  res@cnLevelSpacingF      = .1

; ; now change the size of the tickmark labels
; ;  res@tmXBLabelFontHeightF = 0.025                 ; resize tick labels
; ;  res@tmYLLabelFontHeightF = 0.025
;   res@gsnRightString       =  "1" ; "W/m~S~2~N~"

;   do i = 0,n-1
;      res@gsnLeftString     = "Solin(case "+tostring(i)+")"
;      plot(i) = gsn_csm_contour_map_ce(wks,X_avg(i,:,:),res)
;   end do

; ;************************************************
; ; create panel
; ;************************************************
;   resP                     = True
;   resP@gsnFrame            = False
;   resP@gsnMaximize         = True
;   resP@gsnPanelLabelBar    = True                ; add common colorbar
;   gsn_panel(wks,plot,(/3,3/),resP)
;   frame(wks)

; ;************************************************
; ; second page
; ;************************************************
;   do i = 0,n-1
; ;     res@gsnLeftString     = "SSRD("+warm_s(i)+"-"+warm_e(i)+")"
; ;     plot(i) = gsn_csm_contour_map_ce(wks,x1_w(i,:,:),res)
;   end do
;   gsn_panel(wks,plot,(/3,3/),resP)
;   frame(wks)

; ;************************************************
; ; third page
; ;************************************************
;    do i = 0,n-1
; ;      res@gsnLeftString     = "SSRU("+base_s(i)+"-"+base_e(i)+")"
; ;      plot(i) = gsn_csm_contour_map_ce(wks,y1_b(i,:,:),res)
;    end do
;    gsn_panel(wks,plot,(/3,3/),resP)
;    frame(wks)

; ;************************************************
; ; forth page
; ;************************************************
;   do i = 0,n-1
; ;     res@gsnLeftString     = "SSRU("+warm_s(i)+"-"+warm_e(i)+")"
; ;     plot(i) = gsn_csm_contour_map_ce(wks,y1_w(i,:,:),res)
;   end do
; ;  gsn_panel(wks,plot,(/3,3/),resP)
;   frame(wks)

; ;************************************************
; ; third page
; ;************************************************
;   gsn_define_colormap(wks,"BlueWhiteOrangeRed")

; ;  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; ;  res@cnMinLevelValF       = -.1
; ;  res@cnMaxLevelValF       = .1
; ;  res@cnLevelSpacingF      =   0.01

;   do i = 0,n-1
; ;     res@gsnLeftString     = "Delta SSRD"+warm_s(i)+"/"+base_s(i)+")"
; ;     plot(i) = gsn_csm_contour_map_ce(wks,d_x1(i,:,:),res)
;   end do
; ;  gsn_panel(wks,plot,(/3,3/),resP)
;   frame(wks)

;   do i = 0,n-1
; ;     res@gsnLeftString     = "Delta SSRU"+warm_s(i)+"/"+base_s(i)+")"
; ;     plot(i) = gsn_csm_contour_map_ce(wks,d_y1(i,:,:),res)
;   end do
; ;  gsn_panel(wks,plot,(/3,3/),resP)
;   frame(wks)


; delete(wks)
end
