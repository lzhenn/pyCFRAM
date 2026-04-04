load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;--------------------read co2----------------------------
  fname    = "/home/ys17-23/Extension2/zyh/AEROSOL/python_data/co2.txt"
  dataset  = asciiread(fname,-1,"string")
  url = "/home/ys17-23/Extension2/zyh/AEROSOL/python_data/"
  delim    = " "
  tfs      = str_fields_count(dataset, delim)

  year     = new((/dimsizes(tfs)/),"integer")
  co2      = new((/dimsizes(tfs)/),"double")
  year     = stringtointeger(str_get_field(dataset, 1, delim))
  co2      = stringtodouble(str_get_field(dataset, 2, delim))
  
  co2!0    = "year"
  co2&year = year

  print(year+" "+co2)
;-------------------data output-------------------------
  co2_b = avg(co2({year|1986:1995}))
  co2_w = avg(co2({year|2004:2013}))

;---Remove file just in case
  file1 = "co2_b.dat"
  system("rm -f " + file1)
;---Write a header to the file
;  fbindirwrite(file1,doubletofloat(co2_b))
 fbindirwrite(url+file1, co2_b)
 system("ls -l " + file1)

;---Remove file just in case
  file2 = "co2_w.dat"
  system("rm -f " + file2)
;---Write a header to the file
  fbindirwrite(url+file2,co2_w)
  system("ls -l " + file2)

; co2_b_read = fbindirread(file1,0,1,"double")
; print(co2_b_read)
;  print(co2_w)
; ;------------------figure output------------------------

;   wks      = gsn_open_wks("eps","CO2")

;   res                  = True                     ; plot mods desired
;   res@tiMainString     = "Global Annual Mean CO2 Concentration"          ; add title
;   res@tiYAxisString    = "ppm"
;   res@tiXAxisString    = "Year"

;   plot  = gsn_csm_xy (wks,year,co2,res) ; create plot

end


