load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;************************************************
; read data
;************************************************
  fin_lw_base = addfile("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/IN_AER_RRTM_LW.nc","r")
  fin_sw_base = addfile("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/IN_AER_RRTM_SW.nc","r")
  aod_lw_base = fin_lw_base->aod
  aod_sw_base = fin_sw_base->aod
  ssa_sw_base = fin_sw_base->ssa
  g_sw_base   = fin_sw_base->g

  print(aod_sw_base(:,0))

  fin_lw_warm = addfile("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/IN_AER_RRTM_LW.ncc","r")
  fin_sw_warm = addfile("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/IN_AER_RRTM_SW.nc","r")
  aod_lw_warm = fin_lw_warm->aod
  aod_sw_warm = fin_sw_warm->aod
  ssa_sw_warm = fin_sw_warm->ssa
  g_sw_warm   = fin_sw_warm->g
  
  url = "/home/ys17-23/Extension2/zyh/AEROSOL/python_data/"
;---Remove file just in case
  file1 = "aerosol_aod_lw_base.dat"
  file2 = "aerosol_aod_sw_base.dat"
  file3 = "aerosol_ssa_sw_base.dat"
  file4 = "aerosol_g_sw_base.dat"

  file5 = "aerosol_aod_lw_warm.dat"
  file6 = "aerosol_aod_sw_warm.dat"
  file7 = "aerosol_ssa_sw_warm.dat"
  file8 = "aerosol_g_sw_warm.dat"

  system("rm -f " + file1)
  system("rm -f " + file2)
  system("rm -f " + file3)
  system("rm -f " + file4)
  system("rm -f " + file5)
  system("rm -f " + file6)
  system("rm -f " + file7)
  system("rm -f " + file8)

;---Write first variable to file
  fbindirwrite(url+file1,flt2dble(aod_lw_base))
  fbindirwrite(url+file2,flt2dble(aod_sw_base))
  fbindirwrite(url+file3,flt2dble(ssa_sw_base))
  fbindirwrite(url+file4,flt2dble(g_sw_base))
  fbindirwrite(url+file5,flt2dble(aod_lw_warm))
  fbindirwrite(url+file6,flt2dble(aod_sw_warm))
  fbindirwrite(url+file7,flt2dble(ssa_sw_warm))
  fbindirwrite(url+file8,flt2dble(g_sw_warm))

  system("ls -l " + file1)
  system("ls -l " + file2)
  system("ls -l " + file3)
  system("ls -l " + file4)
  system("ls -l " + file5)
  system("ls -l " + file6)
  system("ls -l " + file7)
  system("ls -l " + file8)

end
