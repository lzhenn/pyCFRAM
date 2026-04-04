load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;************************************************
; read data
;************************************************
ntime = 1
nlev = 37
nlat = 1
nlon = 1

recl = ntime*(nlev+1)*nlat*nlon
dT_co2    = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_co2.dat",recl,"double")
dT_albedo = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_albedo.dat",recl,"double")
dT_cloud  = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_cloud.dat",recl,"double")
dT_o3     = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_o3.dat",recl,"double")
dT_q      = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_q.dat",recl,"double")
dT_solar  = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_solar.dat",recl,"double")
dT_warm   = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_warm.dat",recl,"double")
dT_ts     = fbinread("/home/ys17-23/Extension2/zyh/AEROSOL/python_data/data_output/dT_ts.dat",recl,"double")

dT_sum    = dT_co2 + dT_albedo + dT_cloud + dT_o3 + dT_q + dT_solar + dT_ts
dT_res    = dT_warm - dT_sum
dT_dyn    = -dT_warm

print(dT_warm)
print(dT_res)
print(dT_dyn)
print(dT_q)

end
