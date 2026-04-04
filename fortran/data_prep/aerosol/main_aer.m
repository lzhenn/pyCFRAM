clear;clc;
%This is the main program for running aerosol calculations. 
%It generates aerosol files by calling the longwave and shortwave aer_writer modules.
%%
direc = '/home/ys17-23/Extension2/zyh/AEROSOL/';%The path for storing the optbands files and output aerosol files
addpath('/home/yangsong3/yh/4_model/CFRAM_A/matlab/');%path of the functions

datainfo = ncinfo('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4');%path of aerosol data
lon = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','lon');
lat = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','lat');
tie = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','time');

%---------get the relevant variable fields
%type of aerosols: (units:kg/kg)
%BC:blackcarbon         DU:dust
%OC: organic carbon     SS:sea salt
%SU: sulfate
%other variables:
%dens:Air Density(kg/m3)
%rh:Relative Humidity
%delp:Pressure Difference(Pa)
dens = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','AIRDENS', ...
       [289 181 1 3],[1 1 Inf 1]);
rh   = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','RH', ...
       [289 181 1 3],[1 1 Inf 1]);
delp = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DELP', ...
       [289 181 1 3],[1 1 Inf 1]);

DU001 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DU001', ...
       [289 181 1 3],[1 1 Inf 1]);
DU002 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DU002', ...
       [289 181 1 3],[1 1 Inf 1]);
DU003 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DU003', ...
       [289 181 1 3],[1 1 Inf 1]);
DU004 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DU004', ...
       [289 181 1 3],[1 1 Inf 1]);
DU005 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','DU005', ...
       [289 181 1 3],[1 1 Inf 1]);
SS001 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SS001', ...
       [289 181 1 3],[1 1 Inf 1]);
SS002 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SS002', ...
       [289 181 1 3],[1 1 Inf 1]);
SS003 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SS003', ...
       [289 181 1 3],[1 1 Inf 1]);
SS004 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SS004', ...
       [289 181 1 3],[1 1 Inf 1]);
SS005 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SS005', ...
       [289 181 1 3],[1 1 Inf 1]);
SO4 = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','SO4', ...
       [289 181 1 3],[1 1 Inf 1]);
BCPHOBIC = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','BCPHOBIC', ...
       [289 181 1 3],[1 1 Inf 1]);
BCPHILIC = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','BCPHILIC', ...
       [289 181 1 3],[1 1 Inf 1]);
OCPHOBIC = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','OCPHOBIC', ...
       [289 181 1 3],[1 1 Inf 1]);
OCPHILIC = ncread('/home/ys17-23/Extension2/zyh/AEROSOL/MERRA2_400.inst3_3d_aer_Nv.20211201.nc4','OCPHILIC', ...
       [289 181 1 3],[1 1 Inf 1]);

faer.DU001 = squeeze(DU001);
faer.DU002 = squeeze(DU002);
faer.DU003 = squeeze(DU003);
faer.DU004 = squeeze(DU004);
faer.DU005 = squeeze(DU005);
faer.SS001 = squeeze(SS001);
faer.SS002 = squeeze(SS002);
faer.SS003 = squeeze(SS003);
faer.SS004 = squeeze(SS004);
faer.SS005 = squeeze(SS005);
faer.SO4 = squeeze(SO4);
faer.BCPHOBIC = squeeze(BCPHOBIC);
faer.BCPHILIC = squeeze(BCPHILIC);
faer.OCPHILIC = squeeze(OCPHILIC);
faer.OCPHOBIC = squeeze(OCPHOBIC);

%---------Call the function to generate the file.
aer_writer_lw(direc,faer,squeeze(rh),squeeze(delp),squeeze(dens));
aer_writer_sw(direc,faer,squeeze(rh),squeeze(delp),squeeze(dens));

%%