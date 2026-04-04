%%%%%%%Call the aer_opt file to compute aerosol-related properties,
%%%%%%%and store the computation results in a netCDF file.
function aer_writer_lw(direc,faer,rh,delp,dens)

%mdpn = [180,425,565,665,760,900,1030,1130,1285,1435,1640,1940,2165,2315,2490,2925];
%mdpl = 1./mdpn.*1e-2;

thick = delp./(dens.*9.8);

[aod_du1] = aer_opt_lw(direc,faer.DU001,'opticsBands_DU.v15_3.RRTMG.nc',1,rh,thick,dens);
[aod_du2] = aer_opt_lw(direc,faer.DU002,'opticsBands_DU.v15_3.RRTMG.nc',2,rh,thick,dens);
[aod_du3] = aer_opt_lw(direc,faer.DU003,'opticsBands_DU.v15_3.RRTMG.nc',3,rh,thick,dens);
[aod_du4] = aer_opt_lw(direc,faer.DU004,'opticsBands_DU.v15_3.RRTMG.nc',4,rh,thick,dens);
[aod_du5] = aer_opt_lw(direc,faer.DU005,'opticsBands_DU.v15_3.RRTMG.nc',5,rh,thick,dens);
[aod_sa]    = aer_opt_lw(direc,faer.SO4, 'opticsBands_SU.v1_3.RRTMG.nc', 0.16*1e-6, rh,thick,dens);
[aod_bco] = aer_opt_lw(direc,faer.BCPHOBIC, 'opticsBands_BC.v1_3.RRTMG.nc', 1, rh,thick,dens);
[aod_bci] = aer_opt_lw(direc,faer.BCPHILIC, 'opticsBands_BC.v1_3.RRTMG.nc', 2, rh,thick,dens);
[aod_oco] = aer_opt_lw(direc,faer.OCPHOBIC, 'opticsBands_OC.v1_3.RRTMG.nc', 1, rh,thick,dens);
[aod_oci] = aer_opt_lw(direc,faer.OCPHILIC, 'opticsBands_OC.v1_3.RRTMG.nc', 2, rh,thick,dens);
[aod_ss1] = aer_opt_lw(direc,faer.SS001, 'opticsBands_SS.v3_5.RRTMG.nc', 1, rh,thick,dens); 
[aod_ss2] = aer_opt_lw(direc,faer.SS002, 'opticsBands_SS.v3_5.RRTMG.nc', 2, rh,thick,dens); 
[aod_ss3] = aer_opt_lw(direc,faer.SS003, 'opticsBands_SS.v3_5.RRTMG.nc', 3, rh,thick,dens);  
[aod_ss4] = aer_opt_lw(direc,faer.SS004, 'opticsBands_SS.v3_5.RRTMG.nc', 4, rh,thick,dens);  
[aod_ss5] = aer_opt_lw(direc,faer.SS005, 'opticsBands_SS.v3_5.RRTMG.nc', 5, rh,thick,dens);  

aod = cat(3,aod_du1,aod_du2,aod_du3,aod_du4,aod_du5,aod_ss1,aod_ss2,aod_ss3,aod_ss4,aod_ss5,...
           aod_bco,aod_bci,aod_oco,aod_oci,aod_sa);

aodt = sum(aod,3);
% save the data into netcdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- DEFINE THE FILE ----------------------------
  outfile = [direc,'IN_AER_RRTM_LW.nc'];
  ncid = netcdf.create(outfile,'CLOBBER');
  %ncid = netcdf.open('olr_nc_test.nc','WRITE');
  %nccreate('olr_nc_test.nc','bdps_olr');
  %-----------------------------DEFINE DIMENSION----------------------------
  dimidz = netcdf.defDim(ncid,'level',size(aodt,2));  
  dimidw = netcdf.defDim(ncid,'wavebands',size(aodt,1));
  %----------------------------DEFINE NEW VARIABLES-------------------------
  levvarid  = netcdf.defVar(ncid,'mdl_level','NC_INT',dimidz);
  banvarid  = netcdf.defVar(ncid,'wavebands','NC_INT',dimidw);
  varid1 = netcdf.defVar(ncid,'aod','NC_FLOAT',[dimidw dimidz]);
  %------------------------------DEFINE ATTRIBUTE---------------------------
  netcdf.putAtt(ncid,levvarid,'units','model_level.');
  netcdf.putAtt(ncid,banvarid,'units','wavebands.');
  netcdf.putAtt(ncid,varid1,'units','dimensionless.');
   %---------------------------GIVE VALUES TO VARIABLES-----------------------
  netcdf.endDef(ncid)
   
  netcdf.putVar(ncid,levvarid,1:72);
  netcdf.putVar(ncid,banvarid,1:16);
  netcdf.putVar(ncid,varid1,aodt);

  netcdf.close(ncid);

%{
% generate INPUT file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 aid = fopen('C:/Users/lwx/Desktop/IN_AER_RRTM1', 'w');
 
 fprintf(aid, '   71\n');
 % levls，AOD、SSA and g
 num = 1:71;
 for l = 1:71
     NLAY = 1;
     IAOD = 1;
     fprintf(aid, '   %2d    %d    %d    %d\n', NLAY, IAOD, IAOD, IAOD);
 % 3. 每层的 AOD
     fprintf(aid, '  %3d', num(l));
     for n = 1:16
            % 注意：MATLAB 索引从 1 开始，与 Python 从 0 开始不同
            % 原代码中 aodt[n,-1*l-1] 对应从底层到顶层的顺序
         fprintf(aid, '%7.4f', aodt(n, 73-l));  % 72-l 实现从底层到顶层的索引
     end
     fprintf(aid, '\n');
end
% 关闭文件
fclose(aid);
%}
