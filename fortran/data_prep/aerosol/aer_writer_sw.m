%%%%%%%Call the aer_opt file to compute aerosol-related properties,
%%%%%%%and store the computation results in a netCDF file.
function aer_writer_sw(direc,faer,rh,delp,dens)

%mdpn = [180,425,565,665,760,900,1030,1130,1285,1435,1640,1940,2165,2315,2490,2925];
%mdpl = 1./mdpn.*1e-2;

thick = delp./(dens.*9.8);

[aod_du1,ssa_du1,g_du1] = aer_opt_sw(direc,faer.DU001,'opticsBands_DU.v15_3.RRTMG.nc',1,rh,thick,dens);
[aod_du2,ssa_du2,g_du2] = aer_opt_sw(direc,faer.DU002,'opticsBands_DU.v15_3.RRTMG.nc',2,rh,thick,dens);
[aod_du3,ssa_du3,g_du3] = aer_opt_sw(direc,faer.DU003,'opticsBands_DU.v15_3.RRTMG.nc',3,rh,thick,dens);
[aod_du4,ssa_du4,g_du4] = aer_opt_sw(direc,faer.DU004,'opticsBands_DU.v15_3.RRTMG.nc',4,rh,thick,dens);
[aod_du5,ssa_du5,g_du5] = aer_opt_sw(direc,faer.DU005,'opticsBands_DU.v15_3.RRTMG.nc',5,rh,thick,dens);
[aod_sa,ssa_sa,g_sa]    = aer_opt_sw(direc,faer.SO4, 'opticsBands_SU.v1_3.RRTMG.nc', 0.16*1e-6, rh,thick,dens);
[aod_bco,ssa_bco,g_bco] = aer_opt_sw(direc,faer.BCPHOBIC, 'opticsBands_BC.v1_3.RRTMG.nc', 1, rh,thick,dens);
[aod_bci,ssa_bci,g_bci] = aer_opt_sw(direc,faer.BCPHILIC, 'opticsBands_BC.v1_3.RRTMG.nc', 2, rh,thick,dens);
[aod_oco,ssa_oco,g_oco] = aer_opt_sw(direc,faer.OCPHOBIC, 'opticsBands_OC.v1_3.RRTMG.nc', 1, rh,thick,dens);
[aod_oci,ssa_oci,g_oci] = aer_opt_sw(direc,faer.OCPHILIC, 'opticsBands_OC.v1_3.RRTMG.nc', 2, rh,thick,dens);
[aod_ss1,ssa_ss1,g_ss1] = aer_opt_sw(direc,faer.SS001, 'opticsBands_SS.v3_5.RRTMG.nc', 1, rh,thick,dens); 
[aod_ss2,ssa_ss2,g_ss2] = aer_opt_sw(direc,faer.SS002, 'opticsBands_SS.v3_5.RRTMG.nc', 2, rh,thick,dens); 
[aod_ss3,ssa_ss3,g_ss3] = aer_opt_sw(direc,faer.SS003, 'opticsBands_SS.v3_5.RRTMG.nc', 3, rh,thick,dens);  
[aod_ss4,ssa_ss4,g_ss4] = aer_opt_sw(direc,faer.SS004, 'opticsBands_SS.v3_5.RRTMG.nc', 4, rh,thick,dens);  
[aod_ss5,ssa_ss5,g_ss5] = aer_opt_sw(direc,faer.SS005, 'opticsBands_SS.v3_5.RRTMG.nc', 5, rh,thick,dens);  

aod = cat(3,aod_du1,aod_du2,aod_du3,aod_du4,aod_du5,aod_ss1,aod_ss2,aod_ss3,aod_ss4,aod_ss5,...
           aod_bco,aod_bci,aod_oco,aod_oci,aod_sa);
ssa = cat(3,ssa_du1,ssa_du2,ssa_du3,ssa_du4,ssa_du5,ssa_ss1,ssa_ss2,ssa_ss3,ssa_ss4,ssa_ss5,...
           ssa_bco,ssa_bci,ssa_oco,ssa_oci,ssa_sa);
g = cat(3,g_du1,g_du2,g_du3,g_du4,g_du5,g_ss1,g_ss2,g_ss3,g_ss4,g_ss5, g_bco,g_bci,g_oco,g_oci,g_sa);

aodt = sum(aod,3);
ssat = sum(ssa.*aod,3)./sum(aod,3);
gt = sum(g.*ssa.*aod,3)./sum(ssa.*aod,3);
% save the data into netcdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- DEFINE THE FILE ----------------------------
  outfile = [direc,'IN_AER_RRTM_SW.nc'];
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
  varid2 = netcdf.defVar(ncid,'ssa','NC_FLOAT',[dimidw dimidz]);
  varid3 = netcdf.defVar(ncid,'g','NC_FLOAT',[dimidw dimidz]);
  %------------------------------DEFINE ATTRIBUTE---------------------------
  netcdf.putAtt(ncid,levvarid,'units','model_level.');
  netcdf.putAtt(ncid,banvarid,'units','wavebands');
  netcdf.putAtt(ncid,varid1,'units','dimensionless');
  netcdf.putAtt(ncid,varid2,'units','dimensionless');
  netcdf.putAtt(ncid,varid3,'units','dimensionless');
   %---------------------------GIVE VALUES TO VARIABLES-----------------------
  netcdf.endDef(ncid)
   
  netcdf.putVar(ncid,levvarid,1:72);
  netcdf.putVar(ncid,banvarid,1:14);
  netcdf.putVar(ncid,varid1,aodt);
  netcdf.putVar(ncid,varid2,ssat);
  netcdf.putVar(ncid,varid3,gt);

  netcdf.close(ncid);

%{
% generate INPUT file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 aid = fopen('C:/Users/lwx/Desktop/IN_AER_RRTM2', 'w');
 
 fprintf(aid, '   71\n');
 % levls，AOD、SSA and g
 num = 1:71;
 for l = 1:71
     NLAY = 1;
     IAOD = 1;
     fprintf(aid, '   %2d    %d    %d    %d\n', NLAY, IAOD, IAOD, IAOD);
 % 3. 每层的 AOD
     fprintf(aid, '  %3d', num(l));
     for n = 1:14
            % 注意：MATLAB 索引从 1 开始，与 Python 从 0 开始不同
            % 原代码中 aodt[n,-1*l-1] 对应从底层到顶层的顺序
         fprintf(aid, '%7.4f', aodt(n, 73-l));  % 72-l 实现从底层到顶层的索引
     end
     fprintf(aid, '\n');
     % 4. SSA
     if aodt(10, 72-l) < 0.00005  % 第10个波段
        for n = 1:14
            fprintf(aid, '%5.2f', 1);
        end
     else
        for n = 1:14
            fprintf(aid, '%5.2f', ssat(n, 73-l));
        end
     end
     fprintf(aid, '\n');
     % 5. 不对称因子 g
     for n = 1:14
         fprintf(aid, '%5.2f', gt(n, 73-l));
     end
     fprintf(aid, '\n');
end
% 关闭文件
fclose(aid);
%}
