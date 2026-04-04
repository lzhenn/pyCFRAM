function [aod1]=aer_opt_lw(direc,faer,name_opt,rd,rh,thick,dens)

%mdpn = [180,425,565,665,760,900,1030,1130,1285,1435,1640,1940,2165,2315,2490,2925];
%mdpl = 1./mdpn.*1e-2;
lwbs = 15:30;

q = faer;%get the 'name_m2' aersol
ds = [direc,name_opt];%The path for storing the optbands files
bext = ncread(ds,'bext');%Mass Extinction Coefficient

if rd < 1 %rd directly represents the effective radius.
    radius = ncread(ds,'radius');
    Lr = find(abs(radius-rd)==min(abs(radius-rd)));
    kext0 = bext(lwbs,:,Lr);  %unit: m2/kg
else  %easier to distinguish sulfate and others,rd represents the index of the array.
    kext0 = bext(lwbs,:,rd);  %at 1um
end
    %find out the nearest rh
rh0 = ncread(ds,'rh');
kext1 = zeros(length(lwbs),72);
for i = 1:72
    rh1 = abs(rh0-rh(i)); %levels of rh0 diff from rh[i]
    plc = find(rh1==min(rh1));
    kext1(:,i) = kext0(:,plc);
end

aod1 = q'.*dens'.*1e3.*thick'.*kext1.*1e-3;
aod1(isnan(aod1)) = 0;





