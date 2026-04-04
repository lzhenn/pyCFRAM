function [aod1,ssa1,g1]=aer_opt_sw(direc,faer,name_opt,rd,rh,thick,dens)

%mdpn = [2925,3625,4325,4900,5650,6925,7875,10450,14425,19325,25825,33500,44000,1710];
%mdpl = 1./mdpn.*1e-2;
swbs = 1:14;

q = faer;%get the 'name_m2' aersol
ds = [direc,name_opt];%The path for storing the optbands files  
bext = ncread(ds,'bext');%Mass Extinction Coefficient
bsca = ncread(ds,'bsca');%Mass Scattering Coefficient
qext = ncread(ds,'qext');%Extinction Efficiency
qsca = ncread(ds,'qsca');%Scattering Efficiency
g = ncread(ds,'g');%Asymmetry Factor

if rd < 1%rd directly represents the effective radius.
    radius = ncread(ds,'radius');
    Lr = find(abs(radius-rd)==min(abs(radius-rd)));
    kext0 = bext(swbs,:,Lr);  %unit: m2/kg
    ssa0 = qsca(swbs,:,Lr)./qext(swbs,:,Lr); %Single Scattering Albedo;Or bsca./bext
    g0 = g(swbs,:,Lr);
else  %easier to distinguish sulfate and others,rd represents the index of the array.
    kext0 = bext(swbs,:,rd);  %at 1um
    ssa0 = qsca(swbs,:,rd)./qext(swbs,:,rd); %Or bsca./bext
    g0 = g(swbs,:,rd);
end
    %find out the nearest rh
rh0 = ncread(ds,'rh');
kext1 = zeros(length(swbs),72);
ssa1 = zeros(length(swbs),72);
g1 = zeros(length(swbs),72);
for i = 1:72
    rh1 = abs(rh0-rh(i)); %levels of rh0 diff from rh[i]
    plc = find(rh1==min(rh1));
    kext1(:,i) = kext0(:,plc);
    ssa1(:,i) = ssa0(:,plc);
    g1(:,i) = g0(:,plc);
end

aod1 = q'.*dens'.*1e3.*thick'.*kext1.*1e-3;
aod1(isnan(aod1)) = 0;
ssa1(isnan(ssa1)) = 0;
g1(isnan(g1)) = 0;




