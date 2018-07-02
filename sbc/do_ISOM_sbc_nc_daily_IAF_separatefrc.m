%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title =([RunName,'surface heat/salt fluxes and wind stress']); % Make title from model name variable

fluxes_file=matfile([RunName,'_air_sea_fluxes_daily.mat']);
ustress_file=matfile('ustress_grid_model.mat');
vstress_file=matfile('vstress_grid_model.mat');
frcname_orig=frcname;

force_id = {'su','sv','shf','swf'}; %load and interp ecco2 data

for tt = 1:4
disp('loading file')
clear swf shf su sv
if tt==1, su = ustress_file.u_stress_All; 
elseif tt==2, sv = vstress_file.v_stress_All;
elseif tt==3, shf = fluxes_file.shfluxGrid;
elseif tt==4, swf = fluxes_file.ssfluxGrid;
end


frcname = [frcname_orig(1:end-6),force_id{tt},'_sbc.nc'];
disp([' Creating the file : ',frcname,' at ' datestr(now)])

time_big = 1:365*(MaxYear-MinYear+1); %total length of forcing days
cycle = time_big(end); %last day after which cycle back to first day
%time_short = time_big((tt-MinYear)*365+1:(tt-MinYear)*365+1+364); %subset of total forcing over which this file forces
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;

id = netcdf.create(frcname, 'clobber');
%
%  Create dimensions
%
    disp('Create dimension')

xi_u_dim = netcdf.defDim(id, 'xi_u', L);
eta_u_dim = netcdf.defDim(id, 'eta_u', Mp);
xi_v_dim = netcdf.defDim(id, 'xi_v', Lp);
eta_v_dim = netcdf.defDim(id, 'eta_v', M);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', Lp);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', Mp);
xi_psi_dim = netcdf.defDim(id, 'xi_psi', L);
eta_psi_dim = netcdf.defDim(id, 'eta_psi', M);
sms_time_dim = netcdf.defDim(id, 'sms_time', length(time_big));
shf_time_dim = netcdf.defDim(id, 'shf_time', length(time_big));
swf_time_dim = netcdf.defDim(id, 'swf_time', length(time_big));
%
%  Create variables and attributes
%
    disp('Create variable attribute')
if tt==1 | tt==2
sms_time_id = netcdf.defVar(id, 'sms_time', 'double', sms_time_dim);
netcdf.putAtt(id, sms_time_id, 'long_name', 'surface momentum stress time');
netcdf.putAtt(id, sms_time_id, 'units', 'day');
netcdf.putAtt(id, sms_time_id, 'cycle_length', cycle);
elseif tt==3
shf_time_id = netcdf.defVar(id, 'shf_time', 'double', shf_time_dim);
netcdf.putAtt(id, shf_time_id, 'long_name', 'surface heat flux time');
netcdf.putAtt(id, shf_time_id, 'units', 'day');
netcdf.putAtt(id, shf_time_id, 'cycle_length', cycle);
elseif tt==4
swf_time_id = netcdf.defVar(id, 'swf_time', 'double', swf_time_dim);
netcdf.putAtt(id, swf_time_id, 'long_name', 'surface freshwater flux time');
netcdf.putAtt(id, swf_time_id, 'units', 'day');
netcdf.putAtt(id, swf_time_id, 'cycle_length', cycle);
end
if tt==1
sustr_id = netcdf.defVar(id, 'sustr', 'double', [xi_u_dim, eta_u_dim sms_time_dim]);
netcdf.putAtt(id, sustr_id, 'long_name', 'surface u-momentum stress');
netcdf.putAtt(id, sustr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, sustr_id, 'missing_value', -1e34);
netcdf.putAtt(id, sustr_id, '_FillValue', -1e34);
elseif tt==2
svstr_id = netcdf.defVar(id, 'svstr', 'double', [xi_v_dim, eta_v_dim sms_time_dim]);
netcdf.putAtt(id, svstr_id, 'long_name', 'surface v-momentum stress');
netcdf.putAtt(id, svstr_id, 'units', 'Newton meter-2');
netcdf.putAtt(id, svstr_id, 'missing_value', -1e34);
netcdf.putAtt(id, svstr_id, '_FillValue', -1e34);
elseif tt==3
shflux_id = netcdf.defVar(id, 'shflux', 'double', [xi_rho_dim, eta_rho_dim shf_time_dim]);
netcdf.putAtt(id, shflux_id, 'long_name', 'surface net heat flux');
netcdf.putAtt(id, shflux_id, 'units', 'Watts meter-2');
netcdf.putAtt(id, shflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, shflux_id, '_FillValue', -1e34);
elseif tt==4
swflux_id = netcdf.defVar(id, 'swflux', 'double', [xi_rho_dim, eta_rho_dim swf_time_dim]);
netcdf.putAtt(id, swflux_id, 'long_name', 'surface freshwater flux (E-P)');
netcdf.putAtt(id, swflux_id, 'units', 'centimetre day-1');
netcdf.putAtt(id, swflux_id, 'missing_value', -1e34);
netcdf.putAtt(id, swflux_id, '_FillValue', -1e34);
netcdf.putAtt(id, swflux_id, 'positive', 'net evaporation');
netcdf.putAtt(id, swflux_id, 'negative', 'net precipitation');
end
%
% Create global attributes
%
    disp('Create global attribute')
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', title);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'clim_file', frcname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');
netcdf.endDef(id);


    disp('Write time variable')

swft = time_big;  %data at each day of cycle
shft = time_big; 

if tt==1 | tt==2,netcdf.putVar(id, sms_time_id, swft);
elseif tt==3, netcdf.putVar(id, shf_time_id, shft);
elseif tt==4, netcdf.putVar(id, swf_time_id, swft); end

    disp('inpainting nans in wind fields...')
addpath('/ds/projects/iomp/matlab_scripts')
for j = 1:length(time_big);
if tt==3,if any(any(isnan(shf(j,:,:)))), shf(j,:,:)= inpaint_nans(squeeze( shf(j,:,:)),1); end, end
if tt==4,if any(any(isnan(swf(j,:,:)))), swf(j,:,:)= inpaint_nans(squeeze( swf(j,:,:)),1); end, end
if tt==1,if any(any(isnan(su(j,:,:)))), su(j,:,:) = inpaint_nans(squeeze(su(j,:,:)),1); end, end
if tt==2,if any(any(isnan(sv(j,:,:)))), sv(j,:,:) = inpaint_nans(squeeze(sv(j,:,:)),1); end, end
if ~rem(j,length(time_big)/20)
disp([num2str(j/length(time_big)*100),'%'])
end
end
rmpath('/ds/projects/iomp/matlab_scripts')


    disp('correct for large summer heat flux')

% Correct Heat flux for large summer values (as per NCEP2 heat flux data)
if tt==3,shf(shf > 0) = shf(shf > 0)*0.5; end

refSalt = 34.4; % Reference salinity


if tt==3,
temporary=permute(shf,[3 2 1]); 
pause(15), ('saving shflux to netcdf')
netcdf.putVar(id, shflux_id, temporary);% (W/m^2)
elseif tt==4,
temporary=permute(swf./refSalt*100,[3 2 1]);
pause(15)
netcdf.putVar(id, swflux_id, temporary); 
end% kg/m^2->cm day^-1

    disp('Write wind variables')

if tt==1,
%temporary=nan(size(su,3),size(su,2),size(su,1)-1);
temporary=permute((su(:,:,1:end-1)+su(:,:,2:end))./2,[3 2 1]); 
pause(15)
netcdf.putVar(id, sustr_id, temporary);
elseif tt==2,
%temporary=nan(size(sv,3),size(sv,2)-1,size(sv,1));
temporary=permute((sv(:,1:end-1,:)+sv(:,2:end,:))./2,[3 2 1]); 
pause(15)
netcdf.putVar(id, svstr_id, temporary); 
end
clear temporary

netcdf.close(id);

    disp('DONE')

end
