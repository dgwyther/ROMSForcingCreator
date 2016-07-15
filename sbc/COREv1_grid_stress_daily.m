%%CORE_grid_stress_daily_to_monthly
% 
% script to load CORE data and convert to Monthly and interpolate to the model grid
%
% BGF wrote original script
% DEG major update (Dec 2015)
% DEG rewrite for TOK forcing files (Jul 2016)

% Load model grid:
lon_rho = ncread(grdname,'lon_rho')';
lat_rho = ncread(grdname,'lat_rho')';
angle = ncread(grdname,'angle')';
mask_rho = ncread(grdname,'mask_rho')';


% MinYear = 1992;
% MaxYear = 2007; %1993; % This is loaded from make_sbc.m
NumYears = MaxYear-MinYear+1;
SamRate = 4;% Sampling per day

%%% Load data
U_10_data = '/ds/projects/iomp/obs/COREv1_23JUN2016/CNYF/u_10_mod.clim.nc';
V_10_data = '/ds/projects/iomp/obs/COREv1_23JUN2016/CNYF/v_10_mod.clim.nc';

U_10_MOD = ncread(U_10_data,'U_10_MOD',[Imin_wind J_min_wind 1],[Imax_wind-Imin_wind+1 Jmax_wind-Jmin_wind+1 Inf]); %loads as lon,lat,time
V_10_MOD = ncread(V_10_data,'V_10_MOD',[Imin_wind J_min_wind 1],[Imax_wind-Imin_wind+1 Jmax_wind-Jmin_wind+1 Inf]); %loads as lon,lat,time
LONS = ncread(U_10_data,'LON');
LATS = ncread(V_10_data,'LAT');
[LONS,LATS]=meshgrid(LON(Imin_wind:Imax_wind),LAT(Jmin_wind:Jmax_wind));


%% convert vel -> stress
rhoAir= 1.3;
Cd = 1.4e-3;

taux = rhoAir*Cd.*abs(U_10_MOD).*U_10_MOD;
tauy = rhoAir*Cd.*abs(V_10_MOD).*V_10_MOD;


if 0 % daily->monthly data
% Makes monthly climatologies from daily climatologies: 
MDM = [31 28 31 30 31 30 31 31 30 31 30 31]; % matrix number of day in each month, without leap years
kk=1;
for ii = 1:12
        uwm_stress(ii,:,:) = nanmean(squeeze(taux(kk:kk+SamRate*MDM(ii)-1,:,:)));
        uws_stress(ii,:,:) = nanstd(squeeze(taux(kk:kk+SamRate*MDM(ii)-1,:,:)));
        vwm_stress(ii,:,:) = nanmean(squeeze(tauy(kk:kk+SamRate*MDM(ii)-1,:,:)));
        vws_stress(ii,:,:) = nanstd(squeeze(tauy(kk:kk+SamRate*MDM(ii)-1,:,:)));
        kk = kk+SamRate*MDM(ii);
end


AISuwm_stress=[];
AISvwm_stress=[];

 % Interpolate each months data to your grid
for i = 1:12;i,
    AISuwm_stress(i,:,:) = griddata(LONS,LATS,squeeze(uwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
    AISvwm_stress(i,:,:) = griddata(LONS,LATS,squeeze(vwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
end


end
if 1 % daily data

uwm_stress=nanmean(reshape(taux,Imax_wind-Imin_wind+1,Jmax_wind-Jmin_wind+1,SamRate,[]),3);
vwm_stress=nanmean(reshape(tauy,Imax_wind-Imin_wind+1,Jmax_wind-Jmin_wind+1,SamRate,[]),3);

for ii = 1:size(uwm_stress,3);
AISuwm_stress(ii,:,:) = griddata(LONS,LATS,squeeze(uwm_stress(:,:,ii)),lon_rho,lat_rho,'cubic');
AISvwm_stress(ii,:,:) = griddata(LONS,LATS,squeeze(vwm_stress(:,:,ii)),lon_rho,lat_rho,'cubic');
if ~rem(ii,round(size(uwm_stress,3)/10)), disp([num2str(ii/size(uwm_stress,3)*100),'% interpolated']), end
end
end


%% Rotates currents from model domain XI, ETA to north, south grid.
[p q r] = size(AISuwm_stress); %[time lon lat]
urot_stress = nan(p,q,r);
vrot_stress = nan(p,q,r);


for i = 1:q;
    for j = 1:r;
        
    uIn_stress = squeeze(AISuwm_stress(:,i,j));
    vIn_stress = squeeze(AISvwm_stress(:,i,j));
    
%[uOut_stress vOut_stress] = rotate_model_currents([uIn_stress vIn_stress],-angle(i,j));

uOut_stress = uIn_stress;

vOut_stress = vIn_stress;

urot_stress(:,i,j) = uOut_stress;
vrot_stress(:,i,j) = vOut_stress;
    end
end

u_stress_All = urot_stress;
v_stress_All = vrot_stress;
disp('Saving u_stress and v_stress')
% save u and v component with the model grid, and a vector rotation
nameval=['ustress_grid_model.mat'];
save(nameval,'u_stress_All','-v7.3');
nameval=['vstress_grid_model.mat'];
save(nameval,'v_stress_All','-v7.3');
disp('Saved.')

