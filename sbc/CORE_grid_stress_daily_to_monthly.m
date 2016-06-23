%%CORE_grid_stress_daily_to_monthly
% 
% script to load CORE data and convert to Monthly and interpolate to the model grid
%
% BGF wrote original script
% DEG major update (Dec 2015)

% Load model grid:
lon_rho = ncread(grdname,'lon_rho')';
lat_rho = ncread(grdname,'lat_rho')';
angle = ncread(grdname,'angle')';
mask_rho = ncread(grdname,'mask_rho')';


% MinYear = 1992;
% MaxYear = 2007; %1993; % This is loaded from make_sbc.m
NumYears = MaxYear-MinYear+1;
SamRate = 4;% Sampling per day
dpy=365; % Days per year: 366 for leap year (1992,1996,2000,2004,2008), 365 plain year.

uwm_stress = nan(NumYears*12,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);
uws_stress = nan(NumYears*12,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);
vwm_stress = nan(NumYears*12,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);
vws_stress = nan(NumYears*12,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);

uwndall = nan(NumYears*dpy*SamRate,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);
vwndall = nan(NumYears*dpy*SamRate,Jmax_wind-Jmin_wind+1,Imax_wind-Imin_wind+1);
i = 1;
j = 1;

%%% For 'new' spinup procedure, loop 1992 forcing 16x.
%for Year = MinYear.*ones(1,NumYears);
for Year = MinYear:MaxYear
  Year,
eval(['filenameU =[''/ds/projects/iomp/obs/COREv2_26JAN2010/CIAF/u_10.'',num2str(Year),''.26JAN2010.nc''];']);
eval(['filenameV =[''/ds/projects/iomp/obs/COREv2_26JAN2010/CIAF/v_10.'',num2str(Year),''.26JAN2010.nc''];']);
U_10_MOD = permute(ncread(filenameU,'U_10_MOD'),[3 2 1]);
V_10_MOD = permute(ncread(filenameV,'V_10_MOD'),[3 2 1]);
LON = ncread(filenameU,'LON')';
LAT = ncread(filenameU,'LAT')';

  uwndall(i:i+length(U_10_MOD)-1,:,:) = squeeze(U_10_MOD(:,Jmin_wind:Jmax_wind,Imin_wind:Imax_wind));
  vwndall(j:j+length(V_10_MOD)-1,:,:) = squeeze(V_10_MOD(:,Jmin_wind:Jmax_wind,Imin_wind:Imax_wind));

  i = i+length(U_10_MOD); % I removed a '-1' here. Seemed to be overwriting data unneccessarily
  j = j+length(V_10_MOD); % It overwrote last data of each year, with first of next: DG. 28/7/2011
end

%make CORE grid
[LONS,LATS]=meshgrid(LON(Imin_wind:Imax_wind),LAT(Jmin_wind:Jmax_wind));

% transform wind velocity data (=CORE data, 10m above the sea surface) in
% stress wind (wind on sea surface)
%%
signu = sign(uwndall);
signv = sign(vwndall);

rhoAir = 1.3;
Cd = 1.4e-3;

taux = rhoAir*Cd.*uwndall.^2.*signu;
tauy = rhoAir*Cd.*vwndall.^2.*signv;


% Makes monthly climatologies from daily climatologies: 
MDM = [31 28 31 30 31 30 31 31 30 31 30 31]; % matrix number of day in each month, without leap years
k = 1; i = 1;
for i = 1:12:12*NumYears
    ii = i;
    for j = 1:12
        uwm_stress(ii,:,:) = nanmean(squeeze(taux(k:k+SamRate*MDM(j)-1,:,:)));
        uws_stress(ii,:,:) = nanstd(squeeze(taux(k:k+SamRate*MDM(j)-1,:,:)));
        vwm_stress(ii,:,:) = nanmean(squeeze(tauy(k:k+SamRate*MDM(j)-1,:,:)));
        vws_stress(ii,:,:) = nanstd(squeeze(tauy(k:k+SamRate*MDM(j)-1,:,:)));
        
        ii = ii+1;
        k = k+SamRate*MDM(j);
    end
end

AISuwm_stress=[];
AISvwm_stress=[];

 % Interpolate each months data to your grid
for i = 1:NumYears*12;i,
    AISuwm_stress(i,:,:) = griddata(LONS,LATS,squeeze(uwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
    AISvwm_stress(i,:,:) = griddata(LONS,LATS,squeeze(vwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
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

