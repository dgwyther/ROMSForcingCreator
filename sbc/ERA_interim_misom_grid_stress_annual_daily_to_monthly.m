%% ERA_interim_misom_grid_stress_annual.m
%
% script to load CORE data and convert to Monthly and interpolate to the
% model grid with the katabatic mask --> for wind stress (sea surface)
%
% E Cougnon - Wrote original script for loading ERA data (July 2014)
% D Gwyther - adapted for 1-daily forcing, some corrections and alterations to fit in my BC creator framework (dec 2014)
%	may 2017	Updates to remove old functions
%
%%

addpath('/ds/projects/iomp/matlab_scripts')
% Load model grid (grdname from make_sbc.m):
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
mask_rho=ncread(grdname,'mask_rho')';
angle=ncread(grdname,'angle')';


NumYears = MaxYear-MinYear+1;
SamRate = 1;% Sampling per day -- 2daily for ERA interim

Imin_ERAi = Imin_wind; %lonmin
Imax_ERAi = Imax_wind;
Jmin_ERAi = Jmin_wind; %latmin
Jmax_ERAi = Jmax_wind;



disp('loading ERA-interim data')
u10=ncread('/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_1992_2014.1daily.u10.nc','u10',[Imin_ERAi Jmin_ERAi 1],[Imax_ERAi-Imin_ERAi+1 Jmax_ERAi-Jmin_ERAi+1 Inf]);
u10 = permute(u10,[3 2 1]);
v10=ncread('/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_1992_2014.1daily.v10.nc','v10',[Imin_ERAi Jmin_ERAi 1],[Imax_ERAi-Imin_ERAi+1 Jmax_ERAi-Jmin_ERAi+1 Inf]);
v10 = permute(v10,[3 2 1]);
longitude=double(ncread('/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_1992_2014.1daily.u10.nc','longitude',[Imin_ERAi],[Imax_ERAi-Imin_ERAi+1]));
latitude=double(ncread('/ds/projects/iomp/obs/ERA_Interim/ERA_Interim_1992_2014.1daily.u10.nc','latitude', [Jmin_ERAi],[Jmax_ERAi-Jmin_ERAi+1]));


uwndall=[];
vwndall=[];
u_stress_All = [];
v_stress_All = [];
u_index=1;
v_index=1;

LeapYears = [1992:4:2040]; %leap years til 2040
dpm    = [31,28,31,30,31,30,31,31,30,31,30,31];
dpm_ly = [31,29,31,30,31,30,31,31,30,31,30,31];
dpm_c    = [0,cumsum(dpm)];
dpm_ly_c = [0,cumsum(dpm_ly)];

for YearInd = MinYear:MaxYear; 
disp(num2str(YearInd))
clear uwnd vwnd uw_stress vw_stress signu signv
 if any(YearInd == LeapYears)
          % allocate space for matrix for daily data for one year
  uwnd = nan(366*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vwnd = nan(366*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
%  uw_stress = nan(366,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
%  vw_stress = nan(366,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);          
 else
          % allocate space for matrix for daily data for one year
  uwnd = nan(365*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
  vwnd = nan(365*SamRate,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
%  uw_stress = nan(365,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
%  vw_stress = nan(365,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
 end

% uw_stress = nan(12,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);
% vw_stress = nan(12,Jmax_ERAi-Jmin_ERAi+1,Imax_ERAi-Imin_ERAi+1);

 uwnd(:,:,:) = squeeze(u10(u_index:u_index+size(uwnd,1)-1,:,:));
 vwnd(:,:,:) = squeeze(v10(v_index:v_index+size(vwnd,1)-1,:,:));
disp('upscaling daily -> monthly')
for ii=1:12;
 if any(YearInd == LeapYears)
  uwnd2(ii,:,:) = nanmean(uwnd(dpm_ly_c(ii)+1:dpm_ly_c(ii+1),:,:),1);
  vwnd2(ii,:,:) = nanmean(vwnd(dpm_ly_c(ii)+1:dpm_ly_c(ii+1),:,:),1);
 else
  uwnd2(ii,:,:) = nanmean(uwnd(dpm_c(ii)+1:dpm_c(ii+1),:,:),1);
  vwnd2(ii,:,:) = nanmean(vwnd(dpm_c(ii)+1:dpm_c(ii+1),:,:),1);
 end
end
   signu = sign(uwnd2);
   signv = sign(vwnd2);

   rhoAir = 1.3;
   Cd = 1.4e-3;

   taux = rhoAir*Cd.*uwnd2.^2.*signu;
   tauy = rhoAir*Cd.*vwnd2.^2.*signv;

   % Want daily data, so far data are 2 daily -- sample one every 2 data to get the data
   % taken from 12.00
if 1 %daily data at 12:00 only
       uw_stress=taux;
       vw_stress=tauy;
elseif 0 % daily data at 00:00 and 12:00
   k=1;
   for i=1:2:length(uwnd)-1;
       uw_stress(k,:,:) = taux(i,:,:);
       vw_stress(k,:,:) = tauy(i,:,:);	
       k=k+1;
   end
end

   AISuw_stress=[];
   AISvw_stress=[];

   % Interpolate each daily data to ROMS grid
disp(['Interpolating daily data to ROMS grid for ',num2str(YearInd)])
   for ii = 1:size(uw_stress,1);
       AISuw_stress(ii,:,:) = griddata(longitude,latitude,squeeze(uw_stress(ii,:,:)),lon_rho,lat_rho,'cubic');
       AISvw_stress(ii,:,:) = griddata(longitude,latitude,squeeze(vw_stress(ii,:,:)),lon_rho,lat_rho,'cubic');
   end

   % %% Rotates currents from model domain XI, ETA to north, south grid.
   % [p q r] = size(AISuw_stress); %[time lon lat]
   % urot_stress = nan(p,q,r);
   % vrot_stress = nan(p,q,r);
   % 
   % for i = 1:q;
   %     for j = 1:r;
   %         
   %     uIn_stress = squeeze(AISuw_stress(:,i,j));
   %     vIn_stress = squeeze(AISvw_stress(:,i,j));
   %      
   %     [uOut_stress vOut_stress] = rotate_model_currents([uIn_stress vIn_stress],-angle(i,j));
   % 
   %     uOut_stress = uIn_stress;
   %     vOut_stress = vIn_stress;
   % 
   %     urot_stress(:,i,j) = uOut_stress;
   %     vrot_stress(:,i,j) = vOut_stress;
   %     end
   % end
   u_stress=AISuw_stress;
   v_stress=AISvw_stress;
   u_stress_All(u_index:u_index+12-1,:,:) = AISuw_stress;
   v_stress_All(v_index:v_index+12-1,:,:) = AISvw_stress;


   disp(['Saving u_stress and v_stress for ' num2str(YearInd) ''])
   % save u and v component with the model grid, and a vector rotation
   nameval=['ustress_grid_model_' num2str(YearInd) '_monthly.mat'];
   save(nameval,'u_stress','-v7.3'); 
   nameval=['vstress_grid_model_' num2str(YearInd) '_monthly.mat'];
   save(nameval,'v_stress','-v7.3');
   disp(['' num2str(YearInd) ' Saved.']),

   u_index = u_index+12;
   v_index = v_index+12;
end

 disp(['Saving u_stress and v_stress for ' num2str(MinYear) ' to ' num2str(MaxYear) ' where daily is upscaled to monthly data'])
   % save u and v component with the model grid, and a vector rotation
   nameval=['ustress_grid_model_',num2str(MinYear),'_',num2str(MaxYear),'_monthly.mat'];
   save(nameval,'u_stress_All','-v7.3');
   nameval=['vstress_grid_model_' num2str(MinYear),'_',num2str(MaxYear),'_monthly.mat'];
   save(nameval,'v_stress_All','-v7.3');

