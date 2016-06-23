% script to load CORE data and convert to Monthly and interpolate to the
% model grid with the katabatic mask --> for wind stress (sea surface)

% Load CORE data:
%addpath(genpath('/u/ameijers/matlab_scripts')) %script to read ncep data
%addpath(genpath('/v/hpclibrary/mertz/mdl/frc/COREv2_26JAN2010/CIAF/'))

% Load model grid:
ncload(grdname,'lon_rho','lat_rho','angle','mask_rho')


% MinYear = 1992;
% MaxYear = 2007; %1993; % This is loaded from make_sbc.m
NumYears = MaxYear-MinYear+1;
SamRate = 4;% Sampling per day
dpy=365; % Days per year: 366 for leap year (1992,1996,2000,2004,2008), 365 plain year.

% % choose a smaller region for the Totten
% Imin_core = 60;
% Imax_core = 80;
% Jmin_core = 10; % Now, this data is loaded from make_sbc.m
% Jmax_core = 18;



uwm_stress = nan(NumYears*12,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);
uws_stress = nan(NumYears*12,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);
vwm_stress = nan(NumYears*12,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);
vws_stress = nan(NumYears*12,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);

uwndall = nan(NumYears*dpy*SamRate,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);
vwndall = nan(NumYears*dpy*SamRate,Jmax_core-Jmin_core+1,Imax_core-Imin_core+1);
i = 1;
j = 1;

%%% For 'new' spinup procedure, loop 1992 forcing 16x.
for Year = MinYear.*ones(1,NumYears);
%for Year = MinYear:MaxYear
  Year,
  eval(['ncload /ds/projects/iomp/obs/COREv2_26JAN2010/CIAF/u_10.' num2str(Year) '.26JAN2010.nc U_10_MOD LON LAT']);
  eval(['ncload /ds/projects/iomp/obs/COREv2_26JAN2010/CIAF/v_10.' num2str(Year) '.26JAN2010.nc V_10_MOD LON LAT']);

  uwndall(i:i+length(U_10_MOD)-1,:,:) = squeeze(U_10_MOD(:,Jmin_core:Jmax_core,Imin_core:Imax_core));
  vwndall(j:j+length(V_10_MOD)-1,:,:) = squeeze(V_10_MOD(:,Jmin_core:Jmax_core,Imin_core:Imax_core));

  i = i+length(U_10_MOD); % I removed a '-1' here. Seemed to be overwriting data unneccessarily
  j = j+length(V_10_MOD); % It overwrote last data of each year, with first of next: DG. 28/7/2011
end

% transform wind velocity data (=CORE data, 10m above the sea surface) in
% stress wind (wind on sea surface)

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
    AISuwm_stress(i,:,:) = griddata(LON(Imin_core:Imax_core),LAT(Jmin_core:Jmax_core),squeeze(uwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
    AISvwm_stress(i,:,:) = griddata(LON(Imin_core:Imax_core),LAT(Jmin_core:Jmax_core),squeeze(vwm_stress(i,:,:)),lon_rho,lat_rho,'cubic');
end

% % % % multiply CORE's data by katamask
% % % 
% % % load /v/hpclibrary/mertz/ana/ecougnon/core_katamask__model_grid/katamask_x__model_grid.mat
% % % load /v/hpclibrary/mertz/ana/ecougnon/core_katamask__model_grid/katamask_y__model_grid.mat
% % % 
% % % [p q r] = size(AISuwm_stress);
% % % AISuwm_kata_stress = nan(p,q,r);
% % % AISvwm_kata_stress = nan(p,q,r);
% % % 
% % % for i = 1:q;
% % %     for j = 1:r;
% % %         
% % %     AISuwm_kata_stress(:,i,j) = AISuwm_stress(:,i,j).*katamask_x__model_grid(i,j);
% % %     AISvwm_kata_stress(:,i,j) = AISvwm_stress(:,i,j).*katamask_y__model_grid(i,j);
% % %     
% % %     end
% % % end

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
disp('Saving u_stress and v_stress')
% save u and v component with the model grid, and a vector rotation
save([WorkingDir,'ustress_grid_model'],'urot_stress')
save([WorkingDir,'vstress_grid_model'],'vrot_stress')
disp('Saved.')
% % % % % Rotates currents from model domain XI, ETA to north, south grid with
% % % % % CORE's data with katamask
% % % % 
% % % % [p q r] = size(AISuwm_kata_stress); %[time lon lat]
% % % % urot_kata_stress = nan(p,q,r);
% % % % vrot_kata_stress = nan(p,q,r);
% % % % 
% % % % for i = 1:q;
% % % %     for j = 1:r;
% % % %         
% % % %     uIn_kata_stress = squeeze(AISuwm_kata_stress(:,i,j));
% % % %     vIn_kata_stress = squeeze(AISvwm_kata_stress(:,i,j));
% % % %     
% % % % [uOut_kata_stress vOut_kata_stress] = rotate_model_currents([uIn_kata_stress vIn_kata_stress],-angle(i,j));
% % % % 
% % % % urot_kata_stress(:,i,j) = uOut_kata_stress;
% % % % vrot_kata_stress(:,i,j) = vOut_kata_stress;
% % % %     end
% % % % end
% % % % 
% % % % % save u and v component with the model grid and the katabatic mask, with a
% % % % % vector rotation
% % % % save ustress_grid_model_kata urot_kata_stress
% % % % save vstress_grid_model_kata vrot_kata_stress
% % % % 
% % % % figure,
% % % % pcolor(lon_rho,lat_rho,squeeze(urot_kata_stress(6,:,:))), 
% % % % shading flat,
% % % % title('u wind (ish)'),
% % % % colorbar,
% % % % figure,
% % % % pcolor(lon_rho,lat_rho,squeeze(vrot_kata_stress(6,:,:))), 
% % % % shading flat,
% % % % title('v wind (ish)'),
% % % % colorbar,

%{
%test the rotation direction
%figure, plot([0 AISuwm_kata(1,150,300)],[0 AISvwm_kata(1,150,300)],'bo-'); 
%hold on,
%plot([0 urot_kata(1,150,300)],[0 vrot_kata(1,150,300)],'ro-'); 

% test plot to AIS (before rotation) see the difference with and rot (next
% plot)
%figure,
%pcolor(lon_rho,lat_rho,squeeze(urot_stress(6,:,:))), 
%shading flat,
%title('u wind (m/s)'),
%colorbar,
%figure,
%pcolor(lon_rho,lat_rho,squeeze(vrot_stress(6,:,:))), 
%shading flat,
%title('v wind (m/s)'),
%colorbar,
%}

