%% READ_TAMURA_DAILY
% Read and save daily ERA-interim data from T. Tamura and save to mat file for ROMS surface forcing file creation
%
% Created by B K. Galton-Fenzi 
% Adapted to Daily forcing by E. Cougnon at /ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily/ (May 2013)
% Poked around a bit by D Gwyther at /ds/projects/iomp/totten/ana/dgwyther/netcdf/forcing_creator/sbc/read_tamura_daily.m (March and Dec 2014)
% major changes by D Gwyther, Jan 2016
%%%%%%

lon_rho=ncread(grdname,'lon_rho')';
lat_rho=ncread(grdname,'lat_rho')';

NumYears = MaxYear-MinYear+1;

shfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,365*NumYears); %a little bit smaller
ssfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,365*NumYears); %due to not including
shflux = nan(365*NumYears,xmax-xmin+1,ymax-ymin+1);    %leap years.
ssflux = nan(365*NumYears,xmax-xmin+1,ymax-ymin+1);

%% Read in land mask
fid3=fopen('/ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily_latest/EASE_landmask_H.data','r');
landmaskNaN = reshape(fread(fid3,721*721*1,'float32=>double'),721,721);
landmaskNaN(landmaskNaN==0)=NaN;

ii = 1;
Month = ['jan';'feb';'mar';'apr';'may';'jun';'jul';'aug';'sep';'oct';'nov';'dec'];
DaysPerMonth=[31,28,31,30,31,30,31,31,30,31,30,31];
LeapYears = [1992:4:2040]; %leap years til 2040
%% Read in data:
for YearInd = MinYear:MaxYear; %1992:2007;
 for MonthInd = 1:12; %Keep track of current month

  %this loop corrects for leap years
  if any(YearInd == LeapYears)
   if MonthInd == 2 %correct feb no.of.days
    dms = 29;
   else
    dms = DaysPerMonth(MonthInd);
   end
  else
   dms = DaysPerMonth(MonthInd);
  end

    display(['Reading in GrADS data ' num2str(YearInd),' ',Month(MonthInd,:),' ...']);
    eval(['fid = fopen(''/ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily_latest/TSDM2hb_' num2str(YearInd),'_',Month(MonthInd,:),'.data'',''r'');']),
    if YearInd < 2002 %See readme file -- data from 1992 to 2002 -- 9 variables (ERA40, ERA-interim and NCEP2) while > 2002 only 6 variables (ERA-interim and NCEP2)
        tmp = reshape(fread(fid,721*721*dms*9,'float32=>double'),721,721,9,dms);
    else
        tmp = reshape(fread(fid,721*721*dms*6,'float32=>double'),721,721,6,dms); 
    end
    
    tmp = bsxfun(@times,tmp,landmaskNaN); % land mask
    
    if YearInd < 2002 % See readme file -- data from 1992 to 2002 -- the ERA-interim data are in variable 4(hf), 5(ip) and 6(sf)
        shfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,4,:)); %[W/m^2] 7:NCEP 4:ERA-interim
        ssfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,6,:)); %[Kg/m^2] 9:NCEP 6:ERA-interim
    else
        shfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,1,:)); %[W/m^2] 4:NCEP 1:ERA-interim
        ssfluxtmp(:,:,ii:ii+dms-1) = squeeze(tmp(xmin:xmax,ymin:ymax,3,:)); %[Kg/m^2] 6:NCEP 3:ERA-interim
    end
    
  ii = ii+dms;
 end
save YearInd.txt YearInd -ascii
end

DayNumber = ii-1;% save the day index to a more descriptive Varname (don't forget to remove 1 to account forfirst day being indexed)

if 0 % do inpainting
 for iii = 1:size(shfluxtmp,3);
    disp(['interpolating NaNs in fluxes...',num2str(iii/size(shfluxtmp,3)*100),'% done'])
    shflux(iii,:,:) = inpaint_nans(squeeze(shfluxtmp(:,:,iii)),2); % Heat: Watts/m^2
    ssflux(iii,:,:) = inpaint_nans(squeeze(ssfluxtmp(:,:,iii)),2); % (m/month)/(s/month)*(kg/m^2)  cm/day
 end
elseif 1
shflux = permute(shfluxtmp,[3 1 2]);
ssflux = permute(ssfluxtmp,[3 1 2]);
end

%%

disp('Saving Heat/Salt fluxes from Takeshi''s data')
save([RunName,'_Takeshi_subset_daily.mat'],'shfluxtmp','ssfluxtmp')
disp('Saved, now gridding to model grid')

fid2 = fopen('/ds/projects/iomp/obs/Tamura_air_sea_fluxes/daily/latlon.data');
B = fread(fid2,721*721*2,'float32=>double');
BB = reshape(B,721,721,2);

display('Making tidy for ROMS grid...');
AllLon = squeeze(BB(xmin:xmax,ymin:ymax,2));
AllLat = squeeze(BB(xmin:xmax,ymin:ymax,1));

ssfluxGrid = nan(DayNumber,size(lat_rho,1),size(lat_rho,2));
shfluxGrid = nan(DayNumber,size(lat_rho,1),size(lat_rho,2));

for j = 1:DayNumber;
    ssfluxGrid(j,:,:) = griddata(AllLon,AllLat,squeeze(ssflux(j,:,:)),lon_rho,lat_rho,'cubic');
    shfluxGrid(j,:,:) = griddata(AllLon,AllLat,squeeze(shflux(j,:,:)),lon_rho,lat_rho,'cubic');
    disp(['gridding ',num2str(j/DayNumber*100),' done.'])
end

disp('Saving gridded Heat/Salt fluxes')
save([RunName,'_air_sea_fluxes_daily.mat'],'shfluxGrid','ssfluxGrid','-v7.3')
disp('Saved.')
