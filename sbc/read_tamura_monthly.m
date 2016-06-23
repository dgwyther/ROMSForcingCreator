%% READ_TAMURA_MONTHLY

% Created by B. K. Galton-Fenzi
% Updated by D. Gwyther (Dec 2015) for use with a high res Amery run

lon_rho=ncread(grdname,'lon_rho')';
lat_rho=ncread(grdname,'lat_rho')';

NumYears = MaxYear-MinYear+1;

shfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,12*NumYears);
ssfluxtmp = nan(xmax-xmin+1,ymax-ymin+1,12*NumYears);
shflux = nan(12*NumYears,xmax-xmin+1,ymax-ymin+1);
ssflux = nan(12*NumYears,xmax-xmin+1,ymax-ymin+1);

%% Read in land mask
fid3=fopen('/ds/projects/iomp/obs/Tamura_air_sea_fluxes/landmask.data','r');
landmaskNaN = reshape(fread(fid3,721*721*1,'float32=>double'),721,721); 
landmaskNaN(landmaskNaN==0)=NaN;


ii = 1;

%% Read in data:
for YearInd = MinYear:MaxYear; 
    
    display(['Reading in GrADS data ' num2str(YearInd) ' ...']);
    eval(['fid = fopen(''/ds/projects/iomp/obs/Tamura_air_sea_fluxes/heat_salt_flux_' num2str(YearInd) '.data'',''r'');']),
    if YearInd < 2002
        tmp = reshape(fread(fid,721*721*12*9,'float32=>double'),721,721,9,12);
    else
        tmp = reshape(fread(fid,721*721*12*6,'float32=>double'),721,721,6,12); 
    end
    
    tmp = bsxfun(@times,tmp,landmaskNaN); % land mask
 
    if YearInd < 2002
        shfluxtmp(:,:,ii:ii+11) = squeeze(tmp(xmin:xmax,ymin:ymax,7,:)); % W/m^2
        ssfluxtmp(:,:,ii:ii+11) = squeeze(tmp(xmin:xmax,ymin:ymax,9,:)); % Kg/m^2
    else
        shfluxtmp(:,:,ii:ii+11) = squeeze(tmp(xmin:xmax,ymin:ymax,4,:)); % W/m^2
        ssfluxtmp(:,:,ii:ii+11) = squeeze(tmp(xmin:xmax,ymin:ymax,6,:)); % Kg/m^2
    end
    
    ii = ii+12;
end

if 0 %do inpainting
 for ii = 1:12*NumYears;
    ii,
    shflux(ii,:,:) = inpaint_nans(squeeze(shfluxtmp(:,:,i)),2); % Heat: Watts/m^2
    ssflux(ii,:,:) = inpaint_nans(squeeze(ssfluxtmp(:,:,i)),2); % (m/month)/(s/month)*(kg/m^2)  cm/day
 end
elseif 1
shflux = permute(shfluxtmp,[3 1 2]);
ssflux = permute(ssfluxtmp,[3 1 2]);
end

disp('Saving Heat/Salt fluxes from Takeshi''s data')
save([RunName,'_Takeshi_monthly_subset.mat'],'shfluxtmp','ssfluxtmp')
disp('Saved, now gridding to model grid')

%%
fid2 = fopen('/ds/projects/iomp/obs/Tamura_air_sea_fluxes/latlon.data');
B = fread(fid2,721*721*2,'float32=>double');
BB = reshape(B,721,721,2);

display('Making tidy ...');
AllLon = squeeze(BB(xmin:xmax,ymin:ymax,2));
AllLat = squeeze(BB(xmin:xmax,ymin:ymax,1));

ssfluxGrid = nan(12*NumYears,size(lat_rho,1),size(lat_rho,2));
shfluxGrid = nan(12*NumYears,size(lat_rho,1),size(lat_rho,2));

for jj = 1:12*NumYears;
    ssfluxGrid(jj,:,:) = griddata(AllLon,AllLat,squeeze(ssflux(jj,:,:)),lon_rho,lat_rho,'cubic');
    shfluxGrid(jj,:,:) = griddata(AllLon,AllLat,squeeze(shflux(jj,:,:)),lon_rho,lat_rho,'cubic');
    
end

%IceMask = zice; IceMask(zice < 0) = 1;
%LandMask = mask_rho;
%WaterMask = zice; WaterMask(zice < 0) = 0; WaterMask(zice == 0) = 1;WaterMask = WaterMask.*LandMask;
%WaterMaskNaN = WaterMask; WaterMaskNaN(WaterMask == 0 ) = NaN;
disp('Saving Heat/Salt fluxes')
save([RunName,'_air_sea_fluxes_monthly.mat'],'shfluxGrid','ssfluxGrid','-v7.3')
%!cp misom015_NCEP.mat /ds/projects/mertz/ana/bkgalton/mer015_IAF/sbc/.
disp('Saved.')
