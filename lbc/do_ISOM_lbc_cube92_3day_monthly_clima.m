% Load ECCO data
load(['cube92_iaf_salt_',RunName,'.mat'])
load(['cube92_iaf_theta_',RunName,'.mat'])
load(['cube92_iaf_uvel_',RunName,'.mat'])
load(['cube92_iaf_vvel_',RunName,'.mat'])

theta(theta < -10) = NaN;
salt(salt < 0) = NaN;

%% Make climatology
disp('making climatology')
theta = squeeze(nanmean(reshape(theta,[12,MaxYear-MinYear+1,size(theta,2),size(theta,3),size(theta,4)]),2));
salt = squeeze(nanmean(reshape(salt,[12,MaxYear-MinYear+1,size(salt,2),size(salt,3),size(salt,4)]),2));
uvel = squeeze(nanmean(reshape(uvel,[12,MaxYear-MinYear+1,size(uvel,2),size(uvel,3),size(uvel,4)]),2));
vvel = squeeze(nanmean(reshape(vvel,[12,MaxYear-MinYear+1,size(vvel,2),size(vvel,3),size(vvel,4)]),2));


%Cube92_3day_DateVec=datestr(datenum(MinYear,1,2):3:datenum(MaxYear+1,1,1),'yyyymmdd');
%CheckDayVec=str2num(datestr(datenum(1993,1,1):1:datenum(1993,12,31),'mmdd')); %non-leap yr
%YearLength=length(CheckDayVec);
%theta_full=theta; salt_full=salt; uvel_full=uvel; vvel_full=vvel; %make backups
%theta = nan(YearLength,size(theta,2),size(theta,3),size(theta,4));
%salt = nan(YearLength,size(salt,2),size(salt,3),size(salt,4));
%uvel = nan(YearLength,size(uvel,2),size(uvel,3),size(uvel,4));
%vvel = nan(YearLength,size(vvel,2),size(vvel,3),size(vvel,4));
%
%for tt=1:YearLength %1:num days in yr (no leap yr in cube92)
%find_ind = find(str2num(Cube92_3day_DateVec(:,5:end))==CheckDayVec(tt));
%theta(tt,:,:,:) = squeeze(nanmean(theta_full(find_ind,:,:,:),1));
%salt(tt,:,:,:) = squeeze(nanmean(salt_full(find_ind,:,:,:),1));
%uvel(tt,:,:,:) = squeeze(nanmean(uvel_full(find_ind,:,:,:),1));
%vvel(tt,:,:,:) = squeeze(nanmean(vvel_full(find_ind,:,:,:),1));
%if ~rem(tt,round(YearLength/10))
%disp([num2str(tt/YearLength*100) '% at ' datestr(now)])
%end
%end

%% remove a day to get no. days to 364......
%theta(365,:,:,:)=[]; 
%salt(365,:,:,:)=[];
%uvel(365,:,:,:)=[];
%vvel(365,:,:,:)=[];

% Find nearest indexes for locations:
Xloc = [xmin,xmax];
Yloc = [ymin,ymax];

LATITUDE_T = ncread('/ds/projects/iomp/obs/ECCO2/cube92_real/THETA.nc/THETA.1440x720x50.19920102.nc','LATITUDE_T',[Yloc(1)],[Yloc(2)-Yloc(1)+1]);
LONGITUDE_T = ncread('/ds/projects/iomp/obs/ECCO2/cube92_real/THETA.nc/THETA.1440x720x50.19920102.nc','LONGITUDE_T',[Xloc(1)],[Xloc(2)-Xloc(1)+1]);
DEPTH_T = ncread('/ds/projects/iomp/obs/ECCO2/cube92_real/THETA.nc/THETA.1440x720x50.19920102.nc','DEPTH_T');

[lons lats] = meshgrid(LONGITUDE_T,LATITUDE_T);

depth = DEPTH_T;

%% Calculate vertical levels:
%ncload(grdname,'h','zice','lat_rho','lon_rho','mask_rho','mask_zice')
h=ncread(grdname,'h')';
zice=ncread(grdname,'zice')';
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
lat_u=ncread(grdname,'lat_u')';
lon_u=ncread(grdname,'lon_u')';
lat_v=ncread(grdname,'lat_v')';
lon_v=ncread(grdname,'lon_v')';
mask_rho=ncread(grdname,'mask_rho')';
mask_zice=ncread(grdname,'mask_zice')';

addpath(genpath('/ds/projects/iomp/matlab_scripts/ROMS_MATLAB/'))

h = h.*mask_rho;
zice = zice.*mask_zice;
hc = 20;
x = lon_rho;
y = lat_rho;
kgrid = 0;
column = 0;
plt = 0;
index = 1;
[z,s,Cs_r]=scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, kgrid, column, index, plt);
kgrid = 1;
[z,s,Cs_w]=scoord(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, 1, column, index, plt);

%% For rho points:

InterpSurfaceIni    = nan(length(Cs_r),size(h,1)*2+size(h,2)-2);
TotalWCT            = [h(2:end,1)',h(end,:),fliplr(h(1:end-1,end)')];
MaskVerSec          = repmat([mask_rho(2:end,1)',mask_rho(end,:),fliplr(mask_rho(1:end-1,end)')],length(Cs_r),1);
ind1 = find(MaskVerSec == 0); MaskVerSecNaN = MaskVerSec;
MaskVerSecNaN(ind1)= NaN;
DepthsInterpSurface = flipud(repmat(TotalWCT,length(Cs_r),1).*repmat(Cs_r',1,size(TotalWCT,2)));

LatInterpSurface    = repmat([lat_rho(2:end,1)',lat_rho(end,:),fliplr(lat_rho(1:end-1,end)')],length(Cs_r),1);
LonInterpSurface    = repmat([lon_rho(2:end,1)',lon_rho(end,:),fliplr(lon_rho(1:end-1,end)')],length(Cs_r),1);

%display(['Regridding ' Months(ind6,:)  ' ...'])
isfc.tmp = nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.slt = nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.u =   nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.v =   nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.dpt = nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));
sfc.prs =  nan(12,size(LonInterpSurface,1),size(LonInterpSurface,2));

[Xm,Zm,Ym] = meshgrid(LONGITUDE_T,-DEPTH_T,LATITUDE_T);

InterpFcn = 'linear';
addpath('/ds/projects/iomp/matlab_scripts')
for ii = 1:12; 
ii,
In1 = squeeze(theta(ii,:,:,:)); 
In2 = squeeze(salt(ii,:,:,:));
In3 = squeeze(uvel(ii,:,:,:));
In4 = squeeze(vvel(ii,:,:,:));
    
isfc.tmp(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In1,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.slt(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In2,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.u(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In3,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.v(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In4,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.dpt(ii,:,:) = -DepthsInterpSurface(:,:);
isfc.prs(ii,:,:) = -DepthsInterpSurface(:,:);
save isfcClima.dat ii -ascii
if ~rem(ii,round(12/10))
disp([num2str(ii/12*100),'%'])
end
end

save isfc_cube92.mat isfc -v7.3

DayPerYear=365;

do_ISOM_lbc_nc_cube92(bryname,grdname,'Lateral Boundaries Salt and Temp',[0 1 1 1],Vtransform, Vstretching, Tcline, theta_s,theta_b,20,N,[(DayPerYear/12/2):(DayPerYear/12):(DayPerYear)-(DayPerYear/12/2)],[DayPerYear*(MaxYear-MinYear+1)],'clobber')


