%% Load ECCO data
load(['cube84_iaf_salt_',RunName,'.mat'])
load(['cube84_iaf_theta_',RunName,'.mat'])
load(['cube84_iaf_uvel_',RunName,'.mat'])
load(['cube84_iaf_vvel_',RunName,'.mat'])

theta(theta < -10) = NaN;
salt(salt < 0) = NaN;

% Load model grid:
%ncload(grdname,'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v')

% Find nearest indexes for locations:
%ncload /u/crcdata/ECCO2/cube84/THETA/THETA.1440x720x50.001.nc LATITUDE_T LONGITUDE_T DEPTH_T
LATITUDE_T = ncread('/ds/projects/iomp/obs/ECCO2/cube78/THETA/THETA.1440x720x50.month001.nc','LATITUDE_T')';
LONGITUDE_T = ncread('/ds/projects/iomp/obs/ECCO2/cube78/THETA/THETA.1440x720x50.month001.nc','LONGITUDE_T')';
DEPTH_T = ncread('/ds/projects/iomp/obs/ECCO2/cube78/THETA/THETA.1440x720x50.month001.nc','DEPTH_T')';

[lons lats] = meshgrid(LONGITUDE_T,LATITUDE_T);

depth = DEPTH_T;
% % For TISOM:
% xmax = 525;
% xmin = 410; % Indices for this are made in make_lbc.m
% ymax = 125;
% ymin = 80;

Xloc = [xmin:xmax];
Yloc = [ymin:ymax];

Var1 = theta;
Var2 = salt;

X = lons(Yloc,Xloc);
Y = lats(Yloc,Xloc);

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
% Vtransform = 1;
% Vstretching = 2;
% theta_s = 0.9;
% theta_b = 4;
%Vtransform = 2;
%Vstretching = 4;
%theta_s = 2;
%theta_b = 2;
%N = 19;
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
Bathy = [];
for i = 1:size(theta,3),
  for j = 1:size(theta,4),

ii = find(isnan(squeeze(theta(1,:,i,j))) == 1);

if length(ii) > 0
Bathy(j,i) = depth(ii(1));
else
Bathy(j,i) = depth(end);
end
end
end

BathyInterp = interp2(X,Y,Bathy,LonInterpSurface,LatInterpSurface,'linear');

L1 = length(mask_rho(2:end,1));
L2 = length(lat_rho(end,:))+L1;
L3 = length(lat_rho(1:end-1,end))+L2;


XX = [];
YY = [];
for i = 1:50;
XX(i,:,:) = X';
YY(i,:,:) = Y';
end

ZZ = [];
for i = 1:size(Bathy,1);
    for j = 1:size(Bathy,2);
ZZ(:,j,i) = -depth;
    end
end

%display(['Regridding ' Months(ind6,:)  ' ...'])
isfc.tmp = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.slt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.u =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.v =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.dpt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
sfc.prs =  nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));

[Xm,Zm,Ym] = meshgrid(LONGITUDE_T(Xloc),-DEPTH_T,LATITUDE_T(Yloc));

InterpFcn = 'linear';
addpath('/ds/projects/iomp/matlab_scripts')
for i = 1:12*(MaxYear-MinYear+1); 
i,
In1 = squeeze(Var1(i,:,:,:)); 
In2 = squeeze(Var2(i,:,:,:));
In3 = squeeze(uvel(i,:,:,:));
In4 = squeeze(vvel(i,:,:,:));
    
isfc.tmp(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In1,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.slt(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In2,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.u(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In3,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.v(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In4,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.dpt(i,:,:) = -DepthsInterpSurface(:,:);
isfc.prs(i,:,:) = -DepthsInterpSurface(:,:);
save isfcMonth.dat i -ascii
end

temp1 = reshape(isfc.tmp,12,[],size(isfc.tmp,2),size(isfc.tmp,3));
temp2 = reshape(isfc.slt,12,[],size(isfc.slt,2),size(isfc.slt,3));
temp3 = reshape(isfc.u,12,[],size(isfc.u,2),size(isfc.u,3));
temp4 = reshape(isfc.v,12,[],size(isfc.v,2),size(isfc.v,3));
temp5 = reshape(isfc.dpt,12,[],size(isfc.dpt,2),size(isfc.dpt,3));
temp6 = reshape(isfc.prs,12,[],size(isfc.prs,2),size(isfc.prs,3));

isfc.tmp = squeeze(nanmean(temp1,2));
isfc.slt = squeeze(nanmean(temp2,2));
isfc.u = squeeze(nanmean(temp3,2));
isfc.v = squeeze(nanmean(temp4,2));
isfc.dpt = squeeze(nanmean(temp5,2));
isfc.prs = squeeze(nanmean(temp6,2));


save isfc_cube84.mat isfc -v7.3

% % % %%TEST
% % % 
% % % lon = LONGITUDE_T(Xloc);
% % % lat = LATITUDE_T(Yloc);
% % % depth = -DEPTH_T;
% % % 
% % % [Xm,Ym,Zm] = meshgrid(lon,depth,lat);
% % % 
% % % ans = interp3(Xm,Ym,Zm,In1,LonInterpSurface,LatInterpSurface,DepthsInterpSurface,InterpFcn);
% % % 
% % % 
% % % figure, flat_pcolor(squeeze(isfc.tmp(1,:,:)));
% % % colorbar,


do_ISOM_lbc_nc_cube84(bryname,grdname,'Lateral Boundaries Salt and Temp',[0 1 1 1],Vtransform, Vstretching, Tcline, theta_s,theta_b,20,N,[(30.4375/2):30.4375:(30.4375*12)-(30.4375/2)],[30.4375*12],'clobber')

