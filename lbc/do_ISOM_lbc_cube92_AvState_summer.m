set(0,'DefaultFigureRenderer','zbuffer')
% ObsDir = '/v/hpchome/bkgalton/massdata/obs/';
% LevitusFileNameTmp = 'ECCO2/otemp.anal1deg.nc';
% LevitusFileNameSal = 'levitus/salt.anal1deg.nc';
% cd /v/hpchome/bkgalton/massdata/obs/tracers
%                                                               %  ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
%
% GrdLocation  = '/v/hpchome/bkgalton/massdata/ais';            % ALL THIS PROBABLY NEEDS TO BE CHANGED!
% 
% RunID = '/ais062_1';                                          %  \/ \/ \/ \/ \/ 

% grdname='/ds/projects/iomp/totten/ana/dgwyther/netcdf/grid/tisom002_grd.n
% c'; % Loaded by make_lbc.m

%% Load ECCO data

load([WorkingDir,'cube92_iaf_salt_',ModelName,RunNo,'.mat'])
load([WorkingDir,'cube92_iaf_theta_',ModelName,RunNo,'.mat'])
load([WorkingDir,'cube92_iaf_uvel_',ModelName,RunNo,'.mat'])
load([WorkingDir,'cube92_iaf_vvel_',ModelName,RunNo,'.mat'])

theta(theta < -10) = NaN;
salt(salt < 0) = NaN;

% Load model grid:
ncload(grdname,'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v')

% Find nearest indexes for locations:
ncload /ds/projects/iomp/obs/ECCO2/cube92_real/THETA_monthly.nc/THETA.1440x720x50.199201.nc LATITUDE_T LONGITUDE_T DEPTH_T
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
ncload(grdname,'h','zice','lat_rho','lon_rho','mask_rho','mask_zice')
h = h.*mask_rho;
zice = zice.*mask_zice;
hc = 20;
x = lon_rho;
 y = lat_rho;
% Vtransform = 1;
% Vstretching = 2;
% theta_s = 0.9;
% theta_b = 4;
Vtransform = 2;
Vstretching = 3;
theta_s = 2;
theta_b = 2;
N = 31;
kgrid = 0;
column = 0;
plt = 0;
index = 1;
[z,s,Cs_r]=scoord_zice(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, kgrid, column, index, plt);
kgrid = 1;
[z,s,Cs_w]=scoord_zice(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
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
tmp = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
slt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
u =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
v =   nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
dpt = nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));
prs =  nan(12*(MaxYear-MinYear+1),size(LonInterpSurface,1),size(LonInterpSurface,2));

[Xm,Zm,Ym] = meshgrid(LONGITUDE_T(Xloc),-DEPTH_T,LATITUDE_T(Yloc));

InterpFcn = 'linear';
for i = 1:12*(MaxYear-MinYear+1); 
i,
In1 = squeeze(Var1(i,:,:,:)); 
    In2 = squeeze(Var2(i,:,:,:));
    In3 = squeeze(uvel(i,:,:,:));
    In4 = squeeze(vvel(i,:,:,:));
    
tmp(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In1,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
slt(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In2,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
u(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In3,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
v(i,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In4,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
dpt(i,:,:) = -DepthsInterpSurface(:,:);
prs(i,:,:) = -DepthsInterpSurface(:,:);
save isfcMonth.dat i -ascii
end

temp1 = reshape(tmp,12,[],size(tmp,2),size(tmp,3));
temp2 = reshape(slt,12,[],size(slt,2),size(slt,3));
temp3 = reshape(u,12,[],size(u,2),size(u,3));
temp4 = reshape(v,12,[],size(v,2),size(v,3));
temp5 = reshape(dpt,12,[],size(dpt,2),size(dpt,3));
temp6 = reshape(prs,12,[],size(prs,2),size(prs,3));

clim.tmp = squeeze(nanmean(temp1,2));
clim.slt = squeeze(nanmean(temp2,2));
clim.u = squeeze(nanmean(temp3,2));
clim.v = squeeze(nanmean(temp4,2));
clim.dpt = squeeze(nanmean(temp5,2));
clim.prs = squeeze(nanmean(temp6,2));

tmp2 = squeeze(nanmean(clim.tmp([1,2,12],:,:),1));
slt2 = squeeze(nanmean(clim.slt([1,2,12],:,:),1));
u2 = squeeze(nanmean(clim.u([1,2,12],:,:),1));
v2 = squeeze(nanmean(clim.v([1,2,12],:,:),1));
dpt2 = squeeze(nanmean(clim.dpt([1,2,12],:,:),1));
prs2 = squeeze(nanmean(clim.prs([1,2,12],:,:),1));

for jj = 1:12%size(isfc.tmp,1);
isfc.tmp(jj,:,:) = tmp2;
isfc.slt(jj,:,:) = slt2;
isfc.u(jj,:,:) = u2;
isfc.v(jj,:,:) = v2;
isfc.dpt(jj,:,:) = dpt2;
isfc.prs(jj,:,:) = prs2;
end



save isfc_AvState.mat isfc -v7.3

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


 do_ISOM_lbc_nc_cube92_AvState(bryname,grdname,'Lateral Boundaries Salt and Temp',[0 1 1 1],...
   5,0.4,20,31,...
 [(30.4375/2):30.4375:(30.4375*12)-(30.4375/2)],[30.4375*12],'clobber')

