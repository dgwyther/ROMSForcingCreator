% Load OKane-CORE data
load(['OKane-COREv1_iaf_salt_',RunName,'.mat'])
load(['OKane-COREv1_iaf_theta_',RunName,'.mat'])
load(['OKane-COREv1_iaf_uvel_',RunName,'.mat'])
load(['OKane-COREv1_iaf_vvel_',RunName,'.mat'])

theta(theta < -10) = NaN;
salt(salt < 0) = NaN;

LATITUDE=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_temp_100E-130E_60S-68S.nc'],'yt_ocean');
LONGITUDE=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_temp_100E-130E_60S-68S.nc'],'xt_ocean');
DEPTHS=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_temp_100E-130E_60S-68S.nc'],'st_ocean');
LATITUDEu=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_u_100E-130E_60S-68S.nc'],'yu_ocean');
LONGITUDEu=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_u_100E-130E_60S-68S.nc'],'xu_ocean');
LATITUDEv=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_v_100E-130E_60S-68S.nc'],'yu_ocean');
LONGITUDEv=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_v_100E-130E_60S-68S.nc'],'xu_ocean');


time=ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_temp_100E-130E_60S-68S.nc'],'time');

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
isfc.tmp = nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.slt = nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.u =   nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.v =   nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));
isfc.dpt = nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));
sfc.prs =  nan(length(time),size(LonInterpSurface,1),size(LonInterpSurface,2));

[Xm,Zm,Ym] = meshgrid(LONGITUDE,-DEPTHS,LATITUDE);
[Xmu,Zmu,Ymu] = meshgrid(LONGITUDEu,-DEPTHS,LATITUDEu);
[Xmv,Zmv,Ymv] = meshgrid(LONGITUDEv,-DEPTHS,LATITUDEv);

InterpFcn = 'linear';
addpath('/ds/projects/iomp/matlab_scripts')
for ii = 1:length(time);
In1 = squeeze(theta(ii,:,:,:)); 
In2 = squeeze(salt(ii,:,:,:));
In3 = squeeze(uvel(ii,:,:,:));
In4 = squeeze(vvel(ii,:,:,:));
    
isfc.tmp(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In1,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.slt(ii,:,:) = inpaint_nans(interp3(Xm,Zm,Ym,In2,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.u(ii,:,:) = inpaint_nans(interp3(Xmu,Zmu,Ymu,In3,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.v(ii,:,:) = inpaint_nans(interp3(Xmv,Zmv,Ymv,In4,LonInterpSurface,DepthsInterpSurface,LatInterpSurface,InterpFcn),2);
isfc.dpt(ii,:,:) = -DepthsInterpSurface(:,:);
isfc.prs(ii,:,:) = -DepthsInterpSurface(:,:);
save isfcClima.dat ii -ascii
if ~rem(ii,round(length(time)/10))
disp([num2str(ii/length(time)*100),'%'])
end
end

save isfc_OKane-COREv1_iaf.mat isfc -v7.3



do_ISOM_lbc_nc_OKane_corev1(bryname,grdname,'Lateral Boundaries Salt and Temp',[0 1 1 1],Vtransform, Vstretching, Tcline, theta_s,theta_b,20,N,[(30.4375/2):30.4375:(30.4375*12)*100-(30.4375/2)],[30.4375*12*100],'clobber')
