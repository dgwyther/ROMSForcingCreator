% make_tides.m
% Ben Galton-Fenzi | 20?
% David Gwyther updated heavily in Oct 2016

grdname = '/ds/projects/iomp/totten/ana/dgwyther/tisom009/grid/tisom008_canal_grd.nc';
frcname = 'tisom009_tides.nc';


%addpath(genpath('/ds/projects/iomp/matlab_scripts/tmd'))
addpath(genpath('/ds/projects/iomp/matlab_scripts/TMD2.03'))
%addpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib')) 

h = ncread(grdname,'h')';
lat_rho=ncread(grdname,'lat_rho')';
lon_rho=ncread(grdname,'lon_rho')';
mask_rho=ncread(grdname,'mask_rho')';
mask_zice=ncread(grdname,'mask_zice')';

[m n] = size(lat_rho);
Hph = NaN(10,m,n);
Gph = NaN(10,m,n);
h_tpxo = NaN(m,n);

%constituents: m2 s2 n2 k2 k1 o1 p1 q1 mf mm

Model = '/ds/projects/iomp/obs/tides/tpxo/t6.2/Model_tpxo6.2';
%[Hh,Gh,h_tpxo,conList]=extract_HC(Model,lat_rho(:),lon_rho(:),'z');
[Hh,Gh,h_tpxo,conList]=tmd_extract_HC(Model,lat_rho(:),lon_rho(:),'z');
Hh = reshape(Hh,10,m,n);
Gh = reshape(Gh,10,m,n);
% Hu = reshape(Hu,10,m,n);
% Gu = reshape(Gu,10,m,n);
% Hv = reshape(Hv,10,m,n);
% Gv = reshape(Gv,10,m,n);
h_tpxo = reshape(h_tpxo,m,n);



frc_title = 'Tidal Forcing';
vars = {'z'}; %,'u','v'};
constituents  = {'m2','s2','n2','k2','k1','o1','p1','q1','mm','mf'};

%[x,y,h,ma,labels] = roms_grid(grid);
[pp qq] = size(h);
n_constituents     = length(constituents);
tidal_constituents = [];
for i=1:n_constituents
  tidal_constituents = [tidal_constituents,constituents{i},','];
end
tidal_constituents = tidal_constituents(1:end-1);

% --------------------------------------------------------------------
% gen NetCDF file:
% --------------------------------------------------------------------
id = netcdf.create(frcname, 'clobber');

% Global attributes:
disp('Create global attribute')
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', frc_title);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'components',tidal_constituents);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');


% Dimensions:
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;
xi_rho_dim = netcdf.defDim(id, 'xi_rho', Lp);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', Mp);
n_constituents_dim = netcdf.defDim(id, 'n_constituents', n_constituents);


% Variables and attributes:
tide_period_id = netcdf.defVar(id, 'tide_period', 'double', n_constituents_dim);
netcdf.putAtt(id, tide_period_id, 'long_name', 'tide angular period');
netcdf.putAtt(id, tide_period_id, 'units', 'hours');
netcdf.putAtt(id, tide_period_id, 'field', 'tide_period, scalar');

tide_Ephase_id = netcdf.defVar(id, 'tide_Ephase', 'double', [xi_rho_dim eta_rho_dim n_constituents_dim]);
netcdf.putAtt(id, tide_Ephase_id, 'long_name', 'tidal elevation phase angle');
netcdf.putAtt(id, tide_Ephase_id, 'units', 'degrees, time of maximum elevation with respect to chosen time origin');
netcdf.putAtt(id, tide_Ephase_id, 'field', 'tide_Ephase, scalar');

tide_Eamp_id = netcdf.defVar(id, 'tide_Eamp', 'double', [xi_rho_dim eta_rho_dim n_constituents_dim]);
netcdf.putAtt(id, tide_Eamp_id, 'long_name', 'tidal elevation amplitude');
netcdf.putAtt(id, tide_Eamp_id, 'units', 'metres');
netcdf.putAtt(id, tide_Eamp_id, 'field', 'tide_Eamp, scalar');
netcdf.endDef(id);


% --------------------------------------------------------------------
% fill forcing file:
% --------------------------------------------------------------------
if 1 % only force boundaries
 tmask = zeros(pp,qq);
 tmask(:,1) = 1; %tmask(:,2) = 1;
 tmask(:,qq) = 1; %tmask(:,qq-1) = 1;
 tmask(1,:) = 1; %tmask(2,:) = 1;
 tmask(pp,:) = 1; %tmask(pp-1,:) = 1;
elseif 0 % force whole surface
mask_zice = ncread(grdname,'mask_zice')';
tmask = ones(pp,qq);
%ma = mask_rho;
mask_rho(mask_zice == 1) = 0;
end

%nc = netcdf(frc_name,'write');
zamp = nan(size(lon_rho,1),size(lon_rho,2),n_constituents);
zpha = nan(size(lon_rho,1),size(lon_rho,2),n_constituents);
Tcons= nan(1,n_constituents);
addpath('/ds/projects/iomp/matlab_scripts/tmd')
for n=1:n_constituents;
  zza = squeeze(Hh(n,:,:));
  zzp = squeeze(Gh(n,:,:));
%   uua = squeeze(Hu(n,:,:)); uua = uua.*h_tpxo./(h.*100); % cm -> m
%   uup = squeeze(Gu(n,:,:)); 
%   vva = squeeze(Hv(n,:,:)); vva = vva.*h_tpxo./(h.*100);
%   vvp = squeeze(Gv(n,:,:)); 

  iza = find(isnan(zza(:)) == 0);
  izx = x(iza); izy =y(iza); izz = zza(iza);
  zamp(:,:,n) = griddata(izx,izy,izz,lon_rho,lat_rho,'nearest').*mask_rho.*tmask;

  izp = find(isnan(zzp(:)) == 0);
  izx = x(izp); izy = y(izp); izz = zzp(izp);
  zpha(:,:,n) = griddata(izx,izy,izz,lon_rho,lat_rho,'nearest').*mask_rho.*tmask;

  Tcons(n) = 1/name2freq(constituents{n});

%     iua = find(isnan(uua(:)) == 0);
%   iux = x(iua); iuy =y(iua); iuz = zza(iua);
%   uamp = griddata(iux,iuy,iuz,x,y,'nearest').*ma.*tmask;
%   
%   iup = find(isnan(uup(:)) == 0);
%   iux = x(iup); iuy =y(iup); iuz = zza(iup);
%   upha = griddata(iux,iuy,iuz,x,y,'nearest').*ma.*tmask;
%   
%      iva = find(isnan(vva(:)) == 0);
%   ivx = x(iva); ivy =y(iva); ivz = zza(iva);
%   vamp = griddata(ivx,ivy,ivz,x,y,'nearest').*ma.*tmask;
%        
%   ivp = find(isnan(vvp(:)) == 0);
%   ivx = x(ivp); ivy =y(ivp); ivz = zza(ivp);
%   vpha = griddata(ivx,ivy,ivz,x,y,'nearest').*ma.*tmask;
%   
%   zamp = ['out.',vars{1},'_',constituents{n},'_amp']; eval(['zamp = ',zamp,';']);
%   zpha = ['out.',vars{1},'_',constituents{n},'_pha']; eval(['zpha = ',zpha,';']);
%   uamp = ['out.',vars{2},'_',constituents{n},'_amp']; eval(['uamp = ',uamp,';']); uamp = uamp/100; % cm -> m
%   upha = ['out.',vars{2},'_',constituents{n},'_pha']; eval(['upha = ',upha,';']);
%   vamp = ['out.',vars{3},'_',constituents{n},'_amp']; eval(['vamp = ',vamp,';']); vamp = vamp/100;
%   vpha = ['out.',vars{3},'_',constituents{n},'_pha']; eval(['vpha = ',vpha,';']);
%  [sema,ecc,inc,pha]=ap2ep(uamp,upha,vamp,vpha);

%  nc{'tide_Ephase'}(n,:) = nan2zero(zpha);
%  nc{'tide_Eamp'}(n,:)   = nan2zero(zamp);
%   nc{'tide_Cphase'}(n,:) = nan2zero(pha);
%   nc{'tide_Cangle'}(n,:) = nan2zero(inc);
%   nc{'tide_Cmin'}(n,:)   = nan2zero(ecc.*sema);
%   nc{'tide_Cmax'}(n,:)   = nan2zero(sema);
%  nc{'tide_period'}(n)   = 1/name2freq(constituents{n});
end
rmpath('/ds/projects/iomp/matlab_scripts/tmd')

zpha(isnan(zpha))=0;
zamp(isnan(zamp))=0;
netcdf.putVar(id, tide_Ephase_id, zpha);
netcdf.putVar(id, tide_Eamp_id, zamp);
netcdf.putVar(id, tide_period_id, Tcons);

close(nc);

