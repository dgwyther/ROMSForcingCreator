grdname = 'aisom004_grd.nc';
ModelName = 'Amery Ice Shelf Ocean Model v10';

xc = [56,85]; 
yc = -[-73.4,-61];
resolution = 1/10; 

plots = 0;

%xc = [110,130]; %xc totten
%yc = -[-67.65,-60];%yc totten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Built on several decades of scripting, from:
%  Pierrick Penven (Original Seagrid scripts from ROMSTOOLS/ROMS_AGRIF)
%  Ben Galton-Fenzi (Major conversion for our application of ROMS)
%  Eva Cougnon (Many bug fixes)
%  David Gwyther (Many bug fixes, conversion to new netcdf library)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------------------------------------
title =([ModelName,' grid file']); % Make title from model name variable
addpath('/ds/projects/iomp/matlab_scripts')
addpath('/ds/projects/iomp/matlab_scripts/ROMS_NetCDF')
addpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib'))

read_asaid_islands_HL

lon = ncread('/ds/projects/iomp/obs/geom/RTopo105/RTopo105_50S.nc','lon')';
lat = ncread('/ds/projects/iomp/obs/geom/RTopo105/RTopo105_50S.nc','lat')';
bathy = double(ncread('/ds/projects/iomp/obs/geom/RTopo105/RTopo105_50S.nc','bathy')');
amask = double(ncread('/ds/projects/iomp/obs/geom/RTopo105/RTopo105_50S.nc','amask')');
draft = double(ncread('/ds/projects/iomp/obs/geom/RTopo105/RTopo105_50S.nc','draft')');

[~,xc1]=min(abs(lon-xc(1))); [~,xc2]=min(abs(lon-xc(2)));
[~,yc1]=min(abs(lat- -yc(1))); [~,yc2]=min(abs(lat- -yc(2)));

Xrng = [floor(0.95*xc1):floor(1.05*xc2)];
Yrng = [floor(0.95*yc1):floor(1.05*yc2)]; %done programmatically now, with 5% leeway

load /ds/projects/iomp/obs/geom/RTopo105b/RTopo105_coast.asc
load /ds/projects/iomp/obs/geom/RTopo105b/RTopo105_gl.asc

[lons lats] = meshgrid(lon,lat);
e = nan(size(lons));
n = nan(size(lats));
for i = 1:size(e,1);
   [t1 t2] = polar_stereo_deluxe(squeeze(lats(i,:)),squeeze(lons(i,:)),0,0,1,-71);

e(i,:) = t1';
n(i,:)= t2';
end

nhind = find(RTopo105_coast(:,2) >= -60);
RTopo105_coast(nhind,2) = NaN;

[ecl ncl] = polar_stereo_deluxe(RTopo105_coast(:,2),RTopo105_coast(:,1),0,0,1,-71);
[egl ngl] = polar_stereo_deluxe(RTopo105_gl(:,2),RTopo105_gl(:,1),0,0,1,-71);

% 

lightgray = 0.8.*ones(3,1);
darkgray = 0.5.*ones(3,1);

% 1) Select join latitude:
JoinLat = pi()/180 * yc(2);%deg2rad(yc(2));
pion2 = pi()/2;
pi2 = 2*pi();
Lats = [];
MinLat = yc(2);
ii = 1;
while MinLat < yc(1)

    Lats(ii) = MinLat;
    MinLat = MinLat+(cosd(MinLat)*resolution);
    ii = ii+1;
end
Lats = fliplr(-Lats);
Lons = [xc(1):resolution:xc(2)];
%Lons(Lons > xc(2)) = []; 
[pg.lon pg.lat] = meshgrid(Lons,Lats);
[pg.x,pg.y]=polar_stereo_deluxe(pg.lat,pg.lon,0,0,1,-71);
pg.x=pg.x*1000; pg.y=pg.y*1000; %km->m

if plots
res = 10;
figure, 
flat_pcolor(e(1:res:end,1:res:end),n(1:res:end,1:res:end),bathy(1:res:end,1:res:end));
hold on, 
res = 20;
plot(pg.x(1:res:end,1:res:end)/1000,pg.y(1:res:end,1:res:end)/1000,'-','color',darkgray);    
plot(pg.x(1:res:end,1:res:end)'/1000,pg.y(1:res:end,1:res:end)'/1000,'-','color',darkgray); 
plot(pg.x([1,end],[1,end])/1000,pg.y([1,end],[1,end])/1000,'k-');  
plot(pg.x([1,end],[1:end])'/1000,pg.y([1,end],[1:end])'/1000,'k-');  
axis([1200 3500 0 2000])
hold on, plot(ecl,ncl,'k.');   
plot(egl,ngl,'k.');

lon_rho=ncread('/ds/projects/iomp/aisom/mdl/aisom001/ais_grdpolarb_truncated.nc','lon_rho')';
lat_rho=ncread('/ds/projects/iomp/aisom/mdl/aisom001/ais_grdpolarb_truncated.nc','lat_rho')';

[ox,oy]=polar_stereo_deluxe(lat_rho,lon_rho,0,0,1,-71);
%ox=ox*1000; oy=oy*1000; %km->m
res = 10;
plot(ox(1:res:end,1:res:end),oy(1:res:end,1:res:end),'b-');
plot(ox(1:res:end,1:res:end)',oy(1:res:end,1:res:end)','b-');
plot(ox([1,end],[1,end]),oy([1,end],[1,end]),'b-');
plot(ox([1,end],[1:end])',oy([1,end],[1:end])','b-');


%Xrng = [1.4e4:1.6e4]; %in km
%Yrng =  [900:1800];   %in km


figure, flat_pcolor(lons(Yrng,Xrng),lats(Yrng,Xrng),double(amask(Yrng,Xrng)))
hold on,
plot(pg.lon([1,end],[1,end]),pg.lat([1,end],[1,end]),'k-');
plot(pg.lon([1,end],[1:end])',pg.lat([1,end],[1:end])','k-');
end

wm = amask(Yrng,Xrng);
wm(wm ~= 0 & wm ~= 2) = NaN;
wm(~isnan(wm)) = 1;
wm(isnan(wm)) = 0;

im = amask(Yrng,Xrng);
im(im ~= 2) = 0; im(im == 2) = 1;

% Load rtopo grounding line information:

%     pg.delr      = fliplr(pg.delr');
%     pg.deltheta  = fliplr(pg.deltheta');
%     pg.cosrx     = fliplr(pg.cosrx');
%     pg.cosry     = fliplr(pg.cosry');
%     pg.costhetax = fliplr(pg.costhetax');
%     pg.costhetay = fliplr(pg.costhetay');
    pg.h1      = double(-interp2(lons(Yrng,Xrng),lats(Yrng,Xrng),bathy(Yrng,Xrng),pg.lon,pg.lat)); % This is the bathymetry data
    pg.zice    = double(interp2(lons(Yrng,Xrng),lats(Yrng,Xrng),draft(Yrng,Xrng),pg.lon,pg.lat)); % This is the bathymetry data

    pg.icemask    = double(interp2(lons(Yrng,Xrng),lats(Yrng,Xrng),im,pg.lon,pg.lat,'nearest')); 
%NEED TO DO
    pg.watermask = double(interp2(lons(Yrng,Xrng),lats(Yrng,Xrng),wm,pg.lon,pg.lat,'nearest')); %water mask
%% Post loading data editing...

theInterpFcn = 'interp2';
theInterpMethod = 'nearest';

xtmp_zappa         = feval(theInterpFcn, pg.x, 1, theInterpMethod);
ytmp_zappa         = feval(theInterpFcn, pg.y, 1, theInterpMethod);


clear xtmp8 ytmp8 xtmp ytmp
for indi=1:size(pg.x,1)
  xtmp8(2*indi-1,:) = pg.x(indi,:);
  xtmp8(2*indi,:) = pg.x(indi,:);
  ytmp8(2*indi-1,:) = pg.y(indi,:);
  ytmp8(2*indi,:) = pg.y(indi,:);
end

for indi=1:size(pg.x,2)
    xtmp(:,2*indi-1) = xtmp8(:,indi);
    xtmp(:,2*indi) = xtmp8(:,indi);
    ytmp(:,2*indi-1) = ytmp8(:,indi);
    ytmp(:,2*indi) = ytmp8(:,indi);
end

xtmp(:,1)=[];
ytmp(:,1)=[];

xtmp(1,:)=[];
ytmp(1,:)=[];


[m, n] = size(xtmp);
if ~rem(m, 2), m = m-1; end   % m, n must be odd.
if ~rem(n, 2), n = n-1; end

i1 = 3:2:m-2; j1 = 3:2:n-2;
i2 = 2:2:m-1; j2 = 2:2:n-1;

xtmp2 = xtmp(i1,j1);
ytmp2 = ytmp(i1,j1);

s.x = xtmp(i2,j2);
s.y = ytmp(i2,j2);
s.zb = griddata(pg.x,pg.y,pg.h1,xtmp2,ytmp2,theInterpMethod);
s.zd = griddata(pg.x,pg.y,pg.zice,xtmp2,ytmp2,theInterpMethod);
s.mw = griddata(pg.x,pg.y,pg.watermask,xtmp2,ytmp2,theInterpMethod);
s.mi = griddata(pg.x,pg.y,pg.icemask,xtmp2,ytmp2,theInterpMethod);

s.mw = logical(s.mw);
s.mi = logical(s.mi);

%  s.zb = pg.h1(i2,j2);
%  s.zd = pg.zice;
%  s.mw = pg.watermask;
%  s.mi = pg.icemask;
% 
 s.zb(1,:) = s.zb(3,:);
s.zd(1,:) = s.zd(3,:);
 s.zb(2,:) = s.zb(4,:);
 s.zd(2,:) = s.zd(4,:);
%tangle = 0.5*(pg.cosrx(2:end,:)+pg.cosrx(1:end-1,:));
s.angle = zeros(size(pg.x(2:end,2:end))); %-0.5*(tangle(:,2:end)+tangle(:,1:end-1));
[m n] = size(s.x);

%Xo = 500000;
%Yo = 10000000;
%C = 146;
%S = 0.9996;
% [s.lon s.lat]=rf(s.x(:),s.y(:),C,S,Xo,Yo,1);
 % 
 s.lat = pg.lat; %reshape(s.lat,m,n);
 s.lon = pg.lon;%reshape(s.lon,m,n);
% 
% clear Xo Yo C S

s.clipping_depths(1) = 20;
s.clipping_depths(2) = 5000;

clear m n
%% Create the Roms NetCDF file.

grid_x = s.x;
grid_y = s.y;


[m, n] = size(grid_x);

geogrid_lon = s.lon;
geogrid_lat = s.lat;

%% Need to calculate dely and delx

%%% NOTE: NOT SURE ABOUT THIS!!
tempdelx=[];
tempdely=[];
s.delx = [];
s.dely = [];
for ind = 1:size(s.x,2)-1
    s.delx(:,ind) = ( (s.x(:,ind+1) - s.x(:,ind)).^2 + (s.y(:,ind+1) - s.y(:,ind)).^2 ).^0.5;
end

ind=[];
for ind = 1:size(s.x,1)-1
    s.dely(ind,:) = ( (s.x(ind+1,:) - s.x(ind,:)).^2 + (s.y(ind+1,:) - s.y(ind,:)).^2 ).^0.5;
end

s.delx(end,:)=[]; %Chop off last row, since the equivalent row is removed in dely
s.dely(:,end)=[]; %Chop off last column, since equivalent column is not used in delx


geometry{1} = s.delx;
geometry{2} = s.dely;

mask = ~~s.mw;
land = mask;
water = ~land;

s.mw00 = s.mw;
s.mw00(isnan(s.mw)) = 0;

projection = 'Mercator';
% ROMS needs radians.
min_depth = s.clipping_depths(1);
max_depth = s.clipping_depths(2);
%min_depth_ice = 50;

% Smooth bathymetry and icedraft
raw_s = s;
resss = 0.5;

%ii = find(s.zd >= -min_depth_ice);
%s.zd(ii) = -min_depth_ice;

%% Bit more cleaning:
%s = raw_s;
bathymetry = s.zb;
s.mi(s.mi+s.mw == 1) = 0;
ice_draft = s.zd.*s.mi;

bathymetry = smoothgrid(bathymetry,ones(size(s.mw)),min_depth,max_depth,resss,1,2).*s.mw;

raw_s.blahblah = bathymetry;

% Clip Bathymetry

bathymetry(find(isnan(bathymetry))) = min_depth;
bathymetry(bathymetry<min_depth) = min_depth;
bathymetry((bathymetry + ice_draft) > max_depth) = max_depth;

%ice_draft = -smoothgrid(-ice_draft,s.mi,50,1000,resss,1,2);
wct = (bathymetry + ice_draft);
bathymetry(wct.*s.mi < 0) = -ice_draft(wct.*s.mi < 0) + min_depth;
wct = (bathymetry + ice_draft).*s.mw00;
ii = find(wct < min_depth);
bathymetry(ii) = -ice_draft(ii) + min_depth;
%pp = isnan(bathymetry);
%qq = isnan(ice_draft);

%ice_draft(qq) = 0;
%bathymetry(pp) = min_depth;
wct = bathymetry + ice_draft;

ice_draft(ice_draft > 0) = 0;  % < ---------------------------- What a strange line!
%bathymetry = wct;
%spaced_x = s.spaced_grids{1};
%spaced_y = s.spaced_grids{2};

% Double the grid-size before proceeding.
%  The grid-cell-centers are termed the "rho" points.

theInterpFcn = 'interp2';
theInterpMethod = 'spline';

grid_x = feval(theInterpFcn, grid_x, 1, theInterpMethod);  
grid_y = feval(theInterpFcn, grid_y, 1, theInterpMethod);
%geogrid_lon = feval(theInterpFcn, geogrid_lon, 1, theInterpMethod); <<<---- BUG FIX
%geogrid_lat = feval(theInterpFcn, geogrid_lat, 1, theInterpMethod); <<<---- DEG 3/2016
[geogrid_lat,geogrid_lon] = inverse_polar_stereo(grid_x/1000,grid_y/1000,0,0,1,-71);
%spaced_x = feval(theInterpFcn, spaced_x, 1, theInterpMethod);
%spaced_y = feval(theInterpFcn, spaced_y, 1, theInterpMethod);

[n, m] = size(grid_x);
xl = max(grid_x(:)) - min(grid_x(:));
el = max(grid_y(:)) - min(grid_y(:));


% Save results in NetCDF file
id = netcdf.create(grdname, 'clobber');

%% Global attributes:

%disp(' ## Defining Global Attributes...')
%
%nc.type = ncchar('Gridpak file');
%nc.gridid = ModelName;
%nc.history = ncchar(['Created by "' mfilename '" on ' datestr(now)]);
%
%nc.CPP_options = ncchar('DCOMPLEX, DBLEPREC, NCARG_32, PLOTS,');
%name(nc.CPP_options, 'CPP-options')

% The SeaGrid is now a full array, whose height
%  and width are odd-valued.  We extract staggered
%  sub-grids for the Roms scheme, ignoring the
%  outermost rows and columns.  Thus, the so-called
%  "rho" points correspond to the even-numbered points
%  in an (i, j) Matlab array.  The "psi" points begin
%  at i = 3 and j = 3.  The whole set is indexed as
%  follows:

% rho (2:2:end-1, 2:2:end-1), i.e. (2:2:m, 2:2:n), etc.
% psi (3:2:end-2, 3:2:end-2)
% u   (2:2:end-1, 3:2:end-2)
% v   (3:2:end-2, 2:2:end-1)

if ~rem(m, 2), m = m-1; end   % m, n must be odd.
if ~rem(n, 2), n = n-1; end

i_rho = 2:2:m-1; j_rho = 2:2:n-1;
i_psi = 3:2:m-2; j_psi = 3:2:n-2;
i_u   = 3:2:m-2; j_u   = 2:2:n-1;
i_v   = 2:2:m-1; j_v   = 3:2:n-2;

% The xi direction (left-right):

LP = (m-1)/2;   % The rho dimension.
L = LP-1;       % The psi dimension.

% The eta direction (up-down):

MP = (n-1)/2;   % The rho dimension.
M = MP-1;       % The psi dimension.

disp(' ## Defining Dimensions...')

xi_psi_dim = netcdf.defDim(id, 'xi_psi', L);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', LP);
xi_u_dim = netcdf.defDim(id, 'xi_u', L);
xi_v_dim = netcdf.defDim(id, 'xi_v', LP);

eta_psi_dim = netcdf.defDim(id, 'eta_psi', M);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', MP);
eta_u_dim = netcdf.defDim(id, 'eta_u', MP);
eta_v_dim = netcdf.defDim(id, 'eta_v', M);

one_dim = netcdf.defDim(id, 'one', 1);
two_dim = netcdf.defDim(id, 'two', 2);
bath_dim = netcdf.defDim(id, 'bath', 0);

%nc('xi_psi') = L;
%nc('xi_rho') = LP;
%nc('xi_u') = L;
%nc('xi_v') = LP;

%nc('eta_psi') = M;
%nc('eta_rho') = MP;
%nc('eta_u') = MP;
%nc('eta_v') = M;

%nc('two') = 2;
%nc('bath') = 0; %% (record dimension)
%% Variables and attributes:

disp(' ## Defining Variables and Attributes...')

xl_id = netcdf.defVar(id, 'xl', 'double', one_dim);
netcdf.putAtt(id, xl_id, 'long_name', 'domain length in the XI-direction');
netcdf.putAtt(id, xl_id, 'units', 'meter');
%nc{'xl'} = ncdouble; %% 1 element.
%nc{'xl'}.long_name = ncchar('domain length in the XI-direction');
%nc{''}.units = ncchar('meter');

el_id = netcdf.defVar(id, 'el', 'double', one_dim);
netcdf.putAtt(id, el_id, 'long_name', 'domain length in the ETA-direction');
netcdf.putAtt(id, el_id, 'units', 'meter');
%nc{'el'} = ncdouble; %% 1 element.
%nc{'el'}.long_name = ncchar('domain length in the ETA-direction');
%nc{'el'}.units = ncchar('meter');

JPRJ_id = netcdf.defVar(id, 'JPRJ', 'char', two_dim);
netcdf.putAtt(id, JPRJ_id, 'long_name', 'Map projection type');

%nc{'JPRJ'} = ncchar('two'); %% 2 elements.
%nc{'JPRJ'}.long_name = ncchar('Map projection type');
%nc{'JPRJ'}.option_ME_ = ncchar('Mercator');
%nc{'JPRJ'}.option_ST_ = ncchar('Stereographic');
%nc{'JPRJ'}.option_LC_ = ncchar('Lambert conformal conic');
%name(nc{'JPRJ'}.option_ME_, 'option(ME)')
%name(nc{'JPRJ'}.option_ST_, 'option(ST)')
%name(nc{'JPRJ'}.option_LC_, 'option(LC)')

%nc{'PLAT'} = ncfloat('two'); %% 2 elements.
%nc{'PLAT'}.long_name = ncchar('Reference latitude(s) for map projection');
%nc{'PLAT'}.units = ncchar('degree_north');

%nc{'PLONG'} = ncfloat; %% 1 element.
%nc{'PLONG'}.long_name = ncchar('Reference longitude for map projection');
%nc{'PLONG'}.units = ncchar('degree_east');

%nc{'ROTA'} = ncfloat; %% 1 element.
%nc{'ROTA'}.long_name = ncchar('Rotation angle for map projection');
%nc{'ROTA'}.units = ncchar('degree');

%nc{'JLTS'} = ncchar('two'); %% 2 elements.
%nc{'JLTS'}.long_name = ncchar('How limits of map are chosen');
%nc{'JLTS'}.option_CO_ = ncchar('P1, .. P4 define two opposite corners ');
%nc{'JLTS'}.option_MA_ = ncchar('Maximum (whole world)');
%nc{'JLTS'}.option_AN_ = ncchar('Angles - P1..P4 define angles to edge of domain');
%nc{'JLTS'}.option_LI_ = ncchar('Limits - P1..P4 define limits in u,v space');
%name(nc{'JLTS'}.option_CO_, 'option(CO)')
%name(nc{'JLTS'}.option_MA_, 'option(MA)')
%name(nc{'JLTS'}.option_AN_, 'option(AN)')
%name(nc{'JLTS'}.option_LI_, 'option(LI)')

%nc{'P1'} = ncfloat; %% 1 element.
%nc{'P1'}.long_name = ncchar('Map limit parameter number 1');
%nc{'P2'} = ncfloat; %% 1 element.
%nc{'P2'}.long_name = ncchar('Map limit parameter number 2');

%nc{'P3'} = ncfloat; %% 1 element.
%nc{'P3'}.long_name = ncchar('Map limit parameter number 3');

%nc{'P4'} = ncfloat; %% 1 element.
%nc{'P4'}.long_name = ncchar('Map limit parameter number 4');

%nc{'XOFF'} = ncfloat; %% 1 element.
%nc{'XOFF'}.long_name = ncchar('Offset in x direction');
%nc{'XOFF'}.units = ncchar('meter');

%nc{'YOFF'} = ncfloat; %% 1 element.
%nc{'YOFF'}.long_name = ncchar('Offset in y direction');
%nc{'YOFF'}.units = ncchar('meter');

depthmin_id = netcdf.defVar(id, 'depthmin', 'double', one_dim);
netcdf.putAtt(id, depthmin_id, 'long_name', 'Shallow bathymetry clipping depth');
netcdf.putAtt(id, depthmin_id, 'units', 'meter');
%nc{'depthmin'} = ncshort; %% 1 element.
%nc{'depthmin'}.long_name = ncchar('Shallow bathymetry clipping depth');
%nc{'depthmin'}.units = ncchar('meter');
depthmax_id = netcdf.defVar(id, 'depthmax', 'double', one_dim);
netcdf.putAtt(id, depthmax_id, 'long_name', 'Deepth bathymetry clipping depth');
netcdf.putAtt(id, depthmax_id, 'units', 'meter');
%nc{'depthmax'} = ncshort; %% 1 element.
%nc{'depthmax'}.long_name = ncchar('Deep bathymetry clipping depth');
%nc{'depthmax'}.units = ncchar('meter');

spherical_id = netcdf.defVar(id, 'spherical', 'char', one_dim);
netcdf.putAtt(id, depthmin_id, 'long_name', 'Grid type logical switch');
%nc{'spherical'} = ncchar; %% 1 element.
%nc{'spherical'}.long_name = ncchar('Grid type logical switch');
%nc{'spherical'}.option_T_ = ncchar('spherical');
%nc{'spherical'}.option_F_ = ncchar('Cartesian');
%name(nc{'spherical'}.option_T_, 'option(T)')
%name(nc{'spherical'}.option_F_, 'option(F)')
%nc{'hraw'} = ncdouble('bath', 'eta_rho', 'xi_rho'); %% 0 elements.
%nc{'hraw'}.long_name = ncchar('Working bathymetry at RHO-points');
%nc{'hraw'}.units = ncchar('meter');
%nc{'hraw'}.field = ncchar('bath, scalar');

h_id = netcdf.defVar(id, 'h', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, h_id, 'long_name', 'Final bathymetry at RHO-points');
netcdf.putAtt(id, h_id, 'units', 'meter');
%nc{'h'} = ncdouble('eta_rho', 'xi_rho');
%nc{'h'}.long_name = ncchar('Final bathymetry at RHO-points');
%nc{'h'}.units = ncchar('meter');
%nc{'h'}.field = ncchar('bath, scalar');

zice_id = netcdf.defVar(id, 'zice', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, zice_id, 'long_name', 'Final ice draft at RHO-points');
netcdf.putAtt(id, zice_id, 'units', 'meter');
%nc{'zice'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'zice'}.long_name = ncchar('Final ice draft at RHO-points');
%nc{'zice'}.units = ncchar('meter');
%nc{'zice'}.field = ncchar('draft, scalar');

f_id = netcdf.defVar(id, 'f', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, f_id, 'long_name', 'Coriolis parameter at RHO-points');
netcdf.putAtt(id, f_id, 'units', 'second-1');
%nc{'f'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'f'}.long_name = ncchar('Coriolis parameter at RHO-points');
%nc{'f'}.units = ncchar('second-1');
%nc{'f'}.field = ncchar('Coriolis, scalar');

pm_id = netcdf.defVar(id, 'pm', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, pm_id, 'long_name', 'curvilinear coordinate metric in XI');
netcdf.putAtt(id, pm_id, 'units', 'meter-1');
%nc{'pm'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
%nc{'pm'}.units = ncchar('meter-1');
%nc{'pm'}.field = ncchar('pm, scalar');

pn_id = netcdf.defVar(id, 'pn', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, pn_id, 'long_name', 'curvilinear coordinate metric in ETA');
netcdf.putAtt(id, pn_id, 'units', 'meter-1');
%nc{'pn'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
%nc{'pn'}.units = ncchar('meter-1');
%nc{'pn'}.field = ncchar('pn, scalar');

dndx_id = netcdf.defVar(id, 'dndx', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, dndx_id, 'long_name', 'xi derivative of inverse metric factor pn');
netcdf.putAtt(id, dndx_id, 'units', 'meter-1');
%nc{'dndx'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'dndx'}.long_name = ncchar('xi derivative of inverse metric factor pn');
%nc{'dndx'}.units = ncchar('meter');
%nc{'dndx'}.field = ncchar('dndx, scalar');

dmde_id = netcdf.defVar(id, 'dmde', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, dmde_id, 'long_name', 'eta derivative of inverse metric factor pm');
netcdf.putAtt(id, dmde_id, 'units', 'meter-1');
%nc{'dmde'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'dmde'}.long_name = ncchar('eta derivative of inverse metric factor pm');
%nc{'dmde'}.units = ncchar('meter');
%nc{'dmde'}.field = ncchar('dmde, scalar');

x_rho_id = netcdf.defVar(id, 'x_rho', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, x_rho_id, 'long_name', 'x location of RHO-points');
netcdf.putAtt(id, x_rho_id, 'units', 'meter');
%nc{'x_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'x_rho'}.long_name = ncchar('x location of RHO-points');
%nc{'x_rho'}.units = ncchar('meter');

y_rho_id = netcdf.defVar(id, 'y_rho', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, y_rho_id, 'long_name', 'y location of RHO-points');
netcdf.putAtt(id, y_rho_id, 'units', 'meter');
%nc{'y_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'y_rho'}.long_name = ncchar('y location of RHO-points');
%nc{'y_rho'}.units = ncchar('meter');

x_psi_id = netcdf.defVar(id, 'x_psi', 'double', [xi_psi_dim eta_psi_dim]);
netcdf.putAtt(id, x_psi_id, 'long_name', 'x location of PSI-points');
netcdf.putAtt(id, x_psi_id, 'units', 'meter');
%hc{'x_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
%nc{'x_psi'}.long_name = ncchar('x location of PSI-points');
%nc{'x_psi'}.units = ncchar('meter');

y_psi_id = netcdf.defVar(id, 'y_psi', 'double', [xi_psi_dim eta_psi_dim]);
netcdf.putAtt(id, y_psi_id, 'long_name', 'y location of PSI-points');
netcdf.putAtt(id, y_psi_id, 'units', 'meter');
%nc{'y_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
%nc{'y_psi'}.long_name = ncchar('y location of PSI-points');
%nc{'y_psi'}.units = ncchar('meter');

x_u_id = netcdf.defVar(id, 'x_u', 'double', [xi_u_dim eta_u_dim]);
netcdf.putAtt(id, x_u_id, 'long_name', 'x location of U-points');
netcdf.putAtt(id, x_u_id, 'units', 'meter');
%nc{'x_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
%nc{'x_u'}.long_name = ncchar('x location of U-points');
%nc{'x_u'}.units = ncchar('meter');

y_u_id = netcdf.defVar(id, 'y_u', 'double', [xi_u_dim eta_u_dim]);
netcdf.putAtt(id, y_u_id, 'long_name', 'y location of U-points');
netcdf.putAtt(id, y_u_id, 'units', 'meter');
%nc{'y_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
%nc{'y_u'}.long_name = ncchar('y location of U-points');
%nc{'y_u'}.units = ncchar('meter');

x_v_id = netcdf.defVar(id, 'x_v', 'double', [xi_v_dim eta_v_dim]);
netcdf.putAtt(id, x_v_id, 'long_name', 'x location of V-points');
netcdf.putAtt(id, x_v_id, 'units', 'meter');
%nc{'x_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
%nc{'x_v'}.long_name = ncchar('x location of V-points');
%nc{'x_v'}.units = ncchar('meter');

y_v_id = netcdf.defVar(id, 'y_v', 'double', [xi_v_dim eta_v_dim]);
netcdf.putAtt(id, y_v_id, 'long_name', 'y location of V-points');
netcdf.putAtt(id, y_v_id, 'units', 'meter');
%nc{'y_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
%nc{'y_v'}.long_name = ncchar('y location of V-points');
%nc{'y_v'}.units = ncchar('meter');

lat_rho_id = netcdf.defVar(id, 'lat_rho', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, lat_rho_id, 'long_name', 'latitude of RHO-points');
netcdf.putAtt(id, lat_rho_id, 'units', 'degree_north');
%nc{'lat_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'lat_rho'}.long_name = ncchar('latitude of RHO-points');
%nc{'lat_rho'}.units = ncchar('degree_north');

lon_rho_id = netcdf.defVar(id, 'lon_rho', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, lon_rho_id, 'long_name', 'longitude of RHO-points');
netcdf.putAtt(id, lon_rho_id, 'units', 'degree_east');
%nc{'lon_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'lon_rho'}.long_name = ncchar('longitude of RHO-points');
%nc{'lon_rho'}.units = ncchar('degree_east');

lat_psi_id = netcdf.defVar(id, 'lat_psi', 'double', [xi_psi_dim eta_psi_dim]);
netcdf.putAtt(id, lat_psi_id, 'long_name', 'latitude of PSI-points');
netcdf.putAtt(id, lat_psi_id, 'units', 'degree_north');
%nc{'lat_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
%nc{'lat_psi'}.long_name = ncchar('latitude of PSI-points');
%nc{'lat_psi'}.units = ncchar('degree_north');

lon_psi_id = netcdf.defVar(id, 'lon_psi', 'double', [xi_psi_dim eta_psi_dim]);
netcdf.putAtt(id, lon_psi_id, 'long_name', 'longitude of PSI-points');
netcdf.putAtt(id, lon_psi_id, 'units', 'degree_east');
%nc{'lon_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
%nc{'lon_psi'}.long_name = ncchar('longitude of PSI-points');
%nc{'lon_psi'}.units = ncchar('degree_east');

lat_u_id = netcdf.defVar(id, 'lat_u', 'double', [xi_u_dim eta_u_dim]);
netcdf.putAtt(id, lat_u_id, 'long_name', 'latitude of U-points');
netcdf.putAtt(id, lat_u_id, 'units', 'degree_north');
%nc{'lat_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
%nc{'lat_u'}.long_name = ncchar('latitude of U-points');
%nc{'lat_u'}.units = ncchar('degree_north');

lon_u_id = netcdf.defVar(id, 'lon_u', 'double', [xi_u_dim eta_u_dim]);
netcdf.putAtt(id, lon_u_id, 'long_name', 'longitude of U-points');
netcdf.putAtt(id, lon_u_id, 'units', 'degree_east');
%nc{'lon_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
%nc{'lon_u'}.long_name = ncchar('longitude of U-points');
%nc{'lon_u'}.units = ncchar('degree_east');

lat_v_id = netcdf.defVar(id, 'lat_v', 'double', [xi_v_dim eta_v_dim]);
netcdf.putAtt(id, lat_v_id, 'long_name', 'latitude of V-points');
netcdf.putAtt(id, lat_v_id, 'units', 'degree_north');
%nc{'lat_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
%nc{'lat_v'}.long_name = ncchar('latitude of V-points');
%nc{'lat_v'}.units = ncchar('degree_north');

lon_v_id = netcdf.defVar(id, 'lon_v', 'double', [xi_v_dim eta_v_dim]);
netcdf.putAtt(id, lon_v_id, 'long_name', 'longitude of V-points');
netcdf.putAtt(id, lon_v_id, 'units', 'degree_east');
%nc{'lon_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
%nc{'lon_v'}.long_name = ncchar('longitude of V-points');
%nc{'lon_v'}.units = ncchar('degree_east');

mask_rho_id = netcdf.defVar(id, 'mask_rho', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, mask_rho_id, 'long_name', 'mask on RHO-points');
netcdf.putAtt(id, mask_rho_id, 'flag_values', [0.0 1.0]);
netcdf.putAtt(id, mask_rho_id, 'flag_meanings', ['land', blanks(1), 'water']);
%nc{'mask_rho'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'mask_rho'}.long_name = ncchar('mask on RHO-points');
%nc{'mask_rho'}.option_0_ = ncchar('land');
%nc{'mask_rho'}.option_1_ = ncchar('water');
%name(nc{'mask_rho'}.option_0_, 'option(0)')
%name(nc{'mask_rho'}.option_1_, 'option(1)')

mask_zice_id = netcdf.defVar(id, 'mask_zice', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, mask_zice_id, 'long_name', 'mask on ice-shelf RHO-points');
netcdf.putAtt(id, mask_zice_id, 'flag_values', [0.0 1.0]);
netcdf.putAtt(id, mask_zice_id, 'flag_meanings', ['not ice', blanks(1), 'ice']);
%nc{'mask_zice'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'mask_zice'}.long_name = ncchar('mask on ice-shelf RHO-points');
%nc{'mask_zice'}.option_0_ = ncchar('land');
%nc{'mask_zice'}.option_1_ = ncchar('water');
%name(nc{'mask_zice'}.option_0_, 'option(0)')
%name(nc{'mask_zice'}.option_1_, 'option(1)')

mask_u_id = netcdf.defVar(id, 'mask_u', 'double', [xi_u_dim eta_u_dim]);
netcdf.putAtt(id, mask_u_id, 'long_name', 'mask on U-points');
netcdf.putAtt(id, mask_u_id, 'flag_values', [0.0 1.0]);
netcdf.putAtt(id, mask_u_id, 'flag_meanings', ['land', blanks(1), 'water']);
%nc{'mask_u'} = ncdouble('eta_u', 'xi_u'); %% 16770 elements.
%nc{'mask_u'}.long_name = ncchar('mask on U-points');
%nc{'mask_u'}.option_0_ = ncchar('land');
%nc{'mask_u'}.option_1_ = ncchar('water');
%name(nc{'mask_u'}.option_0_, 'option(0)')
%name(nc{'mask_u'}.option_1_, 'option(1)')
%                nc{'mask_u'}.FillValue_ = ncdouble(1);

mask_v_id = netcdf.defVar(id, 'mask_v', 'double', [xi_v_dim eta_v_dim]);
netcdf.putAtt(id, mask_v_id, 'long_name', 'mask on V-points');
netcdf.putAtt(id, mask_v_id, 'flag_values', [0.0 1.0]);
netcdf.putAtt(id, mask_v_id, 'flag_meanings', ['land', blanks(1), 'water']);
%nc{'mask_v'} = ncdouble('eta_v', 'xi_v'); %% 16770 elements.
%nc{'mask_v'}.long_name = ncchar('mask on V-points');
%nc{'mask_v'}.option_0_ = ncchar('land');
%nc{'mask_v'}.option_1_ = ncchar('water');
%name(nc{'mask_v'}.option_0_, 'option(0)')
%name(nc{'mask_v'}.option_1_, 'option(1)')
%                nc{'mask_v'}.FillValue_ = ncdouble(1);

mask_psi_id = netcdf.defVar(id, 'mask_psi', 'double', [xi_psi_dim eta_psi_dim]);
netcdf.putAtt(id, mask_psi_id, 'long_name', 'mask on PSI-points');
netcdf.putAtt(id, mask_psi_id, 'flag_values', [0.0 1.0]);
netcdf.putAtt(id, mask_psi_id, 'flag_meanings', ['land', blanks(1), 'water']);
%nc{'mask_psi'} = ncdouble('eta_psi', 'xi_psi'); %% 16641 elements.
%nc{'mask_psi'}.long_name = ncchar('mask on PSI-points');
%nc{'mask_psi'}.option_0_ = ncchar('land');
%nc{'mask_psi'}.option_1_ = ncchar('water');
%name(nc{'mask_psi'}.option_0_, 'option(0)')
%name(nc{'mask_psi'}.option_1_, 'option(1)')
%                nc{'mask_psi'}.FillValue_ = ncdouble(1);

% Now, what about depths: "h" and "hraw".  <== DEPTHS.

% The following seems mistaken:

angle_id = netcdf.defVar(id, 'angle', 'double', [xi_rho_dim eta_rho_dim]);
netcdf.putAtt(id, angle_id, 'long_name', 'angle between xi axis and east');
netcdf.putAtt(id, angle_id, 'units', 'degree');
%nc{'angle'} = ncdouble('eta_rho', 'xi_rho'); %% 16900 elements.
%nc{'angle'}.long_name = ncchar('angle between xi axis and east');
%nc{'angle'}.units = ncchar('degree');

netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', title);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS grid file');
netcdf.endDef(id);



%% Fill the variables with data.

disp(' ## Filling Variables...')
switch lower(projection)
    case 'mercator'
        theProjection = 'ME';
    case 'stereographic'
        theProjection = 'ST';
    case 'lambert conformal conic'
        theProjection = 'LC';
end

% Fill the variables.

% Need (x..., y...) in meters.  Currently, they are
%  in Mercator (projected) units.

netcdf.putVar(id, JPRJ_id, theProjection);
netcdf.putVar(id, spherical_id, 'T');
%nc{'JPRJ'}(:) = theProjection;
%nc{'spherical'}(:) = 'T';   % T or F -- uppercase okay?

netcdf.putVar(id, xl_id, xl);
netcdf.putVar(id, el_id, el);
%nc{'xl'}(:) = xl;
%nc{'el'}(:) = el;

f = 2.*7.29e-5.*sin(geogrid_lat(j_rho, i_rho).*pi./180);
netcdf.putVar(id, f_id, f);
%nc{'f'}(1:MP,1:LP) = f;

netcdf.putVar(id, x_rho_id, grid_x(j_rho, i_rho)');
netcdf.putVar(id, y_rho_id, grid_y(j_rho, i_rho)');
%nc{'x_rho'}(1:MP,1:LP) = grid_x(j_rho, i_rho);
%nc{'y_rho'}(1:MP,1:LP) = grid_y(j_rho, i_rho);

netcdf.putVar(id, x_psi_id, grid_x(j_psi, i_psi)');
netcdf.putVar(id, y_psi_id, grid_y(j_psi, i_psi)');
%nc{'x_psi'}(1:M,1:L) = grid_x(j_psi, i_psi);
%nc{'y_psi'}(1:M,1:L) = grid_y(j_psi, i_psi);

netcdf.putVar(id, x_u_id, grid_x(j_u, i_u)');
netcdf.putVar(id, y_u_id, grid_y(j_u, i_u)');
%nc{'x_u'}(1:MP,1:L) = grid_x(j_u, i_u);
%nc{'y_u'}(1:MP,1:L) = grid_y(j_u, i_u);

netcdf.putVar(id, x_v_id, grid_x(j_v, i_v)');
netcdf.putVar(id, y_v_id, grid_y(j_v, i_v)');
%nc{'x_v'}(1:M,1:LP) = grid_x(j_v, i_v);
%nc{'y_v'}(1:M,1:LP) = grid_y(j_v, i_v);

netcdf.putVar(id, lon_rho_id, geogrid_lon(j_rho, i_rho)');
netcdf.putVar(id, lat_rho_id, geogrid_lat(j_rho, i_rho)');
%nc{'lon_rho'}(1:MP,1:LP) = geogrid_lon(j_rho, i_rho);
%nc{'lat_rho'}(1:MP,1:LP) = geogrid_lat(j_rho, i_rho);

netcdf.putVar(id, lon_psi_id, geogrid_lon(j_psi, i_psi)');
netcdf.putVar(id, lat_psi_id, geogrid_lat(j_psi, i_psi)');
%nc{'lon_psi'}(1:M,1:L) = geogrid_lon(j_psi, i_psi);
%nc{'lat_psi'}(1:M,1:L) = geogrid_lat(j_psi, i_psi);

netcdf.putVar(id, lon_u_id, geogrid_lon(j_u, i_u)');
netcdf.putVar(id, lat_u_id, geogrid_lat(j_u, i_u)');
%nc{'lon_u'}(1:MP,1:L) = geogrid_lon(j_u, i_u);
%nc{'lat_u'}(1:MP,1:L) = geogrid_lat(j_u, i_u);

netcdf.putVar(id, lon_v_id, geogrid_lon(j_v, i_v)');
netcdf.putVar(id, lat_v_id, geogrid_lat(j_v, i_v)');
%nc{'lon_v'}(1:M,1:LP) = geogrid_lon(j_v, i_v);
%nc{'lat_v'}(1:M,1:LP) = geogrid_lat(j_v, i_v);

% Masking.

mask = ~raw_s.mw;
land = mask;
water = ~land;
imask = ~~raw_s.mi;
rmask = water;

imask(rmask == 0) = 0;

%ggfi(imask == 1 | rmask == 0) = 0;

ice_draft2 = ice_draft; % = ice_draft - ggfi;  < --------------------------------- I commented these: DEG.

imask(ice_draft2 < 0) = 1;
imask = ~~imask;

%bathymetry(bathymetry+ice_draft2 < 20) = -ice_draft2(bathymetry+ice_draft2 < 20)+20;

if ~isempty(bathymetry)
netcdf.putVar(id, h_id, bathymetry');
end
if ~isempty(ice_draft)
netcdf.putVar(id, zice_id, ice_draft2');
end
% Calculate other masking arrays.

umask = zeros(size(rmask));
vmask = zeros(size(rmask));
pmask = zeros(size(rmask));

for i = 2:LP
    for j = 1:MP
        umask(j, i-1) = rmask(j, i) * rmask(j, i-1);
    end
end

for i = 1:LP
    for j = 2:MP
        vmask(j-1, i) = rmask(j, i) * rmask(j-1, i);
    end
end
for i = 2:LP
    for j = 2:MP
        pmask(j-1, i-1) = rmask(j, i) * rmask(j, i-1) * rmask(j-1, i) * rmask(j-1, i-1);
    end
end

netcdf.putVar(id, mask_rho_id, double(rmask'));
netcdf.putVar(id, mask_psi_id, double(pmask(1:end-1, 1:end-1)'));
netcdf.putVar(id, mask_u_id, double(umask(1:end, 1:end-1)'));
netcdf.putVar(id, mask_v_id, double(vmask(1:end-1, 1:end)'));
netcdf.putVar(id, mask_zice_id, double(imask'));
%nc{'mask_rho'}(1:MP,1:LP) = rmask;
%nc{'mask_psi'}(1:M,1:L) = pmask(1:end-1, 1:end-1);
%nc{'mask_u'}(1:MP,1:L) = umask(1:end, 1:end-1);
%nc{'mask_v'}(1:M,1:LP) = vmask(1:end-1, 1:end);
%nc{'mask_zice'}(1:MP,1:LP) = imask;
% Average angle -- We should do this via (x, y) components.

% temp = 0.5*(pg.cosrx(1:end-1, :) + pg.cosrx(2:end, :));
% temp = 0.5*(temp(:,1:end-1) + temp(:,2:end));
% ang = zeros(n, m);
% ang(2:2:end,2:2:end) = temp(1:end-1,1:end-1);

%temp = 0.5*(pg.cosrx(1:end-1, :) + pg.cosrx(2:end, :));
%temp = 0.5*(temp(:,1:end-1) + temp(:,2:end));
ang = zeros(n, m);
ang(2:2:end,2:2:end) = s.angle(1:end-1,1:end-1);

netcdf.putVar(id, angle_id, ang(j_rho, i_rho)');
%nc{'angle'}(1:MP,1:LP) = ang(j_rho, i_rho);

% if (0)
%                sx = abs(spaced_x(:, 2:end) - spaced_x(:, 1:end-1));
%                sy = abs(spaced_x(2:end, :) - spaced_x(1:end-1, :));
% elseif (0)
%                sx = abs(spaced_x(1:2:end, 3:2:end) - spaced_x(1:2:end, 1:2:end-2));
%                sy = abs(spaced_y(3:2:end, 1:2:end) - spaced_y(1:2:end-2, 1:2:end));
%                sx = 0.5 * (sx(1:end-1, :) + sx(2:end, :));
%                sy = 0.5 * (sy(:, 1:end-1) + sy(:, 2:end));
% end

% Use geometry from seagrid file.
% Note: need half the number of points.

gx = geometry{1};   % Spherical distances in meters.

gy = geometry{2};

% raw_grid_size = [m, n]

% geometry_sizes = [size(gx) size(gy)]

% sx = 0.5*(gx(1:end-1, :) + gx(2:end, :));
% sy = 0.5*(gy(:, 1:end-1) + gy(:, 2:end));

% raw_s_sizes = [size(sx) size(sy)]

% sx = sx(2:end-1, :);
% sy = sy(:, 2:end-1);

pm = 1 ./ gx;
pn = 1 ./ gy;


% Look for CFL violations

dx =  min(gx,gy);
grav = 9.81;

%WCT = (pg.h1+pg.zice); %use wct and mask!!!!
%wct(:,1)=[];wct(:,end)=[];
%wct(1,:)=[];wct(end,:)=[];

dt = dx./ sqrt(grav*wct);
if plots
figure('renderer','zbuffer'),flat_pcolor(dt.*~mask)
min_dt = min(dt(:));
[min_dt_j,min_dt_i] = find(dt==min_dt);
hold on, plot(min_dt_i,min_dt_j,'ks')
str1 = num2str(min_dt);
str=(['Minimum dt =',str1]);
hold on, text(min_dt_i+5,min_dt_j,str)
end

%


netcdf.putVar(id, pm_id, pm');
netcdf.putVar(id, pn_id, pn');
%nc{'pm'}(1:MP,1:LP) = pm;
%nc{'pn'}(1:MP,1:LP) = pn;

dmde = zeros(size(pm));
dndx = zeros(size(pn));

dmde(2:end-1, :) = 0.5*(1./pm(3:end, :) - 1./pm(1:end-2, :));
dndx(:, 2:end-1) = 0.5*(1./pn(:, 3:end) - 1./pn(:, 1:end-2));

netcdf.putVar(id, dmde_id, dmde');
netcdf.putVar(id, dndx_id, dndx');
%nc{'dmde'}(1:MP,1:LP) = dmde;
%nc{'dndx'}(1:MP,1:LP) = dndx;

netcdf.putVar(id, depthmax_id, max_depth);
netcdf.putVar(id, depthmin_id, min_depth);
%nc{'depthmax'}(:) = max_depth;
%nc{'depthmin'}(:) = min_depth;
% Final size of file:

%s = size(nc);
%disp([' ## Dimensions: ' int2str(s(1))])
%disp([' ## Variables: ' int2str(s(2))])
%disp([' ## Global Attributes: ' int2str(s(3))])
%disp([' ## Record Dimension: ' name(recdim(nc))])
netcdf.close(id);
disp('Done.')
%endef(nc)
%close(nc)
