% make_sbc.m
% David Gwyther | 2011
% updated heavily in Dec 2015

%addpath('/ds/projects/iomp/matlab_scripts')
%addpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib')) %script to read ncep data

 
grdname = '/ds/projects/iomp/totten/ana/dgwyther/tisom009/grid/tisom008_canal_grd.nc'; 
frcname = 'tisom009_dailyclima_sbc.nc';
RunName = 'tisom009';
MinYear = 1992;
MaxYear = 2013;
windbounds = [13 124 20 57]; %lonmin lonmax latmin latmax
Takeshibounds = [410 550 500 640]; %[colmin colmax rowmin rowmax]
% new Amery flux bounds  [100 352 520 720];
% old Totten flux bounds [410 550 500 640];
% new Totten flux bounds [380 550 500 640];
% Totten ERA bounds [13 124 20 57];
% Mertz ERA bounds [86 108 101 107];
% old Totten CORE bounds [60 80 10 18];
% new Amery CORE bounds [29 49 7 18];
FluxType = 1;
% (1) daily fluxes (1992-2013)
% (2) daily fluxes averaged to monthly (1992-2013)
% (3) daily fluxes climatology with TOK anomaly (XX - XX ?)
% (4) monthly fluxes (1992-2007 only)

WindType = 3;
% (1) COREv2 daily->monthly (1992-2007)
% (2) COREv2 daily (1992-2014)
% (3) COREv1 normal year
% (4) ERA-interim daily forcing (1992-2013)
% (5) ERA-interim daily->monthly forcing (1992-2013)

ForcingType = 6;
% (1) Interannual forcing Daily, multiple files
% (2) Interannual forcing Daily, 1 file
% (3) Interannual forcing, Monthly
% (4) Climatology forcing, Monthly
% (5) Climatology forcing (364 days), Monthly
% (6) Climatology forcing (364 days), Daily
% (7) Constant forcing
% (8) Constant forcing (summer)

%%%%% DO NOT EDIT BELOW %%%%%%%%
xmin = Takeshibounds(3);
xmax = Takeshibounds(4);
ymin = Takeshibounds(1);
ymax = Takeshibounds(2);
if FluxType == 1
read_tamura_daily
elseif FluxType == 2
read_tamura_daily_to_monthly
elseif FluxType == 3
read_tamura_daily_clima_TOK
elseif FluxType == 4
read_tamura_monthly
end

Imin_wind = windbounds(1);
Imax_wind = windbounds(2);
Jmin_wind = windbounds(3);
Jmax_wind = windbounds(4);
if WindType == 1;
CORE_grid_stress_daily_to_monthly
elseif WindType == 2;
CORE_stress_daily
elseif WindType == 3
COREv1_stress_normal
elseif WindType == 4;
ERA_interim_misom_grid_stress_annual
elseif WindType == 5;
ERA_interim_misom_grid_stress_annual_daily_to_monthly
end

if ForcingType == 1
do_ISOM_sbc_nc_daily_v2
elseif ForcingType == 2
do_ISOM_sbc_nc_daily_1file_v1
elseif ForcingType == 3
do_ISOM_sbc_nc_monthly_1file_v2
elseif ForcingType == 4
do_ISOM_sbc_nc_monthly_clima
elseif ForcingType == 5
do_ISOM_sbc_nc_monthly_clima_364DayYr
elseif ForcingType == 6
do_ISOM_sbc_nc_daily_clima_364DayYr
elseif ForcingType == 7
do_ISOM_sbc_nc_monthly_AvState
elseif ForcingType == 8
do_ISOM_sbc_nc_monthly_AvState_summer
end

