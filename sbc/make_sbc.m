% make_sbc.m
% David Gwyther | 2011
% updated heavily in Dec 2015

%addpath('/ds/projects/iomp/matlab_scripts')
%addpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib')) %script to read ncep data

 
grdname = '../grid/tisom008_canal_grd.nc'; 
frcname = 'tisom016_1992IAF_sbc.nc';
RunName = 'tisom016';
MinYear = 1992;
MaxYear = 1992;
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
% (1) daily fluxes (1992-2013) (automatically remove leap year data)
% (2) daily fluxes averaged to monthly (1992-2013)
% (3) daily fluxes 1995 
% (4) daily fluxes climatology + TOK COREv2 anomaly (XX - XX ?)
% (5) monthly fluxes (1992-2007 only)
% (6) daily fluxes with leap year data (1992-2013). Be careful converting to climatology

WindType = 4;
% (1) COREv2 daily->monthly (1992-2007)
% (2) COREv2 daily (1992-2014)
% (3) COREv1 normal year (6hrly, 365 days)
% (4) ERA-interim daily forcing (1992-2013)
% (5) ERA-interim daily->monthly forcing (1992-2013)
% (6) ERA-interim daily forcing (1992-2013) with leap year data
% (7) ERA-interim daily forcing (1995)
% (8) ERA-interim daily forcing (1992-2017) with sea ice masking

ForcingType = 3;
% (1) Interannual forcing Daily, multiple files (1 per year)
% (2) Interannual forcing Daily, 1 file
% (3) Interannual forcing Daily, multiple files (1 per forcing type)
% (4) Interannual forcing, Monthly
% (5) Climatology forcing, Monthly  % DEPRECATED. REMOVE. %%%%%%%%%%%%%%%%%%%%%%%
% (6) Climatology forcing (364 days), Monthly
% (7) Climatology forcing (364 days), Daily
% (8) Climatology forcing (365 days), Monthly
% (9) Climatology forcing (365 days), Daily
% (10) Constant forcing
% (11) Constant forcing (summer)

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
read_tamura_1995
elseif FluxType == 4
read_tamura_daily,
read_okane_daily,
do_tamuraclima_plus_okane
elseif FluxType == 5
read_tamura_monthly
elseif FluxType == 6
read_tamura_daily_withleapyears
end

Imin_wind = windbounds(1);
Imax_wind = windbounds(2);
Jmin_wind = windbounds(3);
Jmax_wind = windbounds(4);
if WindType == 1;
COREv2_IAF_stress_daily_to_monthly
elseif WindType == 2;
COREv2_IAF_stress_daily
elseif WindType == 3
COREv1_normalyear_stress_daily
elseif WindType == 4;
ERA_interim_misom_grid_stress_annual
elseif WindType == 5;
ERA_interim_misom_grid_stress_annual_daily_to_monthly
elseif WindType == 6;
ERA_interim_misom_grid_stress_annual_withleapyears
elseif WindType == 7;
ERA_interim_misom_grid_stress_1995
elseif WindType == 8;
Tmin_ERAseaice = 1 % 1992-01-01
Tmax_ERAseaice = daysact('01-jan-1992',  '31-dec-2015')+1 % 2015-12-31 %sea ice data downloaded by DEG (early 2018) had 9497 timesteps: daysact('01-jan-1992',  '31-dec-2017')+1
StressCd='Lupkes'
ERA_interim_grid_stress_daily_seaicemask
end


if ForcingType == 1
do_ISOM_sbc_nc_daily_IAF_yearlyfiles
elseif ForcingType == 2
do_ISOM_sbc_nc_daily_IAF
elseif ForcingType == 3
do_ISOM_sbc_nc_daily_IAF_separatefrc
elseif ForcingType == 4
do_ISOM_sbc_nc_monthly_IAF
elseif ForcingType == 5
do_ISOM_sbc_nc_monthly_clima
elseif ForcingType == 6
do_ISOM_sbc_nc_monthly_clima_364DayYr
elseif ForcingType == 7
do_ISOM_sbc_nc_daily_clima_364DayYr
elseif ForcingType == 8
do_ISOM_sbc_nc_monthly_clima_365DayYr
elseif ForcingType == 9
do_ISOM_sbc_nc_daily_clima_365DayYr
elseif ForcingType == 10
do_ISOM_sbc_nc_monthly_AvState
elseif ForcingType == 11
do_ISOM_sbc_nc_monthly_AvState_summer
end

