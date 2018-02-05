% make_lbc.m
% David Gwyther | 2014-5
% changelog:
% 2014:		 written
% 2015-Sep-08: 	 updated to remove 'interactive' mode. I don't EVER use it.
% 2015-Dec-01:   Updated to remove a bunch of redundant features

%addpath(genpath('/ds/projects/iomp/matlab_scripts'))
grdname = '../../tisom008_canal_grd.nc'; %grid name
%bryname = '/ds/projects/iomp/aisom/ana/dgwyther/grid/aisom002/aisom002_bry.nc'; %output filename.
bryname = 'tisom015_1992-2015clima_bry.nc';
MinYear = 1992;
MaxYear = 2015;
ECCObounds = [410 525 80 125];%as [xmin xmax ymin ymax]; 
RunName = 'tisom015'

Vtransform = 2;
Vstretching = 4;
theta_s = 4;
theta_b = 0.9;
Tcline = 20;
N = 31;

%%%%%%%%%%%%%%%%%%%
% Totten is [410 525 80 125]; %as [xmin xmax ymin ymax];
% Amery (new) is [210 400 50 125];

DataProduct = 'cube92 monthly from 3-daily';
%% options:
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube84 
% (4) cube92 monthly
% (5) cube92 monthly 1995
% (6) cube92 3-daily
% (7) cube92 monthly from 3-daily
% (8) Dinniman ACCIMA 5-km
% (9) O'Kane monthly ocean model (100 normal years) - COREv1 forced
% (10) O'Kane monthly ocean model (1948-2006) - COREv2 forced


ForcingType = 11;
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube92 interannual w/ 365 day year
% (4) cube92 climatology
% (5) cube92 constant forcing
% (6) cube92 climatology w/ 364 day year
% (7) cube92 constant summer forcing
% (8) cube84 climatology
% (9) cube92-3day climatology w/ 364 day year
%(10) cube92-3day climatology
%(11) cube92-3day(monthly avs) climatology w/ 365 day year
%(12) ACCIMA 5-km climatology
%(13) ACCIMA 5-km IAF
%(14) O'Kane-COREv1 monthly interannual (model yrs 1900-2000)
%(15) O'Kane-COREv1 monthly interannual climatology (model yrs 1900-2000)
%(16) O'Kane-COREv2 monthly interannual (1948-2006) - COREv2 forced

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

force_id = {'salt','theta','uvel','vvel'}; %load and interp ecco2 data
xmin = ECCObounds(1); 
xmax = ECCObounds(2);
ymin = ECCObounds(3); 
ymax = ECCObounds(4); 

if strcmp(DataProduct,'cube84 repeated 1992')
    for ii = 1:4
        eval(['do_interp_ecco2_',force_id{ii},'_1yrforcing_A'])
    end
elseif strcmp(DataProduct,'cube84 interannual')
    for ii = 1:4
        eval(['do_interp_ecco2_',force_id{ii},'_A'])
    end
elseif strcmp(DataProduct,'cube84')
    for ii = 1:4
        eval(['do_load_ecco2_',force_id{ii},'_cube84'])
    end
elseif strcmp(DataProduct,'cube92 monthly')
    for ii = 1:4
        disp(['loading ' force_id{ii} ' data'])
        eval(['do_load_ecco2_',force_id{ii},'_cube92'])
    end
elseif strcmp(DataProduct,'cube92 monthly 1995')
    for ii = 1:4
        disp(['loading ' force_id{ii} ' data'])
        eval(['do_load_ecco2_',force_id{ii},'_cube92_1Year'])
    end
elseif strcmp(DataProduct,'cube92 3-daily')
    for ii = 1:4
        disp(['loading ' force_id{ii} ' data'])
        eval(['do_load_ecco2_',force_id{ii},'_cube92_1Year'])
    end
elseif strcmp(DataProduct,'cube92 monthly from 3-daily')
    for ii = 1:4
        disp(['loading ' force_id{ii} ' data at ',datestr(now)])
        eval(['do_load_ecco2_',force_id{ii},'_cube92_3day_monthly'])
    end
elseif strcmp(DataProduct,'Dinniman ACCIMA 5km')
    for ii = 1:4
        eval(['do_load_ecco2_',force_id{ii},'_cube92_3day'])
    end
elseif strcmp(DataProduct,'O''Kane monthly ocean model COREv1')
    for ii = 1:4
        eval(['do_load_ACCIMA_',force_id{ii},'_5km'])
    end
elseif strcmp(DataProduct,'O''Kane monthly ocean model COREv2')
    for ii = 1:4
        eval(['do_load_OKane_',force_id{ii},'_corev1'])
    end
elseif DataProduct == 9
    for ii = 1:4
        eval(['do_load_OKane_',force_id{ii},'_corev2'])
    end
else
disp(['your choice of ',DataProduct,' doesn''t match the list'])
end
disp(['Your data product is ',DataProduct])

%% REGRID ECCO2 TO MODEL & run do_*isom_lbc.m
xmin = ECCObounds(1); 
xmax = ECCObounds(2);
ymin = ECCObounds(3); 
ymax = ECCObounds(4); 

if ForcingType <= 2
do_ISOM_lbc_A
elseif ForcingType == 3
do_ISOM_lbc_cube92
elseif ForcingType == 4
do_ISOM_lbc_cube92_clima THIS OPTION SEEMS TO BE MISSING A FILE?
elseif ForcingType == 5;
do_ISOM_lbc_cube92_AvState
elseif ForcingType == 6;
do_ISOM_lbc_cube92_clima_364DayYr
elseif ForcingType == 7;
do_ISOM_lbc_cube92_AvState_summer
elseif ForcingType == 8;
do_ISOM_lbc_cube84_clima
elseif ForcingType == 9;
do_ISOM_lbc_cube92_3day_clima_364DayYr
elseif ForcingType == 10;
do_ISOM_lbc_cube92_3day_clima
elseif ForcingType == 11;
do_ISOM_lbc_cube92_3day_monthly_clima
elseif ForcingType == 12;
do_ISOM_lbc_ACCIMA_5km_clima
elseif ForcingType == 13;
do_ISOM_lbc_ACCIMA_5km
elseif ForcingType == 14;
do_ISOM_lbc_OKane_corev1_iaf
elseif ForcingType == 15;
do_ISOM_lbc_OKane_corev1_iaf_clima
elseif ForcingType == 16;
do_ISOM_lbc_OKane_corev2_iaf
end



disp('done!')
    
