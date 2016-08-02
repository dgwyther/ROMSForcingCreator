% make_lbc.m
% David Gwyther | 2014-5
% changelog:
% 2014:		 written
% 2015-Sep-08: 	 updated to remove 'interactive' mode. I don't EVER use it.
% 2015-Dec-01:   Updated to remove a bunch of redundant features

%addpath(genpath('/ds/projects/iomp/matlab_scripts'))
grdname = '/mnt/IceOceanVolume/aisom004/grid/aisom004_grd.nc'; %grid name
%bryname = '/ds/projects/iomp/aisom/ana/dgwyther/grid/aisom002/aisom002_bry.nc'; %output filename.
bryname = 'aisom004_cube92-3day_bry.nc';
MinYear = 1992;
MaxYear = 2013;
ECCObounds = [210 400 50 125];%as [xmin xmax ymin ymax]; 
RunName = 'aisom004'

Vtransform = 2;
Vstretching = 4;
theta_s = 4;
theta_b = 0.9;
Tcline = 20;
N = 31;

%%%%%%%%%%%%%%%%%%%
% Totten is [410 525 80 125]; %as [xmin xmax ymin ymax];
% Amery (new) is [210 400 50 125];

DataProduct = 5; 
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube84 
% (4) cube92 monthly
% (5) cube92 3-daily
% (6) Dinniman ACCIMA 5-km
% (7) O'Kane monthly ocean model (100 normal years) - COREv1 forced
% (8) O'Kane monthly ocean model (1948-2006) - COREv2 forced


ForcingType = 10;
% (1) cube84 repeated 1992 ] deprecated
% (2) cube84 interannual   ] deprecated
% (3) cube92 interannual
% (4) cube92 climatology
% (5) cube92 constant forcing
% (6) cube92 climatology w/ 364 day year
% (7) cube92 constant summer forcing
% (8) cube84 climatology
% (9) cube92-3day climatology w/ 364 day year
%(10) cube92-3day climatology
%(11) ACCIMA 5-km climatology
%(12) ACCIMA 5-km IAF
%(13) O'Kane-COREv1 monthly interannual
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

force_id = {'salt','theta','uvel','vvel'}; %load and interp ecco2 data
xmin = ECCObounds(1); 
xmax = ECCObounds(2);
ymin = ECCObounds(3); 
ymax = ECCObounds(4); 

if DataProduct == 1
    for ii = 1:4
        eval(['do_interp_ecco2_',force_id{ii},'_1yrforcing_A'])
    end
elseif DataProduct == 2 
    for ii = 1:4
        eval(['do_interp_ecco2_',force_id{ii},'_A'])
    end  
elseif DataProduct == 3
    for ii = 1:4
        eval(['do_load_ecco2_',force_id{ii},'_cube84'])
    end
elseif DataProduct == 4
    for ii = 1:4
        eval(['do_load_ecco2_',force_id{ii},'_cube92'])
    end
elseif DataProduct == 5
    for ii = 1:4
        eval(['do_load_ecco2_',force_id{ii},'_cube92_3day'])
    end
elseif DataProduct == 6
    for ii = 1:4
        eval(['do_load_ACCIMA_',force_id{ii},'_5km'])
    end
elseif DataProduct == 7
    for ii = 1:4
        eval(['do_load_OKane_',force_id{ii},'_corev1'])
    end
elseif DataProduct == 8
    for ii = 1:4
        eval(['do_load_OKane_',force_id{ii},'_corev2'])
    end

end

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
do_ISOM_lbc_cube92_clima
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
do_ISOM_lbc_ACCIMA_5km_clima
elseif ForcingType == 12;
do_ISOM_lbc_ACCIMA_5km
elseif ForcingType == 13;
do_ISOM_lbc_OKane_corev1_iaf

end



disp('done!')
    
