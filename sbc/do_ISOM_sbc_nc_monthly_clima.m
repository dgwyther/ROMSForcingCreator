%function  do_tisom001_sbc_nc(frcname,grdname,title,time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do_mer015_surf_nc('misom015_sbc.nc','mer012_grd.nc','MISOM surface heat/salt fluxes and wind stress',[1:12*16]);
%       frcname: name of the forcing file
%       grdname: name of the grid file
%       title: title in the netcdf file  
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grdname = '/ds/projects/iomp/totten/ana/dgwyther/netcdf/grid/tisom004_grd.nc';
% frcname = '/ds/projects/iomp/totten/ana/dgwyther/netcdf/iomp_IAF/sbc/tisom004_sbc.nc'; % Passed on by make_sbc.m
%time = [1:(MaxYear-MinYear+1)*12]%[1:192]; % This should be worked out by data given by make_sbc.m

title =([RunName,'surface heat/salt fluxes and wind stress']); % Make title from model name variable

disp('load air sea flux file (big file so load it once and subset later)')
load([RunName,'_air_sea_fluxes_monthly.mat'])
load(['ustress_grid_model.mat'])
load(['vstress_grid_model.mat'])

%frcname = [frcname(1:end-6),'monthly_cyclic_strongerwinter_sbc.nc'];

LeapYears = [1992:4:2040]; %leap years til 2040
%time_big=0;
%for YearInd = MinYear:MaxYear;
  % YearInd
  %    clear nw nc u_stress v_stress shfluxGrid_annual ssfluxGrid_annual shf shft swf swft result
  % nameval=['/ds/projects/iomp/totten/ana/dgwyther/netcdf/tisom007/sbc/tisom007_sbc_' num2str(YearInd) '.nc'];
  %    frcname = nameval;

% if any(YearInd == LeapYears)
%  time = 366;
% else
%  time = 365;
% end
%time_big=time_big+time;
%end
%no_days = time_big;
addpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib/'))
addpath('/ds/projects/iomp/matlab_scripts')
time_big = 12; %12*(MaxYear-MinYear+1);%round(time_big./30.4375);
nc=netcdf(grdname);
L=length(nc('xi_psi'));
M=length(nc('eta_psi'));
result=close(nc);
Lp=L+1;
Mp=M+1;
cycle = [30.4375*12];%*(MaxYear-MinYear+1)];
%cycle_daily   = [30.4375*12*(MaxYear-MinYear+1)];
nw = netcdf(frcname, 'clobber');
result = redef(nw);



%
%  Create dimensions
%
    disp('Create dimension')

nw('xi_u') = L;
nw('eta_u') = Mp;
nw('xi_v') = Lp;
nw('eta_v') = M;
nw('xi_rho') = Lp;
nw('eta_rho') = Mp;
nw('xi_psi') = L;
nw('eta_psi') = M;
nw('sms_time') = time_big;%length(time);
nw('shf_time') = time_big;%length(time);
nw('swf_time') = time_big;%length(time);
%nw('sst_time') = length(time);
%nw('srf_time') = length(time);
%nw('sss_time') = length(time);
%
%  Create variables and attributes
%
    disp('Create variable attribute')

nw{'sms_time'} = ncdouble('sms_time');
nw{'sms_time'}.long_name = ncchar('surface momentum stress time');
nw{'sms_time'}.long_name = 'surface momentum stress time';
nw{'sms_time'}.units = ncchar('days');
nw{'sms_time'}.units = 'days';
nw{'sms_time'}.cycle_length = cycle;

nw{'shf_time'} = ncdouble('shf_time');
nw{'shf_time'}.long_name = ncchar('surface heat flux time');
nw{'shf_time'}.long_name = 'surface heat flux time';
nw{'shf_time'}.units = ncchar('days');
nw{'shf_time'}.units = 'days';
nw{'shf_time'}.cycle_length = cycle;

nw{'swf_time'} = ncdouble('swf_time');
nw{'swf_time'}.long_name = ncchar('surface freshwater flux time');
nw{'swf_time'}.long_name = 'surface freshwater flux time';
nw{'swf_time'}.units = ncchar('days');
nw{'swf_time'}.units = 'days';
nw{'swf_time'}.cycle_length = cycle;

%{
nw{'sst_time'} = ncdouble('sst_time');
nw{'sst_time'}.long_name = ncchar('sea surface temperature time');
nw{'sst_time'}.long_name = 'sea surface temperature time';
nw{'sst_time'}.units = ncchar('days');
nw{'sst_time'}.units = 'days';
nw{'sst_time'}.cycle_length = cycle;

nw{'sss_time'} = ncdouble('sss_time');
nw{'sss_time'}.long_name = ncchar('sea surface salinity time');
nw{'sss_time'}.long_name = 'sea surface salinity time';
nw{'sss_time'}.units = ncchar('days');
nw{'sss_time'}.units = 'days';
nw{'sss_time'}.cycle_length = cycle;

nw{'srf_time'} = ncdouble('srf_time');
nw{'srf_time'}.long_name = ncchar('solar shortwave radiation time');
nw{'srf_time'}.long_name = 'solar shortwave radiation time';
nw{'srf_time'}.units = ncchar('days');
nw{'srf_time'}.units = 'days';
nw{'srf_time'}.cycle_length = cycle;
%}

nw{'sustr'} = ncdouble('sms_time', 'eta_u', 'xi_u');
nw{'sustr'}.long_name = ncchar('surface u-momentum stress');
nw{'sustr'}.long_name = 'surface u-momentum stress';
nw{'sustr'}.units = ncchar('Newton meter-2');
nw{'sustr'}.units = 'Newton meter-2';

nw{'svstr'} = ncdouble('sms_time', 'eta_v', 'xi_v');
nw{'svstr'}.long_name = ncchar('surface v-momentum stress');
nw{'svstr'}.long_name = 'surface v-momentum stress';
nw{'svstr'}.units = ncchar('Newton meter-2');
nw{'svstr'}.units = 'Newton meter-2';

nw{'shflux'} = ncdouble('shf_time', 'eta_rho', 'xi_rho');
nw{'shflux'}.long_name = ncchar('surface net heat flux');
nw{'shflux'}.long_name = 'surface net heat flux';
nw{'shflux'}.units = ncchar('Watts meter-2');
nw{'shflux'}.units = 'Watts meter-2';

nw{'swflux'} = ncdouble('swf_time', 'eta_rho', 'xi_rho');
nw{'swflux'}.long_name = ncchar('surface freshwater flux (E-P)');
nw{'swflux'}.long_name = 'surface freshwater flux (E-P)';
nw{'swflux'}.units = ncchar('centimeter day-1');
nw{'swflux'}.units = 'centimeter day-1';
nw{'swflux'}.positive = ncchar('net evaporation');
nw{'swflux'}.positive = 'net evaporation';
nw{'swflux'}.negative = ncchar('net precipitation');
nw{'swflux'}.negative = 'net precipitation';

nw{'SST'} = ncdouble('sst_time', 'eta_rho', 'xi_rho');
nw{'SST'}.long_name = ncchar('sea surface temperature');
nw{'SST'}.long_name = 'sea surface temperature';
nw{'SST'}.units = ncchar('Celsius');
nw{'SST'}.units = 'Celsius';

nw{'SSS'} = ncdouble('sss_time', 'eta_rho', 'xi_rho');
nw{'SSS'}.long_name = ncchar('sea surface salinity');
nw{'SSS'}.long_name = 'sea surface salinity';
nw{'SSS'}.units = ncchar('PSU');
nw{'SSS'}.units = 'PSU';

nw{'dQdSST'} = ncdouble('sst_time', 'eta_rho', 'xi_rho');
nw{'dQdSST'}.long_name = ncchar('surface net heat flux sensitivity to SST');
nw{'dQdSST'}.long_name = 'surface net heat flux sensitivity to SST';
nw{'dQdSST'}.units = ncchar('Watts meter-2 Celsius-1');
nw{'dQdSST'}.units = 'Watts meter-2 Celsius-1';

nw{'swrad'} = ncdouble('srf_time', 'eta_rho', 'xi_rho');
nw{'swrad'}.long_name = ncchar('solar shortwave radiation');
nw{'swrad'}.long_name = 'solar shortwave radiation';
nw{'swrad'}.units = ncchar('Watts meter-2');
nw{'swrad'}.units = 'Watts meter-2';
nw{'swrad'}.positive = ncchar('downward flux, heating');
nw{'swrad'}.positive = 'downward flux, heating';
nw{'swrad'}.negative = ncchar('upward flux, cooling');
nw{'swrad'}.negative = 'upward flux, cooling';

result = endef(nw);

%
% Create global attributes
%
    disp('Create global attribute')

nw.title = ncchar(title);
nw.title = title;
nw.date = ncchar(date);
nw.date = date;
nw.grd_file = ncchar(grdname);
nw.grd_file = grdname;
nw.type = ncchar('ROMS forcing file');
nw.type = 'ROMS forcing file';

%
% Write time variables
%
    disp('Write time variable')

%--- Quick fix for 1992-2012 takeshi -> yearly takeshi
%shfluxGrid_annual = shfluxGrid(Takeshi_index:Takeshi_index+time-1,:,:);
%ssfluxGrid_annual = ssfluxGrid(Takeshi_index:Takeshi_index+time-1,:,:);
%---

swft = [(30.4375/2):30.4375:(30.4375*12)-(30.4375/2)]; %data at mid each month
shft = [(30.4375/2):30.4375:(30.4375*12)-(30.4375/2)];
%swft = [1:time_big]; %data at each day
%shft = [1:time_big];


%--- For loading all data at once (big files for daily forcing!)
%load([WorkingDir,'ustress_grid_model_',num2str(MinYear),'_',num2str(MaxYear),'.mat'])
%load([WorkingDir,'vstress_grid_model_',num2str(MinYear),'_',num2str(MaxYear),'.mat'])
%urot_stress=u_stress_All;
%vrot_stress=v_stress_All;

%--- For loading year at a time
%load([WorkingDir,'ustress_grid_model_',num2str(YearInd),'.mat'])
%load([WorkingDir,'vstress_grid_model_',num2str(YearInd),'.mat'])


%
% Write time variables
%
%disp('interpolating daily -> monthly')

%dpy    = [31 28 31 30 31 30 31 31 30 31 30 31];
%dpy_ly = [31 29 31 30 31 30 31 31 30 31 30 31];

%days = [dpy_ly,dpy,dpy,dpy,dpy_ly,dpy,dpy,dpy,dpy_ly,dpy,dpy,dpy,dpy_ly,dpy,dpy,dpy,dpy_ly,dpy,dpy,dpy,dpy_ly]; %%%% THIS IS TOTALLY SPECIFIC TO THIS RUN!!!!!!!

%someind=1;
%shf_n = nan([time_big,size(shfluxGrid,2),size(shfluxGrid,3)]);
%swf_n = nan([time_big,size(shfluxGrid,2),size(shfluxGrid,3)]);
%su_n  = nan([time_big,size(shfluxGrid,2),size(shfluxGrid,3)]);
%sv_n  = nan([time_big,size(shfluxGrid,2),size(shfluxGrid,3)]);


%for jj = 1:time_big
%shf_n(jj,:,:) = squeeze(nanmean(shfluxGrid(someind:someind+days(jj)-1,:,:),1));
%swf_n(jj,:,:) = squeeze(nanmean(ssfluxGrid(someind:someind+days(jj)-1,:,:),1));
%su_n(jj,:,:) = squeeze(nanmean(u_stress_All(someind:someind+days(jj)-1,:,:),1));
%sv_n(jj,:,:) = squeeze(nanmean(v_stress_All(someind:someind+days(jj)-1,:,:),1));

%someind=someind+days(jj);
%disp([num2str(jj/time_big*100),'%'])
%end

disp('calculate climatology from monthly data')

%mon_inds = 1:12*(MaxYear-MinYear+1);
%mon_inds = reshape(mon_inds,12,[]);

shfluxMon = reshape(shfluxGrid,12,[],size(shfluxGrid,2),size(shfluxGrid,3));
ssfluxMon = reshape(ssfluxGrid,12,[],size(ssfluxGrid,2),size(ssfluxGrid,3));
u_stress_Mon = reshape(u_stress_All,12,[],size(u_stress_All,2),size(u_stress_All,3));
v_stress_Mon = reshape(v_stress_All,12,[],size(v_stress_All,2),size(v_stress_All,3));

shf = squeeze(nanmean(shfluxMon,2));
swf = squeeze(nanmean(ssfluxMon,2));
su  = squeeze(nanmean(u_stress_Mon,2));
sv  = squeeze(nanmean(v_stress_Mon,2));

disp('inpainting nans in wind fields...')

for j = 1:time_big;
  shf(j,:,:)= inpaint_nans(squeeze( shf(j,:,:)),1);
  swf(j,:,:)= inpaint_nans(squeeze( swf(j,:,:)),1);
  su(j,:,:) = inpaint_nans(squeeze(su(j,:,:)),1);
  sv(j,:,:) = inpaint_nans(squeeze(sv(j,:,:)),1);
disp([num2str(j/time_big*100),'%'])
end

% Correct Heat flux for large summer values (as per NCEP2 heat flux data)

shf(shf > 0) = shf(shf > 0)*0.5;
%swf(swf <=0) = -1e-6;

nw{'sms_time'}(:) = swft;
nw{'shf_time'}(:) = shft;
nw{'swf_time'}(:) = swft;
%nw{'sst_time'}(:) = shft;
%nw{'srf_time'}(:) = shft;
%nw{'sss_time'}(:) = shft;

% Write salt and heat flux variables
    disp('write salt and heat flux variables')

refSalt = 34.4; % Reference salinity
dm = 365.25/12; % Days in a month (365.25days/12 months):

nw{'shflux'}(:,:,:) = shf; % (W/m^2)
nw{'swflux'}(:,:,:) = swf./refSalt/dm*100; % convert swf kg/m2/month to (cm/m2/day)

    disp('Write wind variables')

nw{'sustr'}(:,:,:) = (su(:,:,1:end-1)+su(:,:,2:end))./2;
nw{'svstr'}(:,:,:) = (sv(:,1:end-1,:)+sv(:,2:end,:))./2;

close(nw);
    disp('DONE')

rmpath(genpath('/ds/projects/iomp/matlab_scripts/netcdflib'))
rmpath('/ds/projects/iomp/matlab_scripts')
disp('ALL DONE! YAY')

