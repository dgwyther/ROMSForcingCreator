function do_ISOM_lbc_nc_A(bryname,grdname,title,obc,...
    theta_s,theta_b,hc,N,...
    time,cycle,clobber);

% % % For debugging use!
% % bryname ='tisom001_bry.nc';
% % grdname ='tisom001_grd.nc';
% % title = 'Lateral Boundaries Salt and Temp';
% % obc =[0 1 1 1];
% % theta_s =5;
% % theta_b =0.4;
% % hc =20;
% % N =31;
% % time = [(30.4375/2):30.4375:(30.4375*12*16)-(30.4375/2)];
% % cycle =[30.4375*12*16];
% % clobber ='clobber';
% % % For debugging! 

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function create_bryfile(bryname,grdname,title,obc...
%                          theta_s,theta_b,hc,N,...
%                          time,cycle,clobber);
%
%   This function create the header of a Netcdf climatology
%   file.
%
%   Input:
%
%   bryname      Netcdf climatology file name (character string).
%   grdname      Netcdf grid file name (character string).
%   obc          open boundaries flag (1=open , [S E N W]).
%   theta_s      S-coordinate surface control parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)
%   time         time.(vector)
%   cycle        Length (days) for cycling the climatology.(Real)
%   clobber      Switch to allow or not writing over an existing
%                file.(character string)
%
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
%load mer_ecco2.mat ecco2 %AIS_levitusb.mat
load isfc_AvState.mat
InterpSurface = isfc;
NST = 10;

% for i = 279:283
% InterpSurface.tmp(:,:,i) = InterpSurface.tmp(:,:,284);
% InterpSurface.sal(:,:,i) = InterpSurface.sal(:,:,284);
% end
% for i = 285:290
% InterpSurface.tmp(:,:,i) = InterpSurface.tmp(:,:,284);    
% InterpSurface.sal(:,:,i) = InterpSurface.sal(:,:,284);
% end

% for i = 1:12;
%     for j = 1:length(InterpSurface.tmp);
% InterpSurface.tmp(i,1,j) = 0;
% InterpSurface.sal(i,1,j) = 0;
%     end
% end

%tdays = [(30.4375/2):30.4375:(30.4375*12*16)-(30.4375/2)];

gamma = 0.0001;
hlfreeze = 334000.0;
refSalt = 34.4;
trelax = 30.0 * 86400.0;
saltMax = 35.0;
saltMin = 33.3;
sRateInc = 0.0085;
sRateDec = 0.0283333;
Tmax = -0.9;
Tmin = -1.9;
%tyear = mod(tdays,365.25);

% for i = 1:size(tdays',1)
%     if tyear(i) < 59.0
%         InterpSurface.tmp(i,1,:) = Tmax;
%         InterpSurface.sal(i,1,:) = saltMin;
%     elseif tyear(i) < 259.0
%         InterpSurface.tmp(i,1,:) = Tmin;
%         InterpSurface.sal(i,1,:) = saltMin + sRateInc * (tyear(i) - 59.0);
%     elseif tyear(i) < 319.0
%         InterpSurface.tmp(i,1,:) = Tmin;
%         InterpSurface.sal(i,1,:) = saltMax + sRateDec * (259.0 - tyear(i));
%     else
%         InterpSurface.tmp(i,1,:) = Tmax;
%         InterpSurface.sal(i,1,:) = saltMin;
%     end
% end

disp(' ')
disp([' Creating the file : ',bryname])
disp(' ')
%
%  Read the grid file and check the topography
%
nc = netcdf(grdname, 'nowrite');
h=nc{'h'}(:);
angle=nc{'angle'}(:);
maskr=nc{'mask_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
status=close(nc);
hmin=min(min(h(maskr==1)));
if hc > hmin
    error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
end
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the boundary file
%
type = 'BOUNDARY file' ;
history = 'ROMS' ;
nc = netcdf(bryname,clobber);
result = redef(nc);
%
%  Create dimensions
%
nc('xi_u') = L;
nc('xi_rho') = Lp;
nc('eta_v') = M;
nc('eta_rho') = Mp;
nc('s_rho') = N;
nc('temp_time') = length(time);
nc('salt_time') = length(time);
 for i = 1:NST;
        if i < 10
            eval(['nc(''sand_0' num2str(i) '_time'') = ' num2str(length(time)) ';']);
        else
            eval(['nc(''sand_' num2str(i) '_time'') = ' num2str(length(time)) ';']);
        end
    end
nc('v3d_time') = length(time);
nc('v2d_time') = length(time);
nc('zeta_time') = length(time);

nc('one') = 1;
%
%  Create variables and attributes
%
nc{'theta_s'} = ncdouble('one') ;
nc{'theta_s'}.long_name = ncchar('S-coordinate surface control parameter');
nc{'theta_s'}.long_name = 'S-coordinate surface control parameter';
nc{'theta_s'}.units = ncchar('nondimensional');
nc{'theta_s'}.units = 'nondimensional';
%
nc{'theta_b'} = ncdouble('one') ;
nc{'theta_b'}.long_name = ncchar('S-coordinate bottom control parameter');
nc{'theta_b'}.long_name = 'S-coordinate bottom control parameter';
nc{'theta_b'}.units = ncchar('nondimensional');
nc{'theta_b'}.units = 'nondimensional';
%
nc{'Tcline'} = ncdouble('one') ;
nc{'Tcline'}.long_name = ncchar('S-coordinate surface/bottom layer width');
nc{'Tcline'}.long_name = 'S-coordinate surface/bottom layer width';
nc{'Tcline'}.units = ncchar('meter');
nc{'Tcline'}.units = 'meter';
%
nc{'hc'} = ncdouble('one') ;
nc{'hc'}.long_name = ncchar('S-coordinate parameter, critical depth');
nc{'hc'}.long_name = 'S-coordinate parameter, critical depth';
nc{'hc'}.units = ncchar('meter');
nc{'hc'}.units = 'meter';
%
nc{'sc_r'} = ncdouble('s_rho') ;
nc{'sc_r'}.long_name = ncchar('S-coordinate at RHO-points');
nc{'sc_r'}.long_name = 'S-coordinate at RHO-points';
nc{'sc_r'}.units = ncchar('nondimensional');
nc{'sc_r'}.units = 'nondimensional';
nc{'sc_r'}.valid_min = -1;
nc{'sc_r'}.valid_max = 0;
%
nc{'Cs_r'} = ncdouble('s_rho') ;
nc{'Cs_r'}.long_name = ncchar('S-coordinate stretching curves at RHO-points');
nc{'Cs_r'}.long_name = 'S-coordinate stretching curves at RHO-points';
nc{'Cs_r'}.units = ncchar('nondimensional');
nc{'Cs_r'}.units = 'nondimensional';
nc{'Cs_r'}.valid_min = -1;
nc{'Cs_r'}.valid_max = 0;
%
nc{'temp_time'} = ncdouble('temp_time') ;
nc{'temp_time'}.long_name = ncchar('time for temperature climatology');
nc{'temp_time'}.long_name = 'time for temperature climatology';
nc{'temp_time'}.units = ncchar('day');
nc{'temp_time'}.units = 'day';
nc{'temp_time'}.cycle_length = cycle;
%
nc{'salt_time'} = ncdouble('salt_time') ;
nc{'salt_time'}.long_name = ncchar('time for temperature climatology');
nc{'salt_time'}.long_name = 'time for temperature climatology';
nc{'salt_time'}.units = ncchar('day');
nc{'salt_time'}.units = 'day';
nc{'salt_time'}.cycle_length = cycle;

for i = 1:NST;
    if i < 10
        eval(['nc{''sand_0' num2str(i) '_time''} = ncdouble(''sand_0' num2str(i) '_time'') ;']);
        eval(['nc{''sand_0' num2str(i) '_time''}.long_name = ncchar(''time for frazil climatology'');']);
        eval(['nc{''sand_0' num2str(i) '_time''}.long_name = ''time for frazil climatology'';']);
        eval(['nc{''sand_0' num2str(i) '_time''}.units = ncchar(''day'');']);
        eval(['nc{''sand_0' num2str(i) '_time''}.units = ''day'';']);
        eval(['nc{''sand_0' num2str(i) '_time''}.cycle_length = cycle;']);
    else
        eval(['nc{''sand_' num2str(i) '_time''} = ncdouble(''sand_' num2str(i) '_time'') ;']);
        eval(['nc{''sand_' num2str(i) '_time''}.long_name = ncchar(''time for frazil climatology'');']);
        eval(['nc{''sand_' num2str(i) '_time''}.long_name = ''time for frazil climatology'';']);
        eval(['nc{''sand_' num2str(i) '_time''}.units = ncchar(''day'');']);
        eval(['nc{''sand_' num2str(i) '_time''}.units = ''day'';']);
        eval(['nc{''sand_' num2str(i) '_time''}.cycle_length = cycle;']);
    end
end
%
nc{'v3d_time'} = ncdouble('v3d_time') ;
nc{'v3d_time'}.long_name = ncchar('time for 3-d momentum');
nc{'v3d_time'}.long_name = 'time for temperature climatology';
nc{'v3d_time'}.units = ncchar('day');
nc{'v3d_time'}.units = 'day';
nc{'v3d_time'}.cycle_length = cycle;
%
nc{'v2d_time'} = ncdouble('v2d_time') ;
nc{'v2d_time'}.long_name = ncchar('time for temperature climatology');
nc{'v2d_time'}.long_name = 'time for temperature climatology';
nc{'v2d_time'}.units = ncchar('day');
nc{'v2d_time'}.units = 'day';
nc{'v2d_time'}.cycle_length = cycle;
%
nc{'zeta_time'} = ncdouble('zeta_time') ;
nc{'zeta_time'}.long_name = ncchar('time for sea surface');
nc{'zeta_time'}.long_name = 'time for temperature climatology';
nc{'zeta_time'}.units = ncchar('day');
nc{'zeta_time'}.units = 'day';
nc{'zeta_time'}.cycle_length = cycle;
%
if obc(1)==1
    %
    %   Southern boundary
    %
    nc{'temp_south'} = ncdouble('temp_time','s_rho','xi_rho') ;
    nc{'temp_south'}.long_name = ncchar('southern boundary potential temperature');
    nc{'temp_south'}.long_name = 'southern boundary potential temperature';
    nc{'temp_south'}.units = ncchar('Celsius');
    nc{'temp_south'}.units = 'Celsius';
    %
    nc{'salt_south'} = ncdouble('salt_time','s_rho','xi_rho') ;
    nc{'salt_south'}.long_name = ncchar('southern boundary salinity');
    nc{'salt_south'}.long_name = 'southern boundary salinity';
    nc{'salt_south'}.units = ncchar('PSU');
    nc{'salt_south'}.units = 'PSU';
    %
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_south''} = ncdouble(''sand_0' num2str(i) '_time'',''s_rho'',''xi_rho'') ;']);
            eval(['nc{''sand_0' num2str(i) '_south''}.long_name = ncchar(''southern boundary frazil'');']);
            eval(['nc{''sand_0' num2str(i) '_south''}.long_name = ''southern boundary frazil'';']);
            eval(['nc{''sand_0' num2str(i) '_south''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_0' num2str(i) '_south''}.units = ''kilogram meter-3'';']);
        else
            eval(['nc{''sand_' num2str(i) '_south''} = ncdouble(''sand_' num2str(i) '_time'',''s_rho'',''xi_rho'') ;']);
            eval(['nc{''sand_' num2str(i) '_south''}.long_name = ncchar(''southern boundary frazil'');']);
            eval(['nc{''sand_' num2str(i) '_south''}.long_name = ''southern boundary frazil'';']);
            eval(['nc{''sand_' num2str(i) '_south''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_' num2str(i) '_south''}.units = ''kilogram meter-3'';']);
        end
    end
    %
    nc{'u_south'} = ncdouble('v3d_time','s_rho','xi_u') ;
    nc{'u_south'}.long_name = ncchar('southern boundary u-momentum component');
    nc{'u_south'}.long_name = 'southern boundary u-momentum component';
    nc{'u_south'}.units = ncchar('meter second-1');
    nc{'u_south'}.units = 'meter second-1';
    %
    nc{'v_south'} = ncdouble('v3d_time','s_rho','xi_rho') ;
    nc{'v_south'}.long_name = ncchar('southern boundary v-momentum component');
    nc{'v_south'}.long_name = 'southern boundary v-momentum component';
    nc{'v_south'}.units = ncchar('meter second-1');
    nc{'v_south'}.units = 'meter second-1';
    %
    nc{'ubar_south'} = ncdouble('v2d_time','xi_u') ;
    nc{'ubar_south'}.long_name = ncchar('southern boundary vertically integrated u-momentum component');
    nc{'ubar_south'}.long_name = 'southern boundary vertically integrated u-momentum component';
    nc{'ubar_south'}.units = ncchar('meter second-1');
    nc{'ubar_south'}.units = 'meter second-1';
    %
    nc{'vbar_south'} = ncdouble('v2d_time','xi_rho') ;
    nc{'vbar_south'}.long_name = ncchar('southern boundary vertically integrated v-momentum component');
    nc{'vbar_south'}.long_name = 'southern boundary vertically integrated v-momentum component';
    nc{'vbar_south'}.units = ncchar('meter second-1');
    nc{'vbar_south'}.units = 'meter second-1';
    %
    nc{'zeta_south'} = ncdouble('zeta_time','xi_rho') ;
    nc{'zeta_south'}.long_name = ncchar('southern boundary sea surface height');
    nc{'zeta_south'}.long_name = 'southern boundary sea surface height';
    nc{'zeta_south'}.units = ncchar('meter');
    nc{'zeta_south'}.units = 'meter';
    %
end
%
if obc(2)==1
    %
    %   Eastern boundary
    %
    nc{'temp_east'} = ncdouble('temp_time','s_rho','eta_rho') ;
    nc{'temp_east'}.long_name = ncchar('eastern boundary potential temperature');
    nc{'temp_east'}.long_name = 'eastern boundary potential temperature';
    nc{'temp_east'}.units = ncchar('Celsius');
    nc{'temp_east'}.units = 'Celsius';
    %
    nc{'salt_east'} = ncdouble('salt_time','s_rho','eta_rho') ;
    nc{'salt_east'}.long_name = ncchar('eastern boundary salinity');
    nc{'salt_east'}.long_name = 'eastern boundary salinity';
    nc{'salt_east'}.units = ncchar('PSU');
    nc{'salt_east'}.units = 'PSU';
    %
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_east''} = ncdouble(''sand_0' num2str(i) '_time'',''s_rho'',''eta_rho'') ;']);
            eval(['nc{''sand_0' num2str(i) '_east''}.long_name = ncchar(''eastern boundary frazil'');']);
            eval(['nc{''sand_0' num2str(i) '_east''}.long_name = ''eastern boundary frazil'';']);
            eval(['nc{''sand_0' num2str(i) '_east''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_0' num2str(i) '_east''}.units = ''kilogram meter-3'';']);
        else
            eval(['nc{''sand_' num2str(i) '_east''} = ncdouble(''sand_' num2str(i) '_time'',''s_rho'',''eta_rho'') ;']);
            eval(['nc{''sand_' num2str(i) '_east''}.long_name = ncchar(''eastern boundary frazil'');']);
            eval(['nc{''sand_' num2str(i) '_east''}.long_name = ''eastern boundary frazil'';']);
            eval(['nc{''sand_' num2str(i) '_east''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_' num2str(i) '_east''}.units = ''kilogram meter-3'';']);
        end
    end
    %
    nc{'u_east'} = ncdouble('v3d_time','s_rho','eta_rho') ;
    nc{'u_east'}.long_name = ncchar('eastern boundary u-momentum component');
    nc{'u_east'}.long_name = 'eastern boundary u-momentum component';
    nc{'u_east'}.units = ncchar('meter second-1');
    nc{'u_east'}.units = 'meter second-1';
    %
    nc{'v_east'} = ncdouble('v3d_time','s_rho','eta_v') ;
    nc{'v_east'}.long_name = ncchar('eastern boundary v-momentum component');
    nc{'v_east'}.long_name = 'eastern boundary v-momentum component';
    nc{'v_east'}.units = ncchar('meter second-1');
    nc{'v_east'}.units = 'meter second-1';
    %
    nc{'ubar_east'} = ncdouble('v2d_time','eta_rho') ;
    nc{'ubar_east'}.long_name = ncchar('eastern boundary vertically integrated u-momentum component');
    nc{'ubar_east'}.long_name = 'eastern boundary vertically integrated u-momentum component';
    nc{'ubar_east'}.units = ncchar('meter second-1');
    nc{'ubar_east'}.units = 'meter second-1';
    %
    nc{'vbar_east'} = ncdouble('v2d_time','eta_v') ;
    nc{'vbar_east'}.long_name = ncchar('eastern boundary vertically integrated v-momentum component');
    nc{'vbar_east'}.long_name = 'eastern boundary vertically integrated v-momentum component';
    nc{'vbar_east'}.units = ncchar('meter second-1');
    nc{'vbar_east'}.units = 'meter second-1';
    %
    nc{'zeta_east'} = ncdouble('zeta_time','eta_rho') ;
    nc{'zeta_east'}.long_name = ncchar('eastern boundary sea surface height');
    nc{'zeta_east'}.long_name = 'eastern boundary sea surface height';
    nc{'zeta_east'}.units = ncchar('meter');
    nc{'zeta_east'}.units = 'meter';
    %
end
%
if obc(3)==1
    %
    %   Northern boundary
    %
    nc{'temp_north'} = ncdouble('temp_time','s_rho','xi_rho') ;
    nc{'temp_north'}.long_name = ncchar('northern boundary potential temperature');
    nc{'temp_north'}.long_name = 'northern boundary potential temperature';
    nc{'temp_north'}.units = ncchar('Celsius');
    nc{'temp_north'}.units = 'Celsius';
    %
    nc{'salt_north'} = ncdouble('salt_time','s_rho','xi_rho') ;
    nc{'salt_north'}.long_name = ncchar('northern boundary salinity');
    nc{'salt_north'}.long_name = 'northern boundary salinity';
    nc{'salt_north'}.units = ncchar('PSU');
    nc{'salt_north'}.units = 'PSU';
    %
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_north''} = ncdouble(''sand_0' num2str(i) '_time'',''s_rho'',''xi_rho'') ;']);
            eval(['nc{''sand_0' num2str(i) '_north''}.long_name = ncchar(''northern boundary frazil'');']);
            eval(['nc{''sand_0' num2str(i) '_north''}.long_name = ''northern boundary frazil'';']);
            eval(['nc{''sand_0' num2str(i) '_north''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_0' num2str(i) '_north''}.units = ''kilogram meter-3'';']);
        else
            eval(['nc{''sand_' num2str(i) '_north''} = ncdouble(''sand_' num2str(i) '_time'',''s_rho'',''xi_rho'') ;']);
            eval(['nc{''sand_' num2str(i) '_north''}.long_name = ncchar(''northern boundary frazil'');']);
            eval(['nc{''sand_' num2str(i) '_north''}.long_name = ''northern boundary frazil'';']);
            eval(['nc{''sand_' num2str(i) '_north''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_' num2str(i) '_north''}.units = ''kilogram meter-3'';']);
        end
    end
    %
    nc{'u_north'} = ncdouble('v3d_time','s_rho','xi_u') ;
    nc{'u_north'}.long_name = ncchar('northern boundary u-momentum component');
    nc{'u_north'}.long_name = 'northern boundary u-momentum component';
    nc{'u_north'}.units = ncchar('meter second-1');
    nc{'u_north'}.units = 'meter second-1';
    %
    nc{'v_north'} = ncdouble('v3d_time','s_rho','xi_rho') ;
    nc{'v_north'}.long_name = ncchar('northern boundary v-momentum component');
    nc{'v_north'}.long_name = 'northern boundary v-momentum component';
    nc{'v_north'}.units = ncchar('meter second-1');
    nc{'v_north'}.units = 'meter second-1';
    %
    nc{'ubar_north'} = ncdouble('v2d_time','xi_u') ;
    nc{'ubar_north'}.long_name = ncchar('northern boundary vertically integrated u-momentum component');
    nc{'ubar_north'}.long_name = 'northern boundary vertically integrated u-momentum component';
    nc{'ubar_north'}.units = ncchar('meter second-1');
    nc{'ubar_north'}.units = 'meter second-1';
    %
    nc{'vbar_north'} = ncdouble('v2d_time','xi_rho') ;
    nc{'vbar_north'}.long_name = ncchar('northern boundary vertically integrated v-momentum component');
    nc{'vbar_north'}.long_name = 'northern boundary vertically integrated v-momentum component';
    nc{'vbar_north'}.units = ncchar('meter second-1');
    nc{'vbar_north'}.units = 'meter second-1';
    %
    nc{'zeta_north'} = ncdouble('zeta_time','xi_rho') ;
    nc{'zeta_north'}.long_name = ncchar('northern boundary sea surface height');
    nc{'zeta_north'}.long_name = 'northern boundary sea surface height';
    nc{'zeta_north'}.units = ncchar('meter');
    nc{'zeta_north'}.units = 'meter';
    %
end
%
if obc(4)==1
    %
    %   Western boundary
    %
    nc{'temp_west'} = ncdouble('temp_time','s_rho','eta_rho') ;
    nc{'temp_west'}.long_name = ncchar('western boundary potential temperature');
    nc{'temp_west'}.long_name = 'western boundary potential temperature';
    nc{'temp_west'}.units = ncchar('Celsius');
    nc{'temp_west'}.units = 'Celsius';
    %
    nc{'salt_west'} = ncdouble('salt_time','s_rho','eta_rho') ;
    nc{'salt_west'}.long_name = ncchar('western boundary salinity');
    nc{'salt_west'}.long_name = 'western boundary salinity';
    nc{'salt_west'}.units = ncchar('PSU');
    nc{'salt_west'}.units = 'PSU';
    %
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_west''} = ncdouble(''sand_0' num2str(i) '_time'',''s_rho'',''eta_rho'') ;']);
            eval(['nc{''sand_0' num2str(i) '_west''}.long_name = ncchar(''western boundary frazil'');']);
            eval(['nc{''sand_0' num2str(i) '_west''}.long_name = ''western boundary frazil'';']);
            eval(['nc{''sand_0' num2str(i) '_west''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_0' num2str(i) '_west''}.units = ''kilogram meter-3'';']);
        else
            eval(['nc{''sand_' num2str(i) '_west''} = ncdouble(''sand_' num2str(i) '_time'',''s_rho'',''eta_rho'') ;']);
            eval(['nc{''sand_' num2str(i) '_west''}.long_name = ncchar(''western boundary frazil'');']);
            eval(['nc{''sand_' num2str(i) '_west''}.long_name = ''western boundary frazil'';']);
            eval(['nc{''sand_' num2str(i) '_west''}.units = ncchar(''kilogram meter-3'');']);
            eval(['nc{''sand_' num2str(i) '_west''}.units = ''kilogram meter-3'';']);
        end
    end
    %
    nc{'u_west'} = ncdouble('v3d_time','s_rho','eta_rho') ;
    nc{'u_west'}.long_name = ncchar('western boundary u-momentum component');
    nc{'u_west'}.long_name = 'western boundary u-momentum component';
    nc{'u_west'}.units = ncchar('meter second-1');
    nc{'u_west'}.units = 'meter second-1';
    %
    nc{'v_west'} = ncdouble('v3d_time','s_rho','eta_v') ;
    nc{'v_west'}.long_name = ncchar('western boundary v-momentum component');
    nc{'v_west'}.long_name = 'western boundary v-momentum component';
    nc{'v_west'}.units = ncchar('meter second-1');
    nc{'v_west'}.units = 'meter second-1';
    %
    nc{'ubar_west'} = ncdouble('v2d_time','eta_rho') ;
    nc{'ubar_west'}.long_name = ncchar('western boundary vertically integrated u-momentum component');
    nc{'ubar_west'}.long_name = 'western boundary vertically integrated u-momentum component';
    nc{'ubar_west'}.units = ncchar('meter second-1');
    nc{'ubar_west'}.units = 'meter second-1';
    %
    nc{'vbar_west'} = ncdouble('v2d_time','eta_v') ;
    nc{'vbar_west'}.long_name = ncchar('western boundary vertically integrated v-momentum component');
    nc{'vbar_west'}.long_name = 'western boundary vertically integrated v-momentum component';
    nc{'vbar_west'}.units = ncchar('meter second-1');
    nc{'vbar_west'}.units = 'meter second-1';
    %
    nc{'zeta_west'} = ncdouble('zeta_time','eta_rho') ;
    nc{'zeta_west'}.long_name = ncchar('western boundary sea surface height');
    nc{'zeta_west'}.long_name = 'western boundary sea surface height';
    nc{'zeta_west'}.units = ncchar('meter');
    nc{'zeta_west'}.units = 'meter';
    %
end
%
%
% Create global attributes
%
nc.title = ncchar(title);
nc.title = title;
nc.date = ncchar(date);
nc.date = date;
nc.clim_file = ncchar(bryname);
nc.clim_file = bryname;
nc.grd_file = ncchar(grdname);
nc.grd_file = grdname;
nc.type = ncchar(type);
nc.type = type;
nc.history = ncchar(history);
nc.history = history;
%
% Leave define mode
%
result = endef(nc);
%
% Compute S coordinates
%
ncload(grdname,'h','zice','lat_rho','lon_rho','mask_rho','mask_zice','angle')%Note: must point to correct grid file!

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
[z,sr,Cs_r]=scoord_zice(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, kgrid, column, index, plt);
kgrid = 1;
[z,sw,Cs_w]=scoord_zice(h', zice', x', y', Vtransform, Vstretching, theta_s, theta_b, ...
                 hc, N, 1, column, index, plt);
             
ang = repmat([angle(2:end,1)',angle(end,:),fliplr(angle(1:end-1,end)')],length(Cs_r),1);

for i  = 1:12;
    [up vp] = do_rotate_vels(squeeze(isfc.u(i,:,:)),squeeze(isfc.v(i,:,:)),ang);
    upd(i,:,:) = up;
    vpd(i,:,:) = vp;
end

% Write variables
%
nc{'theta_s'}(:) =  theta_s;
nc{'theta_b'}(:) =  theta_b;
nc{'Tcline'}(:) =  hc;
nc{'hc'}(:) =  hc;
nc{'sc_r'}(:) =  sr;
nc{'Cs_r'}(:) =  Cs_r;
nc{'temp_time'}(:) =  time;
nc{'salt_time'}(:) =  time;
for i = 1:NST;
    if i < 10
        eval(['nc{''sand_0' num2str(i) '_time''}(:) = time;']);
    else
        eval(['nc{''sand_' num2str(i) '_time''}(:) = time ;']);
    end
end
nc{'v3d_time'}(:) =  time;
nc{'v2d_time'}(:) =  time;
nc{'zeta_time'}(:) =  time;
%
if obc(1)==1
    nc{'u_south'}(:) =  0;
    nc{'v_south'}(:) =  0;
    nc{'ubar_south'}(:) =  0;
    nc{'vbar_south'}(:) =  0;
    nc{'zeta_south'}(:) =  0;
    nc{'temp_south'}(:) =  0;
    nc{'salt_south'}(:) =  0;
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_south''}(:) = ' num2str(0) ';']);
        else
            eval(['nc{''sand_' num2str(i) '_south''}(:) = ' num2str(0) ';']);
        end
    end
end
if obc(2)==1
    nc{'u_east'}(:) = flipdim(flipdim(InterpSurface.u(:,:,M+Lp:end),2),3);
    nc{'v_east'}(:) = flipdim(flipdim(InterpSurface.v(:,:,M+Lp:end-1),2),3);
    nc{'ubar_east'}(:) =  0;
    nc{'vbar_east'}(:) =  0;
    nc{'zeta_east'}(:) =  0;
    nc{'temp_east'}(:) =  flipdim(flipdim(InterpSurface.tmp(:,:,M+Lp:end),2),3);
    nc{'salt_east'}(:) =  flipdim(flipdim(InterpSurface.slt(:,:,M+Lp:end),2),3);
  %  nc{'temp_east'}(:) =  flipdim(InterpSurface.tmp(:,:,1:340),2);
  %  nc{'salt_east'}(:) =  flipdim(InterpSurface.sal(:,:,1:340),2);
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_east''}(:) = ' num2str(0) ';']);
        else
            eval(['nc{''sand_' num2str(i) '_east''}(:) = ' num2str(0) ';']);
        end
    end
end
if obc(3)==1
    nc{'u_north'}(:) =  flipdim(InterpSurface.u(:,:,Mp+1:Mp+Lp-1),2);
    nc{'v_north'}(:) =  flipdim(InterpSurface.v(:,:,Mp+1:Mp+Lp),2);
    nc{'ubar_north'}(:) =  0;
    nc{'vbar_north'}(:) =  0;
    nc{'zeta_north'}(:) =  0;
      nc{'temp_north'}(:) =  flipdim(InterpSurface.tmp(:,:,Mp+1:Mp+Lp),2);
      nc{'salt_north'}(:) =  flipdim(InterpSurface.slt(:,:,Mp+1:Mp+Lp),2);
    %nc{'temp_north'}(:) =  flipdim(InterpSurface.tmp(:,:,340:340+162),2);
    %nc{'salt_north'}(:) =  flipdim(InterpSurface.sal(:,:,340:340+162),2);
    for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_north''}(:) = ' num2str(0) ';']);
        else
            eval(['nc{''sand_' num2str(i) '_north''}(:) = ' num2str(0) ';']);
        end
    end
end
if obc(4)==1
    nc{'u_west'}(:) = flipdim(InterpSurface.u(:,:,1:Mp),2);
    nc{'v_west'}(:) = flipdim(InterpSurface.v(:,:,2:Mp),2);
    nc{'ubar_west'}(:) =  0;
    nc{'vbar_west'}(:) =  0;
    nc{'zeta_west'}(:) =  0;
    nc{'temp_west'}(:) =  flipdim(InterpSurface.tmp(:,:,1:Mp),2);
    nc{'salt_west'}(:) =  flipdim(InterpSurface.slt(:,:,1:Mp),2);
    %size(nc{'temp_west'}(:))
    %size(flipdim(flipdim(InterpSurface.tmp(:,:,340+162:end),2),3))
    %nc{'temp_west'}(:) =  flipdim(flipdim(InterpSurface.tmp(:,:,340+162:end),2),3);
    %nc{'salt_west'}(:) =  flipdim(flipdim(InterpSurface.sal(:,:,340+162:end),2),3);
        for i = 1:NST;
        if i < 10
            eval(['nc{''sand_0' num2str(i) '_west''}(:) = ' num2str(0) ';']);
        else
            eval(['nc{''sand_' num2str(i) '_west''}(:) = ' num2str(0) ';']);
        end
    end
end
close(nc)
return


