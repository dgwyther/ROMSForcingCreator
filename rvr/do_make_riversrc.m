grdname='/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc'
RunName='tisom015_sgfw'
frcname='tisom015_river_src_CombinedShtChn.nc';
channel_or_sheet=''
InertTracers = [1]; % vector of 1's equal to number of tracers <---- this defines the number of tracers you want.
ZERO_TRANSPORT=0  %if you want to run with tracers and turn on river sources after a restart. This adds rivers and dye with 0 influence, but fills the netcdf with placeholders.
NO_NEGATIVE_TRANSPORT=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title =([RunName,' river source file']); % Make title from model name variable

N = 31; %number of s-levels
Lpinfo=ncinfo(grdname,'lon_rho'); Lp = Lpinfo.Size(1);
Mpinfo=ncinfo(grdname,'lat_rho'); Mp = Mpinfo.Size(2);
L=Lp-1;
M=Mp-1;


disp('load input file')
if strcmp(channel_or_sheet,'channel') | strcmp(channel_or_sheet,'sheet')
disp(['using ',['/mnt/IceOceanVolume/tisom015_sgfw/rvr/Totten_river_sources_',channel_or_sheet,'_.mat']'])
load('/mnt/IceOceanVolume/tisom015_sgfw/rvr/Totten_river_sources_',channel_or_sheet,'_.mat')
else
disp(['using ','/mnt/IceOceanVolume/tisom015_sgfw/rvr/Totten_river_sources.mat'])
load('/mnt/IceOceanVolume/tisom015_sgfw/rvr/Totten_river_sources.mat')
end

USE_DEBUG=0
if USE_DEBUG
disp('in debug mode')
ROMS_X_pos(11:end)=[];
ROMS_E_pos(11:end)=[];
ROMS_direction(11:end)=[];
ROMS_flux_multiplier(11:end)=[];
ROMS_river_fluxes(11:end)=[];
ROMS_river_fluxes(1:10)=1000;
end

if ZERO_TRANSPORT
ROMS_river_fluxes(1:end)=0;
end

if NO_NEGATIVE_TRANSPORT
ROMS_river_fluxes(ROMS_river_fluxes<0)=0;
end

river_time = [1]; % for fixed time / perpetual calendar %[1:365]/86400;
river_data = ROMS_river_fluxes;
river_data_multiplier = ROMS_flux_multiplier;
river_transport = repmat(river_data.*river_data_multiplier,[length(river_time),1]);
river_Xposition = ROMS_X_pos;
river_Eposition = ROMS_E_pos;
river_direction = ROMS_direction;
river_flag = 3*ones(size(river_data));
river_index = 1:length(river_data);
river_no = length(river_index);

disp('define tracer properties')
% set river tracer props
Salinity = 0;
zice = ncread(grdname,'zice'); lat_rho=ncread(grdname,'lat_rho');
for ss=1:length(river_Xposition)
%if river_direction(ss) == 0 %u-face
location_i = river_Xposition(ss)+1; %applies similarly for either u or v-face
location_j = river_Eposition(ss)+1;
%elseif river_direction(ss) == 1 %v-face
%location_i = river_Xposition(ss)
%location_j = river_Eposition(ss)
%end
Temperature(ss)=gsw_t_freezing(Salinity,gsw_p_from_z(zice(location_i,location_j),lat_rho(location_i,location_j)));
end
river_temp = permute(repmat(Temperature',[1,length(river_time),N]),[3 2 1]);
%river_temp = repmat(0,[length(river_time),N,length(river_data)]);
river_salt = repmat(Salinity,[length(river_time),N,length(river_data)]);
if exist('InertTracers')
IniConc = 1;
for tt=1:length(InertTracers)
if tt<10, dyename = ['0',num2str(tt)]; else, dyename = num2str(tt); end
eval(['river_dye_',dyename,' = repmat(IniConc,[length(river_time),N,length(river_data)]);'])
end
end


%set vertical distribution of water
vshape_dist = [1:N]; %linear dist scaling with highest flux at top.
%vshape_dist = ones(1,N)/N; %even flow out of each layer
vshape_default = vshape_dist/sum(vshape_dist);
river_Vshape = repmat(vshape_default',[1,river_no]);

time_big = river_time; %days
cycle = time_big;

time_perpetual=1; %just use a 'none' calendar.


disp(' ')
disp([' Creating the file : ',frcname])
disp(' ')


id = netcdf.create(frcname, 'clobber');
%
%  Create dimensions
%
    disp('Create dimension')

%xi_u_dim = netcdf.defDim(id, 'xi_u', L);
%eta_u_dim = netcdf.defDim(id, 'eta_u', Mp);
%xi_v_dim = netcdf.defDim(id, 'xi_v', Lp);
%eta_v_dim = netcdf.defDim(id, 'eta_v', M);
xi_rho_dim = netcdf.defDim(id, 'xi_rho', Lp);
eta_rho_dim = netcdf.defDim(id, 'eta_rho', Mp);
s_rho_dim = netcdf.defDim(id, 's_rho', N);
river_dim = netcdf.defDim(id, 'river', river_no);
river_time_dim = netcdf.defDim(id, 'river_time', length(time_big));
%
%  Create variables and attributes
%
    disp('Create variable attribute')
river_id = netcdf.defVar(id, 'river', 'double', river_dim);
netcdf.putAtt(id, river_id, 'long_name', 'river index');

river_Xposition_id = netcdf.defVar(id, 'river_Xposition', 'double', river_dim);
netcdf.putAtt(id, river_Xposition_id, 'long_name', 'River xi index');

river_Eposition_id = netcdf.defVar(id, 'river_Eposition', 'double', river_dim);
netcdf.putAtt(id, river_Eposition_id, 'long_name', 'River eta index');

river_direction_id = netcdf.defVar(id, 'river_direction', 'double', river_dim);
netcdf.putAtt(id, river_direction_id, 'long_name', 'River direction');

river_flag_id = netcdf.defVar(id, 'river_flag', 'double', river_dim);
netcdf.putAtt(id, river_flag_id, 'long_name', 'River flag');
netcdf.putAtt(id, river_flag_id, 'option_0','all tracers are off');
netcdf.putAtt(id, river_flag_id, 'option_1','only temperature is on');
netcdf.putAtt(id, river_flag_id, 'option_2','only salinity is on');
netcdf.putAtt(id, river_flag_id, 'option_3','both temperature and salinity are on');

river_Vshape_id = netcdf.defVar(id, 'river_Vshape', 'double', [river_dim s_rho_dim]);
netcdf.putAtt(id, river_Vshape_id, 'long_name', 'River vertical shape');

river_time_id = netcdf.defVar(id, 'river_time', 'double', river_time_dim);
netcdf.putAtt(id, river_time_id, 'long_name', 'River time');
netcdf.putAtt(id, river_time_id, 'units', 'day');
if time_perpetual
netcdf.putAtt(id, river_time_id, 'calendar','none');
elseif time_cycle
netcdf.putAtt(id, river_time_id, 'cycle_length', cycle);
end
river_transport_id = netcdf.defVar(id, 'river_transport', 'double', [river_dim river_time_dim]);
netcdf.putAtt(id, river_transport_id, 'long_name', 'River transport');
netcdf.putAtt(id, river_transport_id, 'units', 'm3/s');

river_temp_id = netcdf.defVar(id, 'river_temp', 'double', [river_dim s_rho_dim river_time_dim]);
netcdf.putAtt(id, river_temp_id, 'long_name', 'River temperature');
netcdf.putAtt(id, river_temp_id, 'units', 'deg C');

river_salt_id = netcdf.defVar(id, 'river_salt', 'double', [river_dim s_rho_dim river_time_dim]);
netcdf.putAtt(id, river_salt_id, 'long_name', 'River salinity');

if exist('InertTracers')
for tt=1:length(InertTracers)
if tt<10, dyename = ['dye_','0',num2str(tt)]; else, dyename = ['dye_',num2str(tt)]; end
eval(['river_',dyename,'_id = netcdf.defVar(id, ''river_',dyename,''', ''double'', [river_dim s_rho_dim river_time_dim]);'])
eval(['netcdf.putAtt(id, river_',dyename,'_id, ''long_name'', ''River ',dyename,'concentration'');'])
end
end


%
% Create global attributes
%
    disp('Create global attribute')
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'title', title);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'date', date);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'grd_file', grdname);
netcdf.putAtt(id, netcdf.getConstant('NC_GLOBAL'), 'type', 'ROMS forcing file');
netcdf.endDef(id);


    disp('Write time variable')

netcdf.putVar(id, river_time_id, river_time);

    disp('write river data variables')

netcdf.putVar(id, river_id, river_index);
netcdf.putVar(id, river_Xposition_id, river_Xposition);
netcdf.putVar(id, river_Eposition_id, river_Eposition);
netcdf.putVar(id, river_direction_id, river_direction);
netcdf.putVar(id, river_flag_id, river_flag);
netcdf.putVar(id, river_Vshape_id, permute(river_Vshape,[2 1]));
netcdf.putVar(id, river_transport_id, permute(river_transport, [2 1]));
netcdf.putVar(id, river_temp_id, permute(river_temp, [3 2 1]));
netcdf.putVar(id, river_salt_id, permute(river_salt, [3 2 1]));


if exist('InertTracers')
for tt=1:length(InertTracers)
if tt<10, dyename = ['dye_','0',num2str(tt)]; else, dyename = ['dye_',num2str(tt)]; end
eval(['netcdf.putVar(id, river_',dyename,'_id, permute(river_',dyename,', [3 2 1]));'])

end
end


netcdf.close(id);

    disp('DONE')


figure('position',[600 600   1000        1000])
map_rivers(grdname,frcname)
axis square,axis([113.5 120 -67.5 -66.6])
%axis([117 118 -67.2 -66.9])
