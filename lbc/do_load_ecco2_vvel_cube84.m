% Load model grid:
%ncload(grdname,'lon_rho','lat_rho','lon_u','lat_u','lon_v','lat_v')

%lon_rho = ncread(grdname,'lon_rho')';
%lat_rho = ncread(grdname,'lat_rho')';
%lon_u = ncread(grdname,'lon_u')';
%lat_u = ncread(grdname,'lat_u')';
%lon_v = ncread(grdname,'lon_v')';
%lat_v = ncread(grdname,'lat_v')';

% Find nearest indexes for locations:
%ncload /ds/projects/iomp/obs/ECCO2/cube84/VVEL/VVEL.1440x720x50.month001.nc LATITUDE_T LONGITUDE_T DEPTH_T
LATITUDE_T=ncread('/ds/projects/iomp/obs/ECCO2/cube78/VVEL/VVEL.1440x720x50.month001.nc','LATITUDE_T')';
LONGITUDE_T=ncread('/ds/projects/iomp/obs/ECCO2/cube78/VVEL/VVEL.1440x720x50.month001.nc','LONGITUDE_T')';
DEPTH_T=ncread('/ds/projects/iomp/obs/ECCO2/cube78/VVEL/VVEL.1440x720x50.month001.nc','DEPTH_T')';

[lons lats] = meshgrid(LONGITUDE_T,LATITUDE_T);
%% Make domain

% % For MISOM:
% xmax = 525;%650;
% xmin = 410;%450;
% ymax = 125;%140; % <----- is loaded by make_lbc.m
% ymin = 80;%40;

Xloc = [xmin:xmax];
Yloc = [ymin:ymax];



%% Find relevant ECCO2 data.
vvel = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
for TimeInd = 1:12*(MaxYear-MinYear+1)
    
    tosave=TimeInd; tosave
    save TimeInd_vvel.dat  tosave -ascii

       Vin = zeros(50,1440,720); 

        if TimeInd  < 10
        eval(['Vin = shiftdim(reshape(fread(fopen(''/ds/projects/iomp/obs/ECCO2/cube84/VVEL/VVEL.1440x720x50.month00' num2str(TimeInd) ''',''r'',''b''),1440*720*50,''float32=>double''),1440,720,50),50);']);
        elseif TimeInd < 100
            eval(['Vin = shiftdim(reshape(fread(fopen(''/ds/projects/iomp/obs/ECCO2/cube84/VVEL/VVEL.1440x720x50.month0' num2str(TimeInd) ''',''r'',''b''),1440*720*50,''float32=>double''),1440,720,50),50);']);
        else
            eval(['Vin = shiftdim(reshape(fread(fopen(''/ds/projects/iomp/obs/ECCO2/cube84/VVEL/VVEL.1440x720x50.month' num2str(TimeInd) ''',''r'',''b''),1440*720*50,''float32=>double''),1440,720,50),50);']);
        end

            
        vvel(TimeInd,:,:,:) = squeeze(Vin(:,Xloc,Yloc));   
    
end

save(['cube84_iaf_vvel_',RunName,'.mat'],'vvel','-v7.3')
