disp('loading data for file 1')
filename1 = '/home/ubuntu/IceOceanVolume/tisom015_sgfw/ROMSForcingCreator/lbc/tisom015_1992-2015clima_bry.nc';
grdname1 = '/home/ubuntu/IceOceanVolume/tisom015_sgfw/tisom008_canal_grd.nc';
name1 = 'tisom015_1992-2015clima_bry';
%
Ts1 = ncread(filename1,'theta_s');
Tb1 = ncread(filename1,'theta_b');
Tcline1 = ncread(filename1,'Tcline');
hc1 = ncread(filename1,'hc');
sc_r1 = ncread(filename1,'sc_r');
Cs_r1 = ncread(filename1,'Cs_r');
time1 = ncread(filename1,'temp_time');
temp_north1=ncread(filename1,'temp_north');
temp_west1=ncread(filename1,'temp_west');
temp_east1=ncread(filename1,'temp_east');
salt_north1=ncread(filename1,'salt_north');
salt_west1=ncread(filename1,'salt_west');
salt_east1=ncread(filename1,'salt_east');
u_north1=ncread(filename1,'u_north');
u_west1=ncread(filename1,'u_west');
u_east1=ncread(filename1,'u_east');
v_north1=ncread(filename1,'v_north');
v_west1=ncread(filename1,'v_west');
v_east1=ncread(filename1,'v_east');
%%
disp('loading data for file 2')
filename2 = '/home/ubuntu/IceOceanVolume/tisom016/lbc/tisom016_Cube92_1992-2015IAF_bry.nc';
grdname2 = '/home/ubuntu/IceOceanVolume/tisom016/grid/tisom008_canal_grd.nc';
name2 = 'tisom016_1992-2015IAF_lbc';
%
Ts2 = ncread(filename2,'theta_s');
Tb2 = ncread(filename2,'theta_b');
Tcline2 = ncread(filename2,'Tcline');
hc2 = ncread(filename2,'hc');
sc_r2 = ncread(filename2,'sc_r');
Cs_r2 = ncread(filename2,'Cs_r');
time2 = ncread(filename2,'temp_time');
temp_north2=ncread(filename2,'temp_north');
temp_west2=ncread(filename2,'temp_west');
temp_east2=ncread(filename2,'temp_east');
salt_north2=ncread(filename2,'salt_north');
salt_west2=ncread(filename2,'salt_west');
salt_east2=ncread(filename2,'salt_east');
u_north2=ncread(filename2,'u_north');
u_west2=ncread(filename2,'u_west');
u_east2=ncread(filename2,'u_east');
v_north2=ncread(filename2,'v_north');
v_west2=ncread(filename2,'v_west');
v_east2=ncread(filename2,'v_east');
%%
disp('begin plotting')
set(0,'defaulttextinterpreter','none')
%% three subplots; av temp west,north,east
figure('pos',[100 100 2400 1200])
subplot(3,1,1)
ntitle('temp_north')
plot(time1,squeeze(nanmean(nanmean(temp_north1,1),2)))
hold on,ntitle('temp_north')
plot(time2,squeeze(nanmean(nanmean(temp_north2,1),2)))
legend(name1,name2)
subplot(3,1,2)
ntitle('temp_west')
plot(time1,squeeze(nanmean(nanmean(temp_west1,1),2)))
hold on,ntitle('temp_west')
plot(time2,squeeze(nanmean(nanmean(temp_west2,1),2)))
legend(name1,name2)
subplot(3,1,3)
ntitle('temp_east')
plot(time1,squeeze(nanmean(nanmean(temp_east1,1),2)))
hold on,ntitle('temp_east')
plot(time2,squeeze(nanmean(nanmean(temp_east2,1),2)))
legend(name1,name2)

%% compare climatologies
figure('pos',[100 100 2400 1200])
subplot(3,1,1)
ntitle('temp_north')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_north1,size(temp_north1,1),size(temp_north1,2),12,[]),4)),1),2)))
hold on,ntitle('temp_north')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_north2,size(temp_north2,1),size(temp_north2,2),12,[]),4)),1),2)))
legend(name1,name2)
subplot(3,1,2)
ntitle('temp_west')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_west1,size(temp_west1,1),size(temp_west1,2),12,[]),4)),1),2)))
hold on,ntitle('temp_west')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_west2,size(temp_west2,1),size(temp_west2,2),12,[]),4)),1),2)))
legend(name1,name2)
subplot(3,1,3)
ntitle('temp_east')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_east1,size(temp_east1,1),size(temp_east1,2),12,[]),4)),1),2)))
hold on,ntitle('temp_east')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(temp_east2,size(temp_east2,1),size(temp_east2,2),12,[]),4)),1),2)))
legend(name1,name2)

%% three subplots; av salt west,north,east
figure('pos',[100 100 2400 1200])
subplot(3,1,1)
ntitle('salt_north')
plot(time1,squeeze(nanmean(nanmean(salt_north1,1),2)))
hold on,('salt_north')
plot(time2,squeeze(nanmean(nanmean(salt_north2,1),2)))
legend(name1,name2)
subplot(3,1,2)
ntitle('salt_west')
plot(time1,squeeze(nanmean(nanmean(salt_west1,1),2)))
hold on,('salt_west')
plot(time2,squeeze(nanmean(nanmean(salt_west2,1),2)))
legend(name1,name2)
subplot(3,1,3)
ntitle('salt_east')
plot(time1,squeeze(nanmean(nanmean(salt_east1,1),2)))
hold on,('salt_east')
plot(time2,squeeze(nanmean(nanmean(salt_east2,1),2)))
legend(name1,name2)
%% 
figure('pos',[100 100 2400 1200])
subplot(3,1,1)
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_north1,size(salt_north1,1),size(salt_north1,2),12,[]),4)),1),2)))
hold on,ntitle('salt_north')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_north2,size(salt_north2,1),size(salt_north2,2),12,[]),4)),1),2)))
legend(name1,name2)
subplot(3,1,2)
ntitle('salt_west')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_west1,size(salt_west1,1),size(salt_west1,2),12,[]),4)),1),2)))
hold on,ntitle('salt_west')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_west2,size(salt_west2,1),size(salt_west2,2),12,[]),4)),1),2)))
legend(name1,name2)
subplot(3,1,3)
ntitle('salt_east')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_east1,size(salt_east1,1),size(salt_east1,2),12,[]),4)),1),2)))
hold on,ntitle('salt_east')
plot([365/12/2:365/12:365-365/12/2],squeeze(nanmean(nanmean(squeeze(nanmean(reshape(salt_east2,size(salt_east2,1),size(salt_east2,2),12,[]),4)),1),2)))
legend(name1,name2)
%%%




set(0,'defaulttextinterpreter','latex')

