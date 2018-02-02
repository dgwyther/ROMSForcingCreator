filename1_shf = '/home/ubuntu/IceOceanVolume/tisom009/sbc/tisom009_dailyclima_sbc.nc';
filename1_ssf = '/home/ubuntu/IceOceanVolume/tisom009/sbc/tisom009_dailyclima_sbc.nc';
filename1_uwn = '/home/ubuntu/IceOceanVolume/tisom009/sbc/tisom009_dailyclima_sbc.nc';
filename1_vwn = '/home/ubuntu/IceOceanVolume/tisom009/sbc/tisom009_dailyclima_sbc.nc';
filename1_grd = '/home/ubuntu/IceOceanVolume/tisom009/grid/tisom008_canal_grd.nc';
name1 ='tisom009_dailyclima_sbc.nc';

filename2_shf = '/home/ubuntu/IceOceanVolume/tisom016/sbc/tisom016_1992-2015IAF_shf_sbc.nc';
filename2_swf = '/home/ubuntu/IceOceanVolume/tisom016/sbc/tisom016_1992-2015IAF_swf_sbc.nc';
filename2_uwn = '/home/ubuntu/IceOceanVolume/tisom016/sbc/tisom016_1992-2015IAF_su_sbc.nc';
filename2_vwn = '/home/ubuntu/IceOceanVolume/tisom016/sbc/tisom016_1992-2015IAF_sv_sbc.nc';
filename2_grd = '/home/ubuntu/IceOceanVolume/tisom016/grid/tisom008_canal_grd.nc';
name2 = 'tisom016_1992-2015IAF_*_sbc.nc';

%load from filename 1
lon_rho1 = ncread(filename1_grd,'lon_rho');
lat_rho1 = ncread(filename1_grd,'lat_rho');
shflux1 = ncread(filename1_shf,'shflux');
ssflux1 = ncread(filename1_ssf,'ssflux');
uwnd1 = ncread(filename1_su,'uwnd');
vwnd1 = ncread(filename1_sv,'vwnd');
zice1 = ncread(filename1_grd,'zice'); 
h1 = ncread(filename1_grd,'h');
mask_zice1 = ncread(filename1_grd,'mask_zice'); mask_zice1(mask_zice1==0)=NaN;
mask_rho1 = ncread(filename1_grd,'mask_rho'); mask_rho1(mask_rho1==0)=NaN;

%load from filename 2
lon_rho2 = ncread(filename2_grd,'lon_rho');
lat_rho2 = ncread(filename2_grd,'lat_rho');
shflux2 = ncread(filename2_shf,'shflux');
ssflux2 = ncread(filename2_ssf,'ssflux');
uwnd2 = ncread(filename2_su,'uwnd');
vwnd2 = ncread(filename2_sv,'vwnd');
zice2 = ncread(filename2_grd,'zice');
h2 = ncread(filename2_grd,'h');
mask_zice2 = ncread(filename2_grd,'mask_zice'); mask_zice2(mask_zice2==0)=NaN;
mask_rho2 = ncread(filename2_grd,'mask_rho'); mask_rho2(mask_rho2==0)=NaN;


%% Do comparison of different surface forcing conditions
%1. Compare mean of shflux and ssflux
figure
subplot(2,3,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(shflux1,3)).*mask_rho1); colorbar
subplot(2,3,2)
flat_pcolor(lon_rho2,lat_rho2,squeeze(nanmean(shflux2,3)).*mask_rho1); colorbar
subplot(2,3,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(shflux2,3)-nanmean(shflux1,3)).*mask_rho1); colorbar, cmocean('balance','pivot',0)
subplot(2,3,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(ssflux1,3)).*mask_rho2); colorbar
subplot(2,3,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(nanmean(ssflux2,3)).*mask_rho2); colorbar
subplot(2,3,6)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(ssflux2,3)-nanmean(ssflux1,3)).*mask_rho2); colorbar, cmocean('balance','pivot',0)

%2. Compare mean of u-wind and v-wind
figure
subplot(2,3,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(uwnd1,3)).*mask_rho1); colorbar
subplot(2,3,2)
flat_pcolor(lon_rho2,lat_rho2,squeeze(nanmean(uwnd2,3)).*mask_rho1); colorbar
subplot(2,3,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(uwnd2,3)-nanmean(uwnd1,3)).*mask_rho1); colorbar, cmocean('balance','pivot',0)
subplot(2,3,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(vwnd1,3)).*mask_rho2); colorbar
subplot(2,3,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(nanmean(vwnd2,3)).*mask_rho2); colorbar
subplot(2,3,6)
flat_pcolor(lon_rho1,lat_rho1,squeeze(nanmean(vwnd2,3)-nanmean(vwnd1,3)).*mask_rho2); colorbar, cmocean('balance','pivot',0)

%3. Compare climatology (months 1, 4, 7, 10) of shflux and ssflux
shflux1_c = squeeze(nanmean(reshape(shflux1,size(shflux1,1),size(shflux1,2),12,[]),4));
shflux2_c = squeeze(nanmean(reshape(shflux2,size(shflux2,1),size(shflux2,2),12,[]),4));
ssflux1_c = squeeze(nanmean(reshape(ssflux1,size(ssflux1,1),size(ssflux1,2),12,[]),4));
ssflux2_c = squeeze(nanmean(reshape(ssflux2,size(ssflux2,1),size(ssflux2,2),12,[]),4));

figure
subplot(3,4,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux1_c(:,:,1)).*mask_rho1); colorbar
ntitle([name1,' shflux for Jan'])
subplot(3,4,2)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux1_c(:,:,4)).*mask_rho1); colorbar
ntitle([name1,' shflux for Apr'])
subplot(3,4,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux1_c(:,:,7)).*mask_rho1); colorbar
ntitle([name1,' shflux for Jul'])
subplot(3,4,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux1_c(:,:,10)).*mask_rho1); colorbar
ntitle([name1,' shflux for Oct'])
subplot(3,4,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(shflux2_c(:,:,1)).*mask_rho2); colorbar
ntitle([name1,' shflux for Jan'])
subplot(3,4,6)
flat_pcolor(lon_rho2,lat_rho2,squeeze(shflux2_c(:,:,4)).*mask_rho2); colorbar
ntitle([name1,' shflux for Apr'])
subplot(3,4,7)
flat_pcolor(lon_rho2,lat_rho2,squeeze(shflux2_c(:,:,7)).*mask_rho2); colorbar
ntitle([name1,' shflux for Jul'])
subplot(3,4,8)
flat_pcolor(lon_rho2,lat_rho2,squeeze(shflux2_c(:,:,10)).*mask_rho2); colorbar
ntitle([name1,' shflux for Oct'])
subplot(3,4,9)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux2_c(:,:,1)-shflux1_c(:,:,1)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' shflux for Jan'])
subplot(3,4,10)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux2_c(:,:,4)-shflux1_c(:,:,4)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' shflux for Apr'])
subplot(3,4,11)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux2_c(:,:,7)-shflux1_c(:,:,7)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' shflux for Jul'])
subplot(3,4,12)
flat_pcolor(lon_rho1,lat_rho1,squeeze(shflux2_c(:,:,10)-shflux1_c(:,:,10)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' shflux for Oct'])

figure
subplot(3,4,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux1_c(:,:,1)).*mask_rho1); colorbar
ntitle([name1,' ssflux for Jan'])
subplot(3,4,2)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux1_c(:,:,4)).*mask_rho1); colorbar
ntitle([name1,' ssflux for Apr'])
subplot(3,4,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux1_c(:,:,7)).*mask_rho1); colorbar
ntitle([name1,' ssflux for Jul'])
subplot(3,4,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux1_c(:,:,10)).*mask_rho1); colorbar
ntitle([name1,' ssflux for Oct'])
subplot(3,4,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(ssflux2_c(:,:,1)).*mask_rho2); colorbar
ntitle([name1,' ssflux for Jan'])
subplot(3,4,6)
flat_pcolor(lon_rho2,lat_rho2,squeeze(ssflux2_c(:,:,4)).*mask_rho2); colorbar
ntitle([name1,' ssflux for Apr'])
subplot(3,4,7)
flat_pcolor(lon_rho2,lat_rho2,squeeze(ssflux2_c(:,:,7)).*mask_rho2); colorbar
ntitle([name1,' ssflux for Jul'])
subplot(3,4,8)
flat_pcolor(lon_rho2,lat_rho2,squeeze(ssflux2_c(:,:,10)).*mask_rho2); colorbar
ntitle([name1,' ssflux for Oct'])
subplot(3,4,9)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux2_c(:,:,1)-ssflux1_c(:,:,1)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' ssflux for Jan'])
subplot(3,4,10)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux2_c(:,:,4)-ssflux1_c(:,:,4)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' ssflux for Apr'])
subplot(3,4,11)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux2_c(:,:,7)-ssflux1_c(:,:,7)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' ssflux for Jul'])
subplot(3,4,12)
flat_pcolor(lon_rho1,lat_rho1,squeeze(ssflux2_c(:,:,10)-ssflux1_c(:,:,10)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' ssflux for Oct'])


%4. Compare climatology (months 1, 4, 7, 10) of u-wind and v-wind
uwnd1_c = squeeze(nanmean(reshape(uwnd1,size(uwnd1,1),size(uwnd1,2),12,[]),4));
uwnd2_c = squeeze(nanmean(reshape(uwnd2,size(uwnd2,1),size(uwnd2,2),12,[]),4));
vwnd1_c = squeeze(nanmean(reshape(vwnd1,size(vwnd1,1),size(vwnd1,2),12,[]),4));
vwnd2_c = squeeze(nanmean(reshape(vwnd2,size(vwnd2,1),size(vwnd2,2),12,[]),4));


figure
subplot(3,4,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd1_c(:,:,1)).*mask_rho1); colorbar
ntitle([name1,' uwnd for Jan'])
subplot(3,4,2)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd1_c(:,:,4)).*mask_rho1); colorbar
ntitle([name1,' uwnd for Apr'])
subplot(3,4,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd1_c(:,:,7)).*mask_rho1); colorbar
ntitle([name1,' uwnd for Jul'])
subplot(3,4,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd1_c(:,:,10)).*mask_rho1); colorbar
ntitle([name1,' uwnd for Oct'])
subplot(3,4,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(uwnd2_c(:,:,1)).*mask_rho2); colorbar
ntitle([name1,' uwnd for Jan'])
subplot(3,4,6)
flat_pcolor(lon_rho2,lat_rho2,squeeze(uwnd2_c(:,:,4)).*mask_rho2); colorbar
ntitle([name1,' uwnd for Apr'])
subplot(3,4,7)
flat_pcolor(lon_rho2,lat_rho2,squeeze(uwnd2_c(:,:,7)).*mask_rho2); colorbar
ntitle([name1,' uwnd for Jul'])
subplot(3,4,8)
flat_pcolor(lon_rho2,lat_rho2,squeeze(uwnd2_c(:,:,10)).*mask_rho2); colorbar
ntitle([name1,' uwnd for Oct'])
subplot(3,4,9)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd2_c(:,:,1)-uwnd1_c(:,:,1)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' uwnd for Jan'])
subplot(3,4,10)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd2_c(:,:,4)-uwnd1_c(:,:,4)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' uwnd for Apr'])
subplot(3,4,11)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd2_c(:,:,7)-uwnd1_c(:,:,7)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' uwnd for Jul'])
subplot(3,4,12)
flat_pcolor(lon_rho1,lat_rho1,squeeze(uwnd2_c(:,:,10)-uwnd1_c(:,:,10)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' uwnd for Oct'])



figure
subplot(3,4,1)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd1_c(:,:,1)).*mask_rho1); colorbar
ntitle([name1,' vwnd for Jan'])
subplot(3,4,2)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd1_c(:,:,4)).*mask_rho1); colorbar
ntitle([name1,' vwnd for Apr'])
subplot(3,4,3)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd1_c(:,:,7)).*mask_rho1); colorbar
ntitle([name1,' vwnd for Jul'])
subplot(3,4,4)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd1_c(:,:,10)).*mask_rho1); colorbar
ntitle([name1,' vwnd for Oct'])
subplot(3,4,5)
flat_pcolor(lon_rho2,lat_rho2,squeeze(vwnd2_c(:,:,1)).*mask_rho2); colorbar
ntitle([name1,' vwnd for Jan'])
subplot(3,4,6)
flat_pcolor(lon_rho2,lat_rho2,squeeze(vwnd2_c(:,:,4)).*mask_rho2); colorbar
ntitle([name1,' vwnd for Apr'])
subplot(3,4,7)
flat_pcolor(lon_rho2,lat_rho2,squeeze(vwnd2_c(:,:,7)).*mask_rho2); colorbar
ntitle([name1,' vwnd for Jul'])
subplot(3,4,8)
flat_pcolor(lon_rho2,lat_rho2,squeeze(vwnd2_c(:,:,10)).*mask_rho2); colorbar
ntitle([name1,' vwnd for Oct'])
subplot(3,4,9)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd2_c(:,:,1)-vwnd1_c(:,:,1)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' vwnd for Jan'])
subplot(3,4,10)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd2_c(:,:,4)-vwnd1_c(:,:,4)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' vwnd for Apr'])
subplot(3,4,11)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd2_c(:,:,7)-vwnd1_c(:,:,7)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' vwnd for Jul'])
subplot(3,4,12)
flat_pcolor(lon_rho1,lat_rho1,squeeze(vwnd2_c(:,:,10)-vwnd1_c(:,:,10)).*mask_rho2); colorbar, cmocean('balance','pivot',0)
ntitle([name2,'-',name1,' vwnd for Oct'])



