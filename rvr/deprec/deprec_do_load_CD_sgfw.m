%load ROMS data
h=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','h');
lon_rho=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lon_rho');
lat_rho=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lat_rho');
lon_u=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lon_u');
lat_u=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lat_u');
lon_v=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lon_v');
lat_v=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','lat_v');
mask_rho=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','mask_rho');
mask_u=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','mask_u');
mask_v=ncread('/mnt/IceOceanVolume/tisom015/grid/tisom008_canal_grd.nc','mask_v');

%make mask which is shifted E/W (u) and N/S (v)
mask_u1 = (mask_rho(1:end-1,:)+mask_rho(2:end,:))/2;
mask_v1 = (mask_rho(:,1:end-1)+mask_rho(:,2:end))/2;
% mask_u-mask_u1 gives the NEXT row with a change from 0-1.
diff_mask_u=logical(mask_u-mask_u1);
lon_u_bnds=lon_u.*double(diff_mask_u);
lon_u_bnds(lon_u_bnds==0)=NaN;
lat_u_bnds=lat_u.*double(diff_mask_u);
lat_u_bnds(lat_u_bnds==0)=NaN;
lon_u_bnds(isnan(lon_u_bnds))=[];
lat_u_bnds(isnan(lat_u_bnds))=[];
%now same for changes in mask_v
diff_mask_v=logical(mask_v-mask_v1);
lon_v_bnds=lon_v.*double(diff_mask_v);
lon_v_bnds(lon_v_bnds==0)=NaN;
lat_v_bnds=lat_v.*double(diff_mask_v);
lat_v_bnds(lat_v_bnds==0)=NaN;
lon_v_bnds(isnan(lon_v_bnds))=[];
lat_v_bnds(isnan(lat_v_bnds))=[];
bnd_coords_cat = [lon_u_bnds,lon_v_bnds;lat_u_bnds,lat_v_bnds]'; %make into 2 cols, many rows
bnd_coords_uORvindex = [ones(length(lon_u_bnds),1);-1*ones(length(lon_v_bnds),1)];

%remove data within some defined boundaries:
%'ear canal' island
polygon_eci = [[117.5 117 116.3 116.3];[-66.95 -67.1 -67.225 -66.95]]';
in.eci = inpolygon(bnd_coords_cat(:,1),bnd_coords_cat(:,2),polygon_eci(:,1),polygon_eci(:,2));
in.eciu = inpolygon(lon_u_bnds,lat_u_bnds,polygon_eci(:,1),polygon_eci(:,2));
in.eciv = inpolygon(lon_v_bnds,lat_v_bnds,polygon_eci(:,1),polygon_eci(:,2));

bnd_coords_cat(in.eci,:)=[];
bnd_coords_uORvindex(in.eci)=[];
lon_u_bnds(in.eciu)=[];
lat_u_bnds(in.eciu)=[];
lon_v_bnds(in.eciv)=[];
lat_v_bnds(in.eciv)=[];


if 0
%testing
b = bwboundaries(mask_rho,'noholes');
x = b{1}(:,1);
y = b{1}(:,2);
X = reshape(bsxfun(@plus,x,[0 -0.5 0.5]),[],1);
Y = reshape(bsxfun(@plus,y,[0 0.5 -0.5]),[],1);
k = boundary(X,Y,1);
figure,imagesc(mask_rho);
hold on
plot(Y(k),X(k),'g','LineWidth',3,'Marker','o',...
    'MarkerFaceColor','r','MarkerEdgeColor','none')
end

%work out if each u/v face is on left(1)/top(2)/right(3)/bottom(4) of water...
where=pdist2([lon_u_bnds;lat_u_bnds]',[lon_u(:),lat_u(:)]);
for ii=1:length(lon_u_bnds);
[~, where_idx(ii)]= min(where(ii,:));
[grid_u_locs(ii,1),grid_u_locs(ii,2)]=ind2sub(size(lon_u),where_idx(ii));
end
clear where_idx where
where=pdist2([lon_v_bnds;lat_v_bnds]',[lon_v(:),lat_v(:)]);
for ii=1:length(lon_v_bnds);
[~, where_idx(ii)]= min(where(ii,:));
[grid_v_locs(ii,1),grid_v_locs(ii,2)]=ind2sub(size(lon_v),where_idx(ii));
end
clear u_dir_vec
for ii=1:size(grid_u_locs,1)
if mask_rho(grid_u_locs(ii,1)+1,grid_u_locs(ii,2)) ==1 
u_dir_vec(ii)=1;
elseif mask_rho(grid_u_locs(ii,1)-1,grid_u_locs(ii,2)) ==1
u_dir_vec(ii)=3;
elseif mask_rho(grid_u_locs(ii,1)+1,grid_u_locs(ii,2)) ==0
u_dir_vec(ii)=3;
end
end
clear v_dir_vec
for ii=1:size(grid_v_locs,1)
if mask_rho(grid_v_locs(ii,1),grid_v_locs(ii,2)+1) ==1
v_dir_vec(ii)=2;
elseif mask_rho(grid_v_locs(ii,1),grid_v_locs(ii,2)-1) ==1
v_dir_vec(ii)=4;
elseif mask_rho(grid_v_locs(ii,1),grid_v_locs(ii,2)+1) ==0
v_dir_vec(ii)=4;
end
end

uANDv_dir_vec=[u_dir_vec,v_dir_vec];

if 0
%1 rightwards flow, 2=upwards flow, 3=leftwards flow, 4=downwards flow
for ii=1:size(grid_u_locs,1) %plot check
if u_dir_vec(ii)==1
plot(bnd_coords_cat(ii,1),bnd_coords_cat(ii,2),'bs','markersize',10)
end
if u_dir_vec(ii)==3
plot(bnd_coords_cat(ii,1),bnd_coords_cat(ii,2),'mx','markersize',10)
end
end
for ii=1:size(grid_v_locs,1) %plot check
if v_dir_vec(ii)==2
plot(lon_v_bnds(ii),lat_v_bnds(ii),'wo','markersize',10)
end
if v_dir_vec(ii)==4
plot(lon_v_bnds(ii),lat_v_bnds(ii),'wd','markersize',10)
end
end
end

if 0
% Sort to try and draw line
[~,I]=sort([lon_u_bnds,lon_v_bnds]);
bnd_coords_inc1 = bnd_coords_cat(I,1);
bnd_coords_inc2 = bnd_coords_cat(I,2);
bnd_coords = [bnd_coords_inc1,bnd_coords_inc2];
bnd_coords_dist = pdist2(bnd_coords,bnd_coords);
bnd_coords_sort = nan(1,size(bnd_coords,1));
bnd_coords_sort(1)=1;
for ii=2:size(bnd_coords,1)
 bnd_coords_dist(:,bnd_coords_sort(ii-1)) = Inf; %don't go backwards?
 [~,closest_idx] = min(bnd_coords_dist(bnd_coords_sort(ii-1),:));
 bnd_coords_sort(ii)=closest_idx;
end
bnd_coords_final(:,1)=bnd_coords(bnd_coords_sort,1);
bnd_coords_final(:,2)=bnd_coords(bnd_coords_sort,2);
uORv_sorted= bnd_coords_uORvindex(bnd_coords_sort); %1 for u and -1 for v
end

% load river data
channel_or_sheet = 'channel'
FlowType = 'T8e-5';
%%
if strcmp(channel_or_sheet,'channel')
Input_data=csvread('Totten_outflow_channel.csv',0,1); %ignore first col header
elseif strcmp(channel_or_sheet,'sheet')
Input_data=csvread('Totten_outflow_sheet.csv',0,1); %ignore first col header
else 
disp('choose either channel or sheet')
end

river_X=Input_data(1,:); %x loc
river_Y=Input_data(2,:); %y loc
% choose type:
%T1. Tel 1e-4
%T5. Tel 5e-4 
%T8. Tel 8e-5
%B1. Bed 1e-4
%B5. Bed 5e-4
%B8. Bed 8e-5
if strcmp(FlowType,'T1e-4')
river_flux = Input_data(3,:); disp('flow type chosen is Tel 1e-4')
elseif strcmp(FlowType,'T5e-4')
river_flux = Input_data(4,:); disp('flow type chosen is Tel 5e-4')
elseif strcmp(FlowType,'T8e-5')
river_flux = Input_data(5,:); disp('flow type chosen is Tel 8e-5')
elseif strcmp(FlowType,'B1e-4')
river_flux = Input_data(6,:); disp('flow type chosen is Bed 1e-4')
elseif strcmp(FlowType,'B5e-4')
river_flux = Input_data(7,:); disp('flow type chosen is Bed 5e-4')
elseif strcmp(FlowType,'B8e-5')
river_flux = Input_data(8,:); disp('flow type chosen is Bed 8e-5')
else
disp('choose correct flow type')
end

[river_lat,river_lon] = ps2ll(river_X,river_Y);

GlaDS_lonlat=[river_lon',river_lat'];
ROMS_lonlat = bnd_coords_cat;%[[lon_u(:);lon_v(:)],[lat_u(:);lat_v(:)]];


dists = pdist2(GlaDS_lonlat,ROMS_lonlat);
closest_ROMS_lonlat = nan(length(GlaDS_lonlat(:,2)),2);

for i = 1:size(GlaDS_lonlat,1)
    [~, S_idx(i)]= min(dists(i,:));
    closest_ROMS_lonlat(i,:) = ROMS_lonlat(S_idx(i),:);
    uORv_ROMS_lonlat(i) = bnd_coords_uORvindex(S_idx(i)); %1 for u, -1 for v
    dir_ROMS_lonlat(i) = uANDv_dir_vec(S_idx(i)); %directions of ROMS points with flow through them.
end



if 0 %old fig
figure('position',[250 750 1000 1000])
flat_pcolor(lon_rho,lat_rho,h);
axis([113 119 -67.5 -67])
hold on,plot(river_lon,river_lat,'x','markersize',3);

%plot lines between original location and new location
plot(closest_ROMS_lonlat(:,1),closest_ROMS_lonlat(:,2),'ro','markersize',5);
for ii=1:length(river_X)
hold on
plot([GlaDS_lonlat(ii,1),closest_ROMS_lonlat(ii,1)],[GlaDS_lonlat(ii,2),closest_ROMS_lonlat(ii,2)],'r-')
end
end

%So I have the following data:
% closest_ROMS_lonlat is the closest location for each river...
% GlaDS_lonlat is the original river location...
% uORv_ROMS_lonlat is index location for u (+1) or v(-1)
% u_dir_vec is index vector for flow to the east (1) or west (3)
% v_dir_vec is index vector for flow to the north (2) or south (4)
% uANDv_dir_vec is directions of u, then v points.

SUM_SIMILAR_LOCATIONS = 0;
if SUM_SIMILAR_LOCATIONS
disp('This doesn''t calculate fluxes correctly - bugged')
for ii=1:size(closest_ROMS_lonlat,1) %loop through srcs
[locations_logical_lonlat(:,ii)]=ismember(closest_ROMS_lonlat,closest_ROMS_lonlat(ii,:),'rows');
temp1=closest_ROMS_lonlat(:,1); temp2=closest_ROMS_lonlat(:,2); ROMS_river_lonlat(ii,:) = unique([temp1(locations_logical_lonlat(:,ii)),temp2(locations_logical_lonlat(:,ii))],'rows');
ROMS_river_fluxes(ii)=sum(river_flux(locations_logical_lonlat(:,ii)));
ROMS_river_uORv(ii,:) = unique(uORv_ROMS_lonlat(locations_logical_lonlat(:,ii)));
if ROMS_river_uORv(ii,:)==1
ROMS_river_dir(ii) = unique(dir_ROMS_lonlat(locations_logical_lonlat(:,ii)));
elseif ROMS_river_uORv(ii,:)==-1
ROMS_river_dir(ii) = unique(dir_ROMS_lonlat(locations_logical_lonlat(:,ii)));
end
end
else
ROMS_river_lonlat=closest_ROMS_lonlat;
ROMS_river_fluxes=river_flux;
ROMS_river_uORv=uORv_ROMS_lonlat;
ROMS_river_dir=dir_ROMS_lonlat;

end


figure('position',[250 750 1000 1000])
flat_pcolor(lon_rho,lat_rho,h);
axis([113 119 -67.5 -67])
hold on,plot(river_lon,river_lat,'x','markersize',3);

for ii=1:size(closest_ROMS_lonlat,1)
if ROMS_river_dir(ii)==1
plot(ROMS_river_lonlat(ii,1),ROMS_river_lonlat(ii,2),'r>')
elseif ROMS_river_dir(ii)==2
plot(ROMS_river_lonlat(ii,1),ROMS_river_lonlat(ii,2),'r^')
elseif ROMS_river_dir(ii)==3
plot(ROMS_river_lonlat(ii,1),ROMS_river_lonlat(ii,2),'r<')
elseif ROMS_river_dir(ii)==4
plot(ROMS_river_lonlat(ii,1),ROMS_river_lonlat(ii,2),'rv')
end



end



%% outputting
% ROMS wants X and E position of each point
% 0 for u face, 1 for v face
% flux multiplier *1 for north/east, *-1 for west/south
% flux for each point


clear ROMS_X_pos ROMS_E_pos
index_vec_v = 1:size(lon_v,1)*size(lon_v,2); index_vec_v = reshape(index_vec_v,size(lon_v,1),size(lon_v,2));
index_vec_u = 1:size(lon_u,1)*size(lon_u,2); index_vec_u = reshape(index_vec_u,size(lon_u,1),size(lon_u,2));
for ii=1:size(closest_ROMS_lonlat,1)
if ROMS_river_uORv(ii)==1 %u point
testing = griddata(lon_u,lat_u,index_vec_u,closest_ROMS_lonlat(ii,1),closest_ROMS_lonlat(ii,2),'nearest');
[ROMS_X_pos(ii),ROMS_E_pos(ii)]=ind2sub(size(lon_u),testing);
elseif ROMS_river_uORv(ii)==-1 %v point
testing = griddata(lon_v,lat_v,index_vec_v,closest_ROMS_lonlat(ii,1),closest_ROMS_lonlat(ii,2),'nearest');
[ROMS_X_pos(ii),ROMS_E_pos(ii)]=ind2sub(size(lon_v),testing);
end
end

ROMS_direction=ROMS_river_uORv; ROMS_direction(ROMS_direction==-1)=0; ROMS_direction = double(~logical(ROMS_direction)); %0 for u, 1 for v.

ROMS_flux_multiplier = ROMS_river_dir; ROMS_flux_multiplier(ROMS_flux_multiplier==3 | ROMS_flux_multiplier==4)=-1; ROMS_flux_multiplier(ROMS_flux_multiplier==2)=1;%1 for north/east, -1 for west/south

%only plot u's

figure('position',[250 750 1000 1000])
flat_pcolor(lon_rho,lat_rho,h);
axis([113 119 -67.5 -67])
hold on,plot(river_lon,river_lat,'x','markersize',3);
for ii=1:length(uORv_ROMS_lonlat)
if uORv_ROMS_lonlat(ii)==1

%plot lines between original location and new location
plot(closest_ROMS_lonlat(ii,1),closest_ROMS_lonlat(ii,2),'ro','markersize',5);
hold on
plot([GlaDS_lonlat(ii,1),closest_ROMS_lonlat(ii,1)],[GlaDS_lonlat(ii,2),closest_ROMS_lonlat(ii,2)],'k-');
end
end

for ii=1:length(uORv_ROMS_lonlat)
if uORv_ROMS_lonlat(ii)==-1

%plot lines between original location and new location
plot(closest_ROMS_lonlat(ii,1),closest_ROMS_lonlat(ii,2),'ro','markersize',5);
hold on
plot([GlaDS_lonlat(ii,1),closest_ROMS_lonlat(ii,1)],[GlaDS_lonlat(ii,2),closest_ROMS_lonlat(ii,2)],'k-');
end
end


figure('position',[250 750 1000 1000])
flat_pcolor(lon_rho,lat_rho,h),
axis([113 119 -67.5 -67])
hold on,plot(river_lon,river_lat,'x','markersize',3)
for ii=1:length(uORv_ROMS_lonlat)
if uORv_ROMS_lonlat(ii)==1

%plot lines between original location and new location
plot(lon_u(ROMS_X_pos(ii),ROMS_E_pos(ii)),lat_u(ROMS_X_pos(ii),ROMS_E_pos(ii)),'ro','markersize',5);
hold on
plot([GlaDS_lonlat(ii,1),lon_u(ROMS_X_pos(ii),ROMS_E_pos(ii))],[GlaDS_lonlat(ii,2),lat_u(ROMS_X_pos(ii),ROMS_E_pos(ii))],'k-')
end
end
for ii=1:length(uORv_ROMS_lonlat)
if uORv_ROMS_lonlat(ii)==-1

%plot lines between original location and new location
plot(lon_v(ROMS_X_pos(ii),ROMS_E_pos(ii)),lat_v(ROMS_X_pos(ii),ROMS_E_pos(ii)),'ro','markersize',5);
hold on
plot([GlaDS_lonlat(ii,1),lon_v(ROMS_X_pos(ii),ROMS_E_pos(ii))],[GlaDS_lonlat(ii,2),lat_v(ROMS_X_pos(ii),ROMS_E_pos(ii))],'k-')
end
end


%convert between MATLAB and ROMS numbering
% see:
% https://www.myroms.org/forum/viewtopic.php?f=18&t=3809
%Explanation:
%1. Array indices start from 1 in MATLAB.
%2. Array indices start from 0 in Fortran (ROMS).
%3. Array dimensions are flipped between MATLAB (eta, xi) and ROMS (xi, eta) notations.
%4. rho point indices in MATLAB are one more of what they are in ROMS. i.e. [code]ROMS(xi, eta) -> MATLAB(eta+1, xi+1)
%5. For u and v points this is slightly difference since these arrays have one member short in one of their dimensions. ROMS solves hydrodynamic equations over the interior points (i=1:Lm and j=1:Nm) for a full grid size of (i=0:L, j=0:M). i.e. Lm=L-1, Mm=M-1 and full range rho(0:L, 0:M), whereas, u(1:Lm, 0:M), v(0:L,1:M). See here for more details.
%e.g.
%Lm=158 means 1:158 interior points in ROMS which corresponds to 0:159 full range points in ROMS and 1:160 full range points in MATLAB.
%Additional Notes:
%1. Assume that the cell with rho(xi= 80,eta= 54) in the above example is a wet cell neighboring a land masked cell with rho(xi= 79,eta= 54) on the east of it. If defining a river / point source coming in from the east the location will be as river_Xposition=80, river_Eposition=54 in ROMS river input file.
%2. You can test your river / point source indices by running map_rivers.m code, which can be found on this post.

for ii=1:length(ROMS_direction)
if ROMS_direction(ii)==0 %u
ROMS_E_pos(ii) = ROMS_E_pos(ii)-1;
elseif ROMS_direction(ii)==1 %v
ROMS_X_pos(ii) = ROMS_X_pos(ii)-1;
end
end
save('Totten_river_sources.mat','ROMS_X_pos','ROMS_E_pos','ROMS_direction','ROMS_flux_multiplier','ROMS_river_fluxes');

