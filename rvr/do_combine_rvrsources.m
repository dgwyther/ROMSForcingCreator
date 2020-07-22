function do_combine_rvrsources(grdName,fileOne,fileTwo,outName)
%% load rvr sources and locs for FILE 1
load(fileOne,'ROMS_X_pos','ROMS_E_pos','ROMS_direction','ROMS_flux_multiplier','ROMS_river_fluxes','ROMS_river_lonlat_filtered')
ROMS_X_pos_s = ROMS_X_pos; ROMS_E_pos_s = ROMS_E_pos; ROMS_direction_s = ROMS_direction; ROMS_flux_multiplier_s = ROMS_flux_multiplier; ROMS_river_fluxes_s = ROMS_river_fluxes; ROMS_river_lonlat_filtered_s=ROMS_river_lonlat_filtered;
%% load rvr sources and locs for FILE 2
load(fileTwo,'ROMS_X_pos','ROMS_E_pos','ROMS_direction','ROMS_flux_multiplier','ROMS_river_fluxes','ROMS_river_lonlat_filtered')
ROMS_X_pos_c = ROMS_X_pos; ROMS_E_pos_c = ROMS_E_pos; ROMS_direction_c = ROMS_direction; ROMS_flux_multiplier_c = ROMS_flux_multiplier; ROMS_river_fluxes_c = ROMS_river_fluxes; ROMS_river_lonlat_filtered_c=ROMS_river_lonlat_filtered;
%% combine sources and locs
ROMS_X_pos = [ROMS_X_pos_s,ROMS_X_pos];
ROMS_E_pos = [ROMS_E_pos_s,ROMS_E_pos_c];
ROMS_direction = [ROMS_direction_s,ROMS_direction_c];
ROMS_flux_multiplier = [ROMS_flux_multiplier_s,ROMS_flux_multiplier_c];
ROMS_river_fluxes = [ROMS_river_fluxes_s,ROMS_river_fluxes_c];
ROMS_river_lonlat = [ROMS_river_lonlat_filtered_s;ROMS_river_lonlat_filtered_c];
%% remove doubles 
Data_mat = [ROMS_X_pos',ROMS_E_pos',ROMS_direction',ROMS_flux_multiplier',ROMS_river_fluxes'];
[C,ia,ic] = unique(Data_mat(:,1:3),'rows','stable');
Data_mat_filtered = nan(max(ic),5); %5 cols to be X,E,dir,fluxmult,flux
ROMS_river_lonlat_filtered = nan(max(ic),2);
for jjj=1:max(ic)
Data_mat_filtered(jjj,1)=unique(Data_mat(ic==jjj,1));
Data_mat_filtered(jjj,2)=unique(Data_mat(ic==jjj,2));
Data_mat_filtered(jjj,3)=unique(Data_mat(ic==jjj,3));
Data_mat_filtered(jjj,4)=unique(Data_mat(ic==jjj,4));
Data_mat_filtered(jjj,5)=sum(Data_mat(ic==jjj,5));
ROMS_river_lonlat_filtered(jjj,1)=unique(ROMS_river_lonlat(ic==jjj,1));
ROMS_river_lonlat_filtered(jjj,2)=unique(ROMS_river_lonlat(ic==jjj,2));
end
%% write new file
ROMS_X_pos = Data_mat_filtered(:,1)';
ROMS_E_pos = Data_mat_filtered(:,2)';
ROMS_direction = Data_mat_filtered(:,3)';
ROMS_flux_multiplier = Data_mat_filtered(:,4)';
ROMS_river_fluxes = Data_mat_filtered(:,5)';
save(outName,'ROMS_X_pos','ROMS_E_pos','ROMS_direction','ROMS_flux_multiplier','ROMS_river_fluxes','ROMS_river_lonlat_filtered');

nansum(ROMS_river_fluxes_s')
nansum(ROMS_river_fluxes_c')
nansum(ROMS_river_fluxes_s')+nansum(ROMS_river_fluxes_c')
nansum(Data_mat_filtered(:,5))

lon_rho = ncread(grdName,'lon_rho');
lat_rho = ncread(grdName,'lat_rho');
h = ncread(grdName,'h');
figure('position',[250 750 2400 1200]);
subaxis(3,1,1,'margin',0.035);
flat_pcolor(lon_rho,lat_rho,h);
axis([112.9 121 -67.6 -66])
hold on,hhh3=bubbleplot(ROMS_river_lonlat_filtered_s(:,1),ROMS_river_lonlat_filtered_s(:,2),[],ROMS_river_fluxes_s');
set(hhh3,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
hhh4=bubbleplot([118 118 118 118 118],[-67.1,-67.2,-67.3 -67.4 -67.5],[],[max(ROMS_river_fluxes_s') max(ROMS_river_fluxes_s')/10 max(ROMS_river_fluxes_s')/20 max(ROMS_river_fluxes_s')/100 1e-3]);
set(hhh4,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
text([118 118 118 118 118]+.1,[-67.1,-67.2,-67.3 -67.4 -67.5],{num2str(max(ROMS_river_fluxes_s')) num2str(max(ROMS_river_fluxes_s')/10) num2str(max(ROMS_river_fluxes_s')/20) num2str(max(ROMS_river_fluxes_s')/100) '1e-3'},'color',[.8 .8 .8])
ntitle(['fluxes from first file; total = ',num2str(nansum(ROMS_river_fluxes_s')),' m3/s'])
subaxis(3,1,2);
flat_pcolor(lon_rho,lat_rho,h);
axis([112.9 121 -67.6 -66])
hold on,hhh3=bubbleplot(ROMS_river_lonlat_filtered_c(:,1),ROMS_river_lonlat_filtered_c(:,2),[],ROMS_river_fluxes_c');
set(hhh3,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
hhh4=bubbleplot([118 118 118 118 118],[-67.1,-67.2,-67.3 -67.4 -67.5],[],[max(ROMS_river_fluxes_c') max(ROMS_river_fluxes_c')/10 max(ROMS_river_fluxes_c')/20 max(ROMS_river_fluxes_c')/100 1e-3]);
set(hhh4,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
text([118 118 118 118 118]+.1,[-67.1,-67.2,-67.3 -67.4 -67.5],{num2str(max(ROMS_river_fluxes_c')) num2str(max(ROMS_river_fluxes_c')/10) num2str(max(ROMS_river_fluxes_c')/20) num2str(max(ROMS_river_fluxes_c')/100) '1e-3'},'color',[.8 .8 .8])
ntitle(['fluxes from second file; total = ',num2str(nansum(ROMS_river_fluxes_c')),' m3/s'])
subaxis(3,1,3);
flat_pcolor(lon_rho,lat_rho,h);
axis([112.9 121 -67.6 -66])
hold on,hhh3=bubbleplot(ROMS_river_lonlat_filtered(:,1),ROMS_river_lonlat_filtered(:,2),[],Data_mat_filtered(:,5));
set(hhh3,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
hhh4=bubbleplot([118 118 118 118 118],[-67.1,-67.2,-67.3 -67.4 -67.5],[],[max(Data_mat_filtered(:,5)) max(Data_mat_filtered(:,5))/10 max(Data_mat_filtered(:,5))/20 max(Data_mat_filtered(:,5))/100 1e-3]);
set(hhh4,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
text([118 118 118 118 118]+.1,[-67.1,-67.2,-67.3 -67.4 -67.5],{num2str(max(Data_mat_filtered(:,5))) num2str(max(Data_mat_filtered(:,5))/10) num2str(max(Data_mat_filtered(:,5))/20) num2str(max(Data_mat_filtered(:,5))/100) '1e-3'},'color',[.8 .8 .8])
ntitle(['combined fluxes; total = ',num2str(nansum(Data_mat_filtered(:,5))),' m3/s'])

