function do_combine_rvrsources(fileOne,fileTwo,outName)
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

figure('position',[250 750 2400 1200])
flat_pcolor(lon_rho,lat_rho,h),
ntitle(['Volume fluxes (m3/s) on ROMS grid flow'])
axis([112.9 121 -67.6 -66])
hold on,hhh3=bubbleplot(ROMS_river_lonlat_filtered(:,1),ROMS_river_lonlat_filtered(:,2),[],Data_mat_filtered(:,5));
set(hhh3,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
hhh4=bubbleplot([118 118 118 118 118],[-67.35,-67.375,-67.4 -67.425 -67.45],[],[max(Data_mat_filtered(:,5)) 1 .5 .1 1e-2]);
set(hhh4,'MarkerEdgeColor',[0 0 0],...
              'MarkerFaceColor',[1 .7 0],...
              'LineWidth',1.5)
text([118 118 118 118 118]+.1,[-67.35,-67.375,-67.4 -67.425 -67.45],{num2str(max(Data_mat_filtered(:,5))) '1' '5e-1' '1e-1' '1e-2'},'color',[.8 .8 .8])


