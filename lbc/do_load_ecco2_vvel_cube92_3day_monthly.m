%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
vvel = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92
year_ind = MinYear:MaxYear;

for TimeInd = 1:length(year_ind);
vvel((TimeInd-1)*12+1:TimeInd*12,:,:,:) = permute(ncread(['/ds/projects/iomp/obs/ECCO2/cube92_real/VVEL_MonthlyFrom3Day.nc/VVEL_' num2str(year_ind(TimeInd)) '_MonMeans.nc'],'VVEL',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[4 3 1 2]);
end

save(['cube92_iaf_vvel_',RunName,'.mat'],'vvel','-v7.3')

