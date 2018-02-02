%% Make domain
Xloc = [xmin xmax]; %lon
Yloc = [ymin ymax]; %lat



%% Initialise data matrix for ECCO2 data.
theta = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92
year_ind = MinYear:MaxYear;

for TimeInd = 1:length(year_ind);
theta((TimeInd-1)*12+1:TimeInd*12,:,:,:) = permute(ncread(['/ds/projects/iomp/obs/ECCO2/cube92_real/THETA_MonthlyFrom3Day.nc/THETA_' num2str(year_ind(TimeInd)) '_MonMeans.nc'],'THETA',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[4 3 1 2]);
end

save(['cube92_iaf_theta_',RunName,'.mat'],'theta','-v7.3')

