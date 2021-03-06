% Make domain
Xloc = [xmin,xmax]; %lon
Yloc = [ymin,ymax]; %lat

% Initilise time vector (important)
DateIndexVector=datestr(datenum(MinYear,1,2):3:datenum(MaxYear+1,1,1),'yyyymmdd');
TimeLength=length(DateIndexVector);
disp('Errrr will have to check this')
% Initialise data matrix for ECCO2 data.
theta = nan(TimeLength,50,xmax-xmin+1,ymax-ymin+1);
    
% Load cube92-3day data within subset

for TimeInd = 1:TimeLength
    
    
theta(TimeInd,:,:,:) = permute(ncread(['/ds/projects/iomp/obs/ECCO2/cube92_real/THETA.nc/THETA.1440x720x50.' num2str(DateIndexVector(TimeInd,:)) '.nc'],'THETA',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[3 1 2]);

if ~rem(TimeInd,round(TimeLength/10))
disp([num2str(TimeInd/TimeLength*100) '% at ' datestr(now)])
end
end

save(['cube92_3day_clima_theta_',RunName,'.mat'],'theta','-v7.3')
clear theta
