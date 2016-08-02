% make climatology of tamura forcing
shfluxClima_tmp= reshape(shfluxClima_tmp(1:365*(MaxYear-MinYear+1),:,:), [365 (MaxYear-MinYear+1) size(shfluxClima_tmp,2) size(shfluxClima_tmp,3)]);
ssfluxClima_tmp= reshape(ssfluxClima_tmp(1:365*(MaxYear-MinYear+1),:,:), [365 (MaxYear-MinYear+1) size(ssfluxClima_tmp,2) size(ssfluxClima_tmp,3)]);
ssfluxClima=squeeze(nanmean(ssfluxClima_tmp,2)); %has size [365 x y]
shfluxClima=squeeze(nanmean(shfluxClima_tmp,2));


% find mean and anomaly of okane forcing 

% add climatology of tamura forcing to anomaly of okane forcing
