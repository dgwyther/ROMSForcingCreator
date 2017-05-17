%% Make domain
Xloc = [xmin:xmax]; %lon
Yloc = [ymin:ymax]; %lat



%% Initialise data matrix for ECCO2 data.
vvel = nan(12*(MaxYear-MinYear+1),50,xmax-xmin+1,ymax-ymin+1);
    
%Cycle through correct naming convention for cube92
year_ind = MinYear:MaxYear;
mon_ind = 1:12;
loop_ind=1;
clear monstr
for yy=MinYear:MaxYear; %number of years
 for mm=1:12 %number of months
  if mm<10
  monstr(loop_ind)=str2num([num2str(yy),'0',num2str(mm)]);
  elseif mm>=10
  monstr(loop_ind)=str2num([num2str(yy),num2str(mm)]);
  end
  loop_ind=loop_ind+1;
  end
end


for TimeInd = 1:12*(MaxYear-MinYear+1)

vvel(TimeInd,:,:,:) = permute(ncread(['/ds/projects/iomp/obs/ECCO2/cube92_real/VVEL_monthly.nc/VVEL.1440x720x50.' num2str(monstr(TimeInd)) '.nc'],'VVEL',[Xloc(1) Yloc(1) 1 1],[Xloc(2)-Xloc(1)+1 Yloc(2)-Yloc(1)+1 Inf Inf]),[3 1 2]);

end

save([WorkingDir,'cube92_iaf_vvel_',ModelName,RunNo,'.mat'],'vvel','-v7.3')
