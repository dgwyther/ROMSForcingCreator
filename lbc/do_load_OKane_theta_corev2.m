% Make domain

% load Terry O'Kane-COREv1 data

    
   
theta = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE2_data/ocean_month_temp_95E-135E_50S-70S.nc'],'temp'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv2_iaf_theta_',RunName,'.mat'],'theta','-v7.3')
clear theta
