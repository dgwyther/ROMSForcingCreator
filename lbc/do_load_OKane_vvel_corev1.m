% Make domain

% load Terry O'Kane-COREv1 data

    
   
vvel = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_v_100E-130E_60S-68S.nc'],'v'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv1_iaf_vvel_',RunName,'.mat'],'vvel','-v7.3')
clear theta
