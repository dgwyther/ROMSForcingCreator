% Make domain

% load Terry O'Kane-COREv1 data

    
   
vvel = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_v_95E-135E_50S-70S.nc'],'v'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv1_iaf_vvel_',RunName,'.mat'],'vvel','-v7.3')
clear theta
