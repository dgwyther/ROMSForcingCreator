% Make domain

% load Terry O'Kane-COREv1 data

    
   
uvel = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE1_data/ocean_month_1800-1899_u_100E-130E_60S-68S.nc'],'u'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv1_iaf_uvel_',RunName,'.mat'],'uvel','-v7.3')
clear theta
