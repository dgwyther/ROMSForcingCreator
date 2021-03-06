% Make domain

% load Terry O'Kane-COREv2 data
uvel = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE2_data/ocean_month_u_95E-135E_50S-70S.nc'],'u'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv2_iaf_uvel_',RunName,'.mat'],'uvel','-v7.3')
clear uvel
