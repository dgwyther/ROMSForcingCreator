% Make domain

% load Terry O'Kane-COREv2 data
salt = permute(ncread(['/ds/projects/iomp/obs/TOK/CORE2_data/ocean_month_salt_95E-135E_50S-70S.nc'],'salt'),[4 3 1 2]); %format time depth x y


save(['OKane-COREv2_iaf_salt_',RunName,'.mat'],'salt','-v7.3')
clear salt
