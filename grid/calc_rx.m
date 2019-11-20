function [rx0,rx1]=calc_rx(roms_grid,Cs_w,N)
% [rx0,rx1]=calc_rx(roms_grid,Cs_w,N)

h = ncread(roms_grid,'h');
zice = ncread(roms_grid,'zice');
mask_zice= zice; mask_zice(zice>=0)=0; mask_zice(zice<0)=1;
mzN = mask_zice; mzN(mzN==0)=NaN;

h2=h ; %.*mzN;
z2=zice ; %.*mzN;

z_w = bsxfun(@plus,bsxfun(@times,h2+z2,shiftdim(Cs_w,-2)),z2);
for i=2:size(h,1)
 for j=2:size(h,2)
  for k=2:N
  rx0_x(i,j) = abs((z_w(i,j,1)-z_w(i-1,j,1))/(z_w(i,j,1)+z_w(i-1,j,1)));
  rx0_y(i,j) = abs((z_w(i,j,1)-z_w(i,j-1,1))/(z_w(i,j,1)+z_w(i,j-1,1)));
  rx1_x(i,j,k) = abs((z_w(i,j,k)-z_w(i,j-1,k)+z_w(i,j,k-1)-z_w(i,j-1,k-1))/(z_w(i,j,k)+z_w(i,j-1,k)-z_w(i,j,k-1)-z_w(i,j-1,k-1)));
  rx1_y(i,j,k) = abs((z_w(i,j,k)-z_w(i,j-1,k)+z_w(i,j,k-1)-z_w(i,j-1,k-1))/(z_w(i,j,k)+z_w(i,j-1,k)-z_w(i,j,k-1)-z_w(i,j-1,k-1)));
  end
 end
end
rx1 = max(max(rx1_x,rx1_y),[],3);
rx0 = max(rx0_x,rx0_y);


