%% make grid

grdname = 'test_grd.nc';
ModelName = 'Amery Ice Shelf Ocean Model v10';

grid_type = 'geographic_grid'

%%%% start making grid

if strcmp(grid_type,'geographic_grid')
res = 1/30; 
xc = [104.5,130]; 
yc = -[-68,-60];

Lats = [];
MinLat = yc(2);
ii = 1;
while MinLat < yc(1)

    Lats(ii) = MinLat;
    MinLat = MinLat+(cosd(MinLat)*resolution);
    ii = ii+1;
end
Lats = fliplr(-Lats);

Lons = [xc(1):resolution:xc(2)];

[lon lat] = meshgrid(Lons,Lats);
[x y]=polar_stereo_deluxe(lat,lon,0,0,1,-71);
x=x*1000; y=y*1000; %km->m

elseif grid_type,'polarstereo_grid')
res = 1; %km
xmin=2210;
xmax=2467;
ymin=398;
ymax=603;
[x,y] = ndgrid(xmin:res/2:xmax+res/2),ymin:res/2:ymax+res/2);
x_rho = nan(size(x_rho,1)+2,size(x_rho,2)+2); x_rho(2:end-1,2:end-1)=x; BREAK HERE
y_rho = y[1::2,1::2]

%calculate lat,lon points
lat,lon = inverse_polar_stereo(x,y,0,0,-71.0)
%calculate curvilinear coordinate distances at rho points
dx = haversine(lon[1::2,0:-2:2],lat[1::2,0:-2:2],lon[1::2,2::2],lat[1::2,2::2])
dy = haversine(lon[0:-2:2,1::2],lat[0:-2:2,1::2],lon[2::2,1::2],lat[2::2,1::2])

#calculate curvilinear coordinate metrices
pm = 1./dx;
pn = 1./dy;
 
dndx = np.empty_like(pm)
dmde = np.empty_like(pn)

dndx[:,1:-1] = 0.5*(pn[:,2:] - pn[:,:-2])
dmde[1:-1,:] = 0.5*(pm[2:,:] - pm[:-2,:])

dndx[:,0]  = 2*dndx[:,1]  - dndx[:,2]
dndx[:,-1] = 2*dndx[:,-2] - dndx[:,-3]
dmde[0,:]  = 2*dmde[1,:]  - dmde[2,:]
dmde[-1,:] = 2*dmde[-2,:] - dmde[-3,:]

%subset lat and lon at rho, psi, u and v points
lon_rho = lon[1::2,1::2]
lat_rho = lat[1::2,1::2]

lon_psi = lon[2:-1:2,2:-1:2]
lat_psi = lat[2:-1:2,2:-1:2]

lon_u = lon[1::2,2:-1:2]
lat_u = lat[1::2,2:-1:2]

lon_v = lon[2:-1:2,1::2]
end