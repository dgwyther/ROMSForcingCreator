%% do_Read_fastice_ADF
%%
%% read in fastice data from Alex Fraser
%% Fraser, A.D., R.A. Massom, K.J. Michael, B.K. Galton-Fenzi, and J.L. Lieser, 2012: East Antarctic Landfast Sea Ice Distribution and Variability, 2000–08. J. Climate, 25, 1137–1156
%% https://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-10-05032.1
%% 
%% Set the cutoff for how persistent over teh 2000-2008 climatology the ice should be to be considered fast ice (PersistenceCutoff)
%% Set the concentration over the 2000-2008 climatology the ice should be to be considered fast ice (ConcentrationCutoff)
%%
%% DEG 09/03/2018

ClimaCutoff = 0.5 % To be fast ice, ice must be present for this % during each 20-day bracket in the climatology. (%/100)
InterannualCutoff = 0.9 % To be included in the mean distribution must be present at least this % of the total time. (%/100)

%% Read in coordinates:
%
lats = reshape(fread(fopen('/ds/projects/iomp/obs/fastice/lats.img'),4300*425,'float32=>double'),4300,425);
lons= reshape(fread(fopen('/ds/projects/iomp/obs/fastice/lons.img'),4300*425,'float32=>double'),4300,425);
%

MonthNums = ['001-020';...
'021-040';...
'041-060';...
'061-080';...
'081-100';...
'101-120';...
'121-140';...
'141-160';...
'161-180';...
'181-200';...
'201-220';...
'221-240';...
'241-260';...
'261-280';...
'281-300';...
'301-320';...
'321-340';...
'341-365'];

i = 1;
fastice = [];
months = [];
for YearInd = 2000:2008
    i,
    if YearInd == 2000
        MonthInds = [4:18];
    else
        MonthInds = [1:18];
    end

    months = [months, MonthInds];

   for MonthInd = MonthInds

eval(['a = imread(''/ds/projects/iomp/obs/fastice/pngs/' num2str(YearInd) '_' MonthNums(MonthInd,:) 'coloured.png'');']);
%%% RGB indexing:
% ocean/icebergs (light blue) 0   255   255
% method 1 fast-ice (red)   255     0     0
% method 2 fast-ice (yellow) 255   255     0
% method 3 fast-ice (blue)    0     0   255
% land/ice shelves (White) 255   255   255

R = squeeze(a(:,:,1));
G = squeeze(a(:,:,2));
B = squeeze(a(:,:,3));

MaskIce = nan(size(a,1),size(a,2));
MaskIce((R(:) == 0 & G(:) == 255 & B(:) == 255))  = 0;

MaskIce( (R(:) == 255 & G(:) == 0 & B(:) == 0) | ...
    (R(:) == 255 & G(:) == 255 & B(:) == 0) | ...
    (R(:) == 0 & G(:) == 0 & B(:) == 255))  = 1;
fastice(i,:,:) = MaskIce;

     i = i+1;
    end
end

%% First calculate the % persistence for each 20-day bracket. e.g. 70% for 01-20 at a certain location = certain location has fast ice for 70% of the time over days 01-20
fastice_clima = [];
for i = 1:18;
ii = find(months == i);
fastice_clima(i,:,:) = squeeze(nanmean(fastice(ii,:,:)));
end
%% Set a cutoff (e.g. 0.5) which defines fastice if it exists for 50% of the time during each 20-day bracket.
%ClimaCutoff = .5 
fastice_clima_nocutoff=fastice_clima;
fastice_clima(fastice_clima >= ClimaCutoff) = 1; %only considered fast ice if it exists for X% of the each 20-day period (X%=ClimaCutoff) (e.g. ice must exist for 50% of the 20-day bracket each year in teh climatology).
fastice_clima(fastice_clima < ClimaCutoff) = 0;
res = 1;
fastice_mean =squeeze(nanmean(fastice_clima(:,1:res:end,1:res:end),1)); % This gives the % mean coverage over the 2000-2008 dataset.
%% Set a second cutoff (e.g. 0.9) which says that we will include fastice if it exists for 90% of the 2000-2008 climatology
%InterannualCutoff = .9;
fastice_meanMask = fastice_mean;
fastice_meanMask(fastice_meanMask >= InterannualCutoff) = 1; % if greater than InterannualCutoff then fastice in Mask. (e.g. ice must exist for 90% of the year to be fastice).
fastice_meanMask(fastice_meanMask < InterannualCutoff) = 0;

%% plotting
figure('pos',[100 100 1000 1000]);
flat_pcolor(lons(1:res:end,1:res:end),lats(1:res:end,1:res:end),fastice_mean');clim=caxis;
hold on,
contour(lons(1:res:end,1:res:end),lats(1:res:end,1:res:end),fastice_meanMask',[1 1],'k--','LineWidth',2);
caxis(clim)
colorbar,
xlabel('Longitude (^oE)','FontSize',15),
ylabel('Latitude (^oN)','FontSize',15);
set(gca,'FontSize',13,'FontWeight','bold');
xlim([104 130]);
ylim([-67.1 -65]);
strtit={['Mean % fast ice cover from 2000-2008 climatology assuming '],[num2str(ClimaCutoff*100),'% coverage in each 20-day bracket',' (1 == always present)'],['Contour shows limit of ',num2str(InterannualCutoff*100),'% cover for mean distribution']};
ntitle(strtit,'FontSize',16),

%save '/ds/projects/iomp/totten/ana/dgwyther/geom/tisom_fastice_raw.mat' lons lats fastice_mean

figure('pos',[200 200 1000 1000])
flat_pcolor(lons,lats,fastice_mean'),
axis([104 130 -67.5 -64]),
xlabel('Longitude (^oE)','FontSize',15),
ylabel('Latitude (^oN)','FontSize',15),
colorbar,set(gca,'FontSize',13,'FontWeight','bold'),
ntitle('Yearly persistence map of Fast Ice (2000 - 2008)','FontSize',16),

disp('output data for ROMS would be the variable fastice_meanMask')
