function map_rivers(GridFile,RiverFile)
% MAP_RIVERS maps the ROMS river locations, showing
% through which cell face the river enters, flowing
% from 'o' to '*' 
% Usage: map_rivers(GridFile,RiverFile);
%  Inputs: GridFile= ROMS grid file
%          RiverFile= ROMS forcing file containing river forcing

if (exist(GridFile, 'file') == 0)
  disp(['The GridFile=' GridFile ' does not exist']);
  error('Please correct');
end;
if (exist(RiverFile, 'file') == 0)
  disp(['The RiverFile=' RiverFile ' does not exist']);
  error('Please correct');
end;

%
% the grid
%

LON_rho=ncread(GridFile,'lon_rho')';
LAT_rho=ncread(GridFile,'lat_rho')';
[eta_rho, xi_rho]=size(LON_rho);
if(exist('x_label')~=1)
  x_label='Longitude (deg)';
end;
if(exist('y_label')~=1)
  y_label='Latitude (deg)';
end;
mask_rho=ncread(GridFile,'mask_rho')';
h=ncread(GridFile,'h')';
mask_rho(~mask_rho)=NaN;
h=h .* -mask_rho;

%
% the rivers
%

nbRiver=length(ncread(RiverFile,'river'));
%t=ncRiver{'river_time'}(:)+2440000;  % time in seconds
t=ncread(RiverFile,'river_time');
Xpos=ncread(RiverFile,'river_Xposition');
Epos=ncread(RiverFile,'river_Eposition');
riv_dir=ncread(RiverFile,'river_direction');
transport=ncread(RiverFile,'river_transport');
size(transport);transport;
if size(transport,2)==1
trans_avg=transport; %for a single time / perpetual forcing
else
trans_avg=mean(transport,2); 
end
size(trans_avg);trans_avg;
% compute "start" (xp1,Epos1) and "stop" (xp2,Epos2) locations
% for lines to indicate river location and direction
LonStar=zeros(0,1);
LatsStar=zeros(0,1);
nbStar=0;
LonDot=zeros(0,1);
LatDot=zeros(0,1);
nbDot=0;
LonLine=zeros(0,2);
LatLine=zeros(0,2);
nbLine=0;
%
% for debugging purposes
DoRight=1;
DoUp=1;
DoLeft=1;
DoDown=1;
%
nbRight=0;
nbUp=0;
nbLeft=0;
nbDown=0;
disp(['nbRiver=' num2str(nbRiver)]);
for iRiver=1:nbRiver
  eRivDir=riv_dir(iRiver,1);
  eTransAvg=trans_avg(iRiver,1);
  nbChoice=0;
  disp(['iRiver=' num2str(iRiver) ...
	'  eRivDir=' num2str(eRivDir) ...
	'  eTransAvg=' num2str(eTransAvg)]);
  if (eRivDir == 0 && eTransAvg > 0)
    nbChoice=nbChoice+1;
    IsRight=1;
  else
    IsRight=0;
  end;
  if (IsRight == 1 && DoRight == 1)
    nbRight=nbRight+1;
    nbStar=nbStar+1;
    iEta1=Epos(iRiver,1)+1;
    iXi1=Xpos(iRiver,1)+1;
    if (iEta1 <= 0 || iEta1 > eta_rho || ...
	iXi1 <= 0 || iXi1 > xi_rho)
      disp('Your river file appear to be wrong');
      disp(['For iRiver=' num2str(iRiver) ...
	    ' we have a right point']);
      disp(['But iEta1=' num2str(iEta1) '  iXi1=' num2str(iXi1)]);
    end;
    LonStar(nbStar,1)=LON_rho(iEta1, iXi1);
    LatStar(nbStar,1)=LAT_rho(iEta1, iXi1);
    iEta2=iEta1;
    iXi2=iXi1-1;
    if (iXi2 > 0)
      nbDot=nbDot+1;
      LonDot(nbDot,1)=LON_rho(iEta2, iXi2);
      LatDot(nbDot,1)=LAT_rho(iEta2, iXi2);
      nbLine=nbLine+1;
      LonLine(nbLine,:)=[LonStar(nbStar,1) LonDot(nbDot,1)];
      LatLine(nbLine,:)=[LatStar(nbStar,1) LatDot(nbDot,1)];
    end;
  end;
  if (eRivDir == 1 && eTransAvg > 0)
    nbChoice=nbChoice+1;
    IsUp=1;
  else
    IsUp=0;
  end;
  if (IsUp == 1 && DoUp == 1)
    nbUp=nbUp+1;
    nbStar=nbStar+1;
    iEta1=Epos(iRiver,1)+1;
    iXi1=Xpos(iRiver,1)+1;
    if (iEta1 <= 0 || iEta1 > eta_rho || ...
	iXi1 <= 0 || iXi1 > xi_rho)
      disp('Your river file appear to be wrong');
      disp(['For iRiver=' num2str(iRiver) ...
	    ' we have a right point']);
      disp(['But iEta1=' num2str(iEta1) '  iXi1=' num2str(iXi1)]);
    end;
    LonStar(nbStar,1)=LON_rho(iEta1, iXi1);
    LatStar(nbStar,1)=LAT_rho(iEta1, iXi1);
    iEta2=iEta1-1;
    iXi2=iXi1;
    if (iEta2 > 0)
      nbDot=nbDot+1;
      LonDot(nbDot,1)=LON_rho(iEta2, iXi2);
      LatDot(nbDot,1)=LAT_rho(iEta2, iXi2);
      nbLine=nbLine+1;
      LonLine(nbLine,:)=[LonStar(nbStar,1) LonDot(nbDot,1)];
      LatLine(nbLine,:)=[LatStar(nbStar,1) LatDot(nbDot,1)];
    end;
  end;
  if (eRivDir == 0 && eTransAvg < 0)
    nbChoice=nbChoice+1;
    IsLeft=1;
  else
    IsLeft=0;
  end;
  if (IsLeft == 1 && DoLeft == 1)
    nbLeft=nbLeft+1;
    nbStar=nbStar+1;
    iEta1=Epos(iRiver,1)+1;
    iXi1=Xpos(iRiver,1);
    if (iEta1 <= 0 || iEta1 > eta_rho || ...
	iXi1 <= 0 || iXi1 > xi_rho)
      disp('Your river file appear to be wrong');
      disp(['For iRiver=' num2str(iRiver) ...
	    ' we have a right point']);
      disp(['But iEta1=' num2str(iEta1) '  iXi1=' num2str(iXi1)]);
    end;
    LonStar(nbStar,1)=LON_rho(iEta1, iXi1);
    LatStar(nbStar,1)=LAT_rho(iEta1, iXi1);
    iEta2=iEta1;
    iXi2=iXi1+1;
    if (iXi2 <= xi_rho)
      nbDot=nbDot+1;
      LonDot(nbDot,1)=LON_rho(iEta2, iXi2);
      LatDot(nbDot,1)=LAT_rho(iEta2, iXi2);
      nbLine=nbLine+1;
      LonLine(nbLine,:)=[LonStar(nbStar,1) LonDot(nbDot,1)];
      LatLine(nbLine,:)=[LatStar(nbStar,1) LatDot(nbDot,1)];
    end;
  end;
  if (eRivDir == 1 && eTransAvg < 0)
    nbChoice=nbChoice+1;
    IsDown=1;
  else
    IsDown=0;
  end;
  if (IsDown == 1 && DoDown == 1)
    nbDown=nbDown+1;
    nbStar=nbStar+1;
    iEta1=Epos(iRiver,1);
    iXi1=Xpos(iRiver,1)+1;
    if (iEta1 <= 0 || iEta1 > eta_rho || ...
	iXi1 <= 0 || iXi1 > xi_rho)
      disp('Your river file appear to be wrong');
      disp(['For iRiver=' num2str(iRiver) ...
	    ' we have a right point']);
      disp(['But iEta1=' num2str(iEta1) '  iXi1=' num2str(iXi1)]);
    end;
    LonStar(nbStar,1)=LON_rho(iEta1, iXi1);
    LatStar(nbStar,1)=LAT_rho(iEta1, iXi1);
    iEta2=iEta1+1;
    iXi2=iXi1;
    if (iEta2 <= eta_rho)
      nbDot=nbDot+1;
      LonDot(nbDot,1)=LON_rho(iEta2, iXi2);
      LatDot(nbDot,1)=LAT_rho(iEta2, iXi2);
      nbLine=nbLine+1;
      LonLine(nbLine,:)=[LonStar(nbStar,1) LonDot(nbDot,1)];
      LatLine(nbLine,:)=[LatStar(nbStar,1) LatDot(nbDot,1)];
    end;
  end;
  disp(['    LonStar=' num2str(LonStar(nbStar,1)) ...
	'  LatStar=' num2str(LatStar(nbStar,1))]);
end;
disp(['nbRight=' num2str(nbRight) ...
      '  nbUp=' num2str(nbUp) ...
      '  nbLeft=' num2str(nbLeft) ...
      '  nbDown=' num2str(nbDown)]);

% plot bathy with river orientation sticks
%pcolor_center_corr(LON_rho,LAT_rho,h);
flat_pcolor(LON_rho,LAT_rho,h);
%pslice2(LON_rho,LAT_rho,h);
lat=nanmean(LAT_rho(:));
xfac=cos(lat*pi/180);
%set(gca, 'DataAspectRatio', [1 xfac 1] );
hold on
plot(LonLine',LatLine','r-');
plot(LonDot,LatDot,'ro');   % arrow "tail"
plot(LonStar,LatStar,'*');    % arrow "head"
for i=1:nbStar
  ts=sprintf(' %d:(%d,%d) %6.1f',i,Xpos(i,1),Epos(i,1),abs(trans_avg(i)));
  text(LonStar(i,1),LatStar(i,1),ts, 'FontSize', 8)
end
hold off;
