
%[ASAID_islands_HL] = textread('/ds/projects/iomp/obs/ASAID/ASAID_island_HL.txt','headerlines',1,'delimiter',' ');

[ASAID_islands_HL_x,ASAID_islands_HL_y,ASAID_islands_HL_lat,ASAID_islands_HL_lon]=textread('/ds/projects/iomp/obs/ASAID/ASAID_island_HL.txt','%n%n%n%n%*[^\n]','headerlines',1);

