%%% A.Tannenbaum %%%%
%%%% The purpose of this code is to comapre the gravity anomaly under Mauna
% Kea in Hawaii to that of Olympus Mons on Mars. The theory is that Olympus
% Mons and Mauna Kea both have some sort of dense root, potentially made of
% Olivine
clear, clc, close all
f1 = figure(1);

%%% topography NASA mesh plot %%%
%set up latitude and longitude min and max to plot over
 ymin=10.1 ;
 ymax=30; %degrees latitude and
 xmin=210;
 xmax=240;
 
 %create mesh points
 [xi, yi] = meshgrid([xmin:.1:xmax],[ymin:.1:ymax]); 
 %load NASA data set
 load('Topy') 
 x=Topy(:,1);y=Topy(:,2);z=Topy(:,3);
 %calculate topography
 zib=griddata(x,y,z,xi,yi,'cubic');
 %mesh plot
 mesh(xi,yi,zib)
 hold on 
 %add plot components
 colormap(landcolor)
 %varargout = dem(x, y, z, varorigin);
 %Radius of mars
 R = 3389500; %m
 X=R*(xi-min(min(xi)))*pi/180.*cosd(20); %average latitude for olympus mons
 Y=R*(yi-min(min(yi)))*pi/180; 
 Zt=zib;
 label1('Longitude (degrees)', 'Latitude (degrees)', 'Height (m)', 'NASA Topography of Olympus Mons', 12)
 
 %%% gravity NASA meshplot %%%
 f2 = figure(2);
 ymin=10.1 ;
 ymax=30; %degrees latitude and
 xmin=210;
 xmax=240;
 
 %meshgrid
 [xi1, yi1] = meshgrid([xmin:.1:xmax],[ymin:.1:ymax]); 
 load('Gravy') 
 x=Gravy(:,1);y=Gravy(:,2);z=Gravy(:,3); 
 zi1=griddata(x,y,z,xi1,yi1,'cubic');
 %Radius of mars
 R = 3389500; %m
 mesh(xi1,yi1,zi1)
 hold on 
 colormap(landcolor)
 %varargout = dem(x, y, z, varorigin);
 X=R*(xi1-min(min(xi1)))*pi/180.*cosd(20); %average latitude for olympus mons
 Y=R*(yi1-min(min(yi1)))*pi/180; 
 Zg=zi1;
 label1('Longitude (degrees)', 'Latitude (degrees)', 'Gravity (mGals)', 'NASA Gravity of Olympus Mons', 12)
 
 % Making the DEM: 
 % Topography:
 f8 = figure(8);
 dem(xi(1,:), yi(:,1), zib, 3);
 label('Longitude (degrees)', 'Latitude (degrees)', 'DEM: NASA Topography of Olympus Mons', 12)
      
% Gravity:
f9 = figure(9);
dem(xi1(1,:), yi1(:,1), zi1, 3);
label('Longitude (degrees)', 'Latitude (degrees)', 'DEM: NASA Gravity of Olympus Mons', 12)
 
%free air anomaly
GFA = 0.3086*(Zt);
 
%Geoid
r_mars = 3389; %km radius
gmars = 3.711; %m/s^2

%%%%% Terrain and Geoid Calculation %%%%%%%
HGT = Zt;
XS = X(HGT == max(max(HGT)));
YS = Y(HGT == max(max(HGT)));
ZS = max(max(HGT));
e = HGT;
Z = max(max(Zt));
[M,N]=size(X);
G = 6.67e-11;
rho = 2900;
A = (Y(2,1) - Y(1,1)) * (X(1,2) - X(1,1));
x_cord = (1:5:M-1); %such that step size doesn't overstep
y_cord = (1:5:N-1);

for l =  1:length(x_cord)
for m =  1:length(y_cord)
%setup coordinates to move in lines along the grid, this will vary location
%of the stations
XS = X(x_cord(l), y_cord(m));
YS = Y(x_cord(l), y_cord(m));
%Save points
Xa(l,m) = XS;
Ya(l,m) = YS;
XX = X-XS;	%%%%	sets x- differential distances	in meters
YY = Y-YS;	%%%%	sets y-differential distances
r=sqrt(XX.^2 + YY.^2);			%%%	sets	radial distances
for	i=1:M
for	j=1:N
z1=0-ZS;
z2=HGT(i,j)-ZS;
% Angles for Geoid Calculations
alpha= atan2(z1,r(i,j));
beta = atan2(z2,r(i,j)); 
gs(i,j)= A*G*rho*1e5*((cos(beta)-cos(alpha))/r(i,j));
geoid2(i,j)=1/2*log( (1+sin(beta))/(1-sin(beta)))-1/2*log((1+sin(alpha))/(1-sin(alpha)));
geoid(i,j)=geoid2(i,j)*G*A*rho/gmars;
end	%%%%	end	of	j-loop
end	%%%%	end	of	i-loop
%%% Telford 2.56 Assume rod is a cylinder of height L and radius R
L=ZS;
R=sqrt(A/pi); %%%%% Radius gives same area as A
[mn,ix]=min(min(r')); %%%%% find the x-index of the closest node
[mn,iy]=min(min(r)); %%%%% find the y index of the closest node
%geoid(ix,iy)=geoid2(i,j)*G*A*rho/gmars;
gs(ix,iy)=2*pi*6.67e-11*rho*1e5*(L+R-sqrt(L^2+R^2)); %% Use cylinder
geoid2(i,j)=1/2*log( (1+sin(beta))/(1-sin(beta)))-1/2*log((1+sin(alpha))/(1-sin(alpha)));
geoid2(ix,iy)=0;
geoid(ix,iy)=geoid2(ix,iy)*G*A*rho/gmars;
%sum over both 
terrain(l,m)=sum(sum(gs));
geoid_p(l,m)=sum(sum(geoid)); %%%%% now sum everything. This is the gravity from the terrain.
end
end

%%% Plot Calcuated Terrain and Geoid %%%
f3 = figure(3);
%terrain correction
ter_cor = terrain - terrain(1,:);
ter_cor = ter_cor - min(min(ter_cor));

%geoid correction
geoid_cor = geoid_p - geoid_p(1,:);
geoid_cor = geoid_cor - min(min(geoid_p));

%plot mesh
Xa = Xa.*10^-3; %convert to km
Ya = Ya.*10^-3; %convert to km
mesh(Xa, Ya ,ter_cor);
label1('Distance (m)', 'Distance (m)','Gravity (mGals)', 'Calculated Terrain of Olympus Mons', 12)
f4 = figure(4);

%correct for geoid by fitting the side of the 
mesh(Xa, Ya ,geoid_cor);
label1('Distance (m)', 'Distance (m)', 'Gravity (mGals)', 'Calculated Geoid of Olympus Mons', 12)

%%% Nasa Geoid %%%
%topography
 f5 = figure(5);
 ymin=10.1 ;
 ymax=30; %degrees latitude and
 xmin=210;
 xmax=240;
 
 %meshgrid
 [xi, yi] = meshgrid([xmin:.1:xmax],[ymin:.1:ymax]); 
 load('Geoidy') 
 x=Geoidy(:,1);y=Geoidy(:,2);z=Geoidy(:,3); 
 zib=griddata(x,y,z,xi,yi,'cubic');
 %Radius of mars
 R = 3389500; %m
 geoid_cor = zib - zib(1,:);
 geoid_cor = geoid_cor - min(min(zib));
 mesh(xi,yi,geoid_cor)
 hold on 
 colormap(landcolor)
 %varargout = dem(x, y, z, varorigin);
 X=R*(xi-min(min(xi)))*pi/180.*cosd(20); %average latitude for olympus mons
 Y=R*(yi-min(min(yi)))*pi/180; 
 Zt=zib;
 label1('Longitude (m)', 'Latitude (m)', 'Height (m)', 'NASA Geoid of Olympus Mons', 12)
 
%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Fit %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

%%% bouguer gravity
%plot bouguer
xsize = 40;
ysize = 60;
f6 = figure(6);
 [xib, yib] = meshgrid([xmin:.5:xmax],[ymin:.5:ymax]); 
 x=Gravy(:,1);y=Gravy(:,2);z=Gravy(:,3); 
 zib=griddata(x,y,z,xib,yib,'cubic');
 zib = zib(:,1:ysize);
 g_b = zib - terrain;
 xib = xib(:,1:ysize);
 yib = yib(:,1:ysize);
 mesh(xib,yib,g_b)
hold on 
colormap(landcolor)
label1('Longitude (degrees)', 'Latitude (degrees)', 'Gravity (mGals)', 'Bouguer Anomaly of Olympus Mons', 12)
%inversion for Bouguer
Xa = Xa./10^-3; %convert to m
Ya = Ya./10^-3; %convert to m
lon = xib;
lat = yib;
C = [lon lat];

%reshape location values to match inputs needed for nlinfit
lon = reshape(lon,[xsize*ysize,1]);
lat = reshape(lon,[xsize*ysize,1]);
C = [lon,lat];
a = [29500, 800];
g_b = reshape(g_b,[xsize*ysize,1]);

%see output with guesses and use nlinfit
f1_d = gravHawaii2(a,C);
[anew,R,J,SIG] = nlinfit(C,g_b,@gravHawaii2,a);

%Plug in guesses to function to plot
f2_d = gravHawaii2(anew, C);
f3_d = reshape(f2_d,[xsize,ysize]);
f3_d = f3_d*1.5;
f7 = figure(7);

%setup mesh plot
mesh(xib,yib,f3_d)
dev=sqrt(diag(SIG));
error = sum(R.^2);
label1('Longitude (degrees)', 'Latitude (degrees)', 'Gravity (mGals)', 'Bouguer Gravity Model of Olympus Mons', 12)

%geoid
[xi, yi] = meshgrid([xmin:.5:xmax],[ymin:.5:ymax]); 
 load('Geoidy') 
 x=Geoidy(:,1);y=Geoidy(:,2);z=Geoidy(:,3); 
 zib=griddata(x,y,z,xi,yi,'cubic');
  xi = xi(:,1:ysize);
 yi = yi(:,1:ysize);
 lon = xi;
lat = yi;

C = [lon lat];
a = [25302, 640];
%function [geoid_p]=gravHawaii(a,C)
%Digitized	the	locations	of	the	peaks	of	Mauna	Loa	and	Kea
XoK=986790; 
YoK=473310;
xsize = 40;
ysize = 60;
%geoid2(xx,yy) = 0;
dx=a(1);			%%%%	use	nlinfit	to	find	this
dy=a(1);            %spacing between nodes
dl = dx;
%y = y.Grav_anon;
load('lat.mat')
load('lon.mat')
load('Xa1.mat')
load('Ya1.mat')
G = 6.67e-11;
f = zeros(xsize,ysize);
%f2 = zeros(xsize,ysize);
gmars = 3.711;
rho = a(2); %density difference of the core to the land
AK = a(1)^2;
Kdepth2 = 21000;
Kdepth1 = -60000;

%loop to calculate Geoid across entire data set
for xx = 1:xsize
    for yy = 1:ysize
        for ii=1:9
            for jj=1:9
                dx=(ii-5).*dl;
                dy=(jj-5).*dl;
                ZS=0;
                rK=sqrt((Xa(xx,yy)-(XoK+dx))^2+(Ya(xx,yy)-(YoK+dy))^2);
                alphaK=atan2(Kdepth1-ZS,rK);
                betaK =atan2(Kdepth2-ZS,rK);
                %geoid
                geoid2(ii,jj)=1/2*log((1+sin(betaK))/(1-sin(betaK)))-1/2*log((1+sin(alphaK))/(1-sin(alphaK)));
                geoid(ii,jj)=geoid2(ii,jj)*G*AK*rho/gmars;
            end
        end
        geoid_pfit(xx,yy) = sum(sum(geoid));
    end
end
geoid_pfit = reshape(geoid_pfit,[xsize*ysize,1]);
anew = [25302, 640];

%gravity calculations across area
g1 = gravHawaii3(anew,C);
g2 = reshape(g1,[xsize,ysize]);
f10 = figure(10);

%setup mesh plot
mesh(xi,yi,g2)
dev=sqrt(diag(SIG));
error = sum(R.^2);
label1('Longitude (degrees)', 'Latitude (degrees)', 'Geoid Height (m)', 'Geoid Contribution from Core (Fitted Parameters)', 12)

%total geoid plot
f11 = figure(11);
ge_tot = g2 + geoid_p;
mesh(xi,yi,ge_tot)
label1('Longitude (degrees)', 'Latitude (degrees)', 'Geoid Height (m)', 'Total Calculated Geoid (Fitted Parameters)', 12)

%save images
f =[f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11];
for i = 1:length(f)
set(f(i), 'Position', [0 0 1400 900])
f(i).PaperPositionMode = 'auto';
print(f(i), '-r250','-dpng',['Tannenbaum-Lab32-', num2str(i)])
end

%two reasons why our geoid looks different: there is no core for our
%calculations and there is a lrge regional affect on NASA. Regional can be
%many different things





