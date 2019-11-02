function [geoid_p]=gravHawaii(a,C)
%%% this function will calculate the geoid given, a the size of the 
% cylinders to approximate over and C, longitude and latitude
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
f2 = zeros(xsize,ysize);
gmars = 3.711;
rho = a(2); %density difference of the core to the land
AK = a(1)^2;
Kdepth2 = 21000;
Kdepth1 = -60000;
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
        geoid_p(xx,yy) = sum(sum(geoid));
    end
end
geoid_p = reshape(geoid_p,[xsize*ysize,1]);
end

