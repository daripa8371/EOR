% Quarter five-spot------------------------
global miuo miuw sor0 swr0
miuo = 12.6;
miuw = 1.26; 
%waterflood 
c0 = 0; g0 = 0;
sigma = 10.001./(g0+1)-0.001;
[swr,sor] = compres(sigma,u,v,compvis(c0));
% enpoint mobility ratio
M1 = compmob(1-sor,c0,sor,swr,1)/compmob(0.21,c0,sor,swr,0);

% shock front mobility ratio
La11 = compmob(0.6,c0,sor,swr,1); Lo11 = compmob(0.6,c0,sor,swr,0);
La21 = compmob(0.21,c0,sor0,swr0,1); Lo21 = compmob(0.21,c0,sor0,swr0,0);
M11 = (La11 + Lo11)/ (La21 + Lo21);

fprintf('Waterflood: M1 is %1.3f and M11 is %1.3f \n',M1,M11);

%Polymer flood
c0 = 0.05; g0 = 0;
sigma = 10.001./(g0+1)-0.001;
[swr,sor] = compres(sigma,u,v,compvis(c0));
% enpoint mobility ratio
M5 = compmob(1-sor,c0,sor,swr,1)/compmob(0.21,c0,sor,swr,0);

% shock front mobility ratio
La15 = compmob(0.63,c0,sor,swr,1); Lo15 = compmob(0.63,c0,sor,swr,0);
La25 = compmob(0.21,0,sor0,swr0,1); Lo25 = compmob(0.21,0,sor0,swr0,0);
M55 = (La15 + Lo15)/ (La25 + Lo25);

fprintf('Polymer flood: M5 is %1.3f and M55 is %1.3f \n',M5,M55);


%Surfactant flood
c0 = 0; g0 = 0.05;
sigma = 10.001./(g0+1)-0.001;
[swr,sor] = compres(sigma,u,v,compvis(c0));
% enpoint mobility ratio
M2 = compmob(1-sor,c0,sor,swr,1)/compmob(0.21,c0,sor,swr,0);

% shock front mobility ratio
La12 = compmob(0.73,c0,sor,swr,1); Lo12 = compmob(0.73,c0,sor,swr,0);
La22 = compmob(0.21,0,sor0,swr0,1); Lo22 = compmob(0.21,0,sor0,swr0,0);
M22 = (La12 + Lo12)/ (La22 + Lo22);

fprintf('Surfactant flood: M2 is %1.3f and M22 is %1.3f \n',M2,M22);

miuo = 12.6;
miuw = 6.3; 

%SP flood 1
c0 = 0.001; g0 = 0.05;
sigma = 10.001./(g0+1)-0.001;
[swr,sor] = compres(sigma,u,v,compvis(c0));
% enpoint mobility ratio
M3 = compmob(1-sor,c0,sor,swr,1)/compmob(0.21,c0,sor,swr,0);

% shock front mobility ratio
La13 = compmob(0.75,c0,sor,swr,1); Lo13 = compmob(0.75,c0,sor,swr,0);
La23 = compmob(0.21,0,sor0,swr0,1); Lo23 = compmob(0.21,0,sor0,swr0,0);
M33 = (La13 + Lo13)/ (La23 + Lo23);
fprintf('SP flood 1: M3 is %1.3f and M33 is %1.3f \n',M3,M33);

%SP flood 2
c0 = 0.15; g0 = 0.05;
sigma = 10.001./(g0+1)-0.001;
[swr,sor] = compres(sigma,u,v,compvis(c0));
% enpoint mobility ratio
M4 = compmob(1-sor,c0,sor,swr,1)/compmob(0.21,c0,sor,swr,0);

% shock front mobility ratio
La14 = compmob(0.75,c0,sor,swr,1); Lo14 = compmob(0.75,c0,sor,swr,0);
La24 = compmob(0.21,0,sor0,swr0,1); Lo24 = compmob(0.21,0,sor0,swr0,0);
M44 = (La14 + Lo14)/ (La24 + Lo24);

fprintf('SP flood 2: M4 is %1.3f and M44 is %1.3f \n',M4,M44);