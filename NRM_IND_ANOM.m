%Script to compute magnetic anomalies of a block of cratonic mantle;
%For Siberian craton;
%using a Monte-Carlo resampling method;
%by Idoko C. M. 2017
tic
%%%Basic constants
deg2rad = pi/180;
rad2deg = 1/deg2rad;

%block of craton constructed at lat1 now at lat 2;
lat1 = -90:90; % latitudes of magnetization of an old continental lithos;
lat2 = 63;%Present Earth field latitude central Sib. craton;

%Anoamly due to Remanent Magnetization;
% observation points in km
x0 = -3000:30:3000;
y0 = 0*ones(size(x0));
z0 = -460*ones(size(x0));%Positive downwards, negaticve - vertical upwards;
nx0 = length(y0);

%Discretization parameters of a bloack of craton of dimensions: 3000 km by
%2000 km by 105 km;
xmin = -1500; xmax = 1500;dx = 300;
ymin = 0; ymax = 2000;dy = 400;
zmin = 40; zmax = 105;dz = 65;%Considering Moho = 40km,Curie depth = 105km;
nxbox = (xmax-xmin)/dx;
nybox = (ymax-ymin)/dy;
nzbox = (zmax-zmin)/dz;
nbox = nxbox*nybox*nzbox;
ibox = 0;%counts initial box;
store = zeros(nbox,3,2);% creates a 3d matrix for the box parameters

%%Below is a routine that assigns co-ordinatetes to discretized boxes%%
for iz = 1:nzbox%loop over z
    z1 = zmin+(iz-1)*dz;
    z2 = z1+dz;
    for ix = 1:nxbox%loop over x
    x1 = xmin+(ix-1)*dx;
    x2 = x1+dx;
        for iy = 1:nybox%loop over y
        y1 = ymin+(iy-1)*dy;
        y2 = y1+dy;
        ibox = ibox+1;
        store(ibox,1,:) = [x1 x2];
        store(ibox,2,:) = [y1,y2];
        store(ibox,3,:) = [z1,z2];
        end
    end
end

%%%Parameters for calculation of anomaly due to remanence;
fi = rad2deg*atan(2*tan(lat2*deg2rad));%Ambient Field inclination;
fd = 16;% Average Ambent field declin.;
tnrm = zeros(nbox,nx0);%presets values for remanent anomalies;
nmodel = 1:100;%Number of Monte-Carlo Models;
nNM = length(nmodel);
Tnrm = zeros(nNM,nx0);

%Loop to compute the magnetic anomalies due to discretized blocks;
%over the observation points and repeated for nNM number of times;
for in = 1:nNM
 m = lognrnd(-2,0.9,1,nbox);%sets a random nrm folowing a lognormal;  
 mi = [90*rand(1,nbox/2),-90*rand(1,nbox/2)];% random mag inclin. from -90 to 90;
 mi = mi(randperm(length(mi)));
md = 360*rand(1,nbox);%random magnetic field delin. for nboxes ;
theta = 90*ones(1,nbox);%Rotates the cordinates to recommended axis for mbox;
for iob = 1:nx0%loop over observation
    for ibox = 1:nbox%loop over boxes
        %Picks out the appropriate cordinate for each box;
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
       %For NRM
       %Computes the magnetic anomaly due to each block;
        t1 = (mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m(ibox(1,1)),theta(ibox(1,1)))...
           + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m(ibox(1,1)),theta(ibox(1,1))));
        tnrm(ibox,iob) = t1;%stores the results of each t1;
    end 
    tnrm = sum(tnrm);%Adds the results for all the observations in the boxes;
end
Tnrm(in,:) = tnrm;%Stores anomalies from all Monte-Carlo calculations;
end

% Anomaly Due to Induced Magnetization
%%observation points
x0 = -3000:30:3000;
y0 = 0*ones(size(x0));
z0 = [-460*ones(size(x0));-460*ones(size(x0));-460*ones(size(x0));-460*ones(size(x0))];
nx0 = length(y0);
%discretization parameters
xmin = -1500;xmax = 1500;dx = 300;
ymin = 0; ymax = 2000;dy = 400;
zmin = 40; zmax =160;dz = 30; %Moho = 40km. 160km thick lithosphere;
nxbox = (xmax-xmin)/dx;
nybox = (ymax-ymin)/dy;
nzbox = (zmax-zmin)/dz;
nbox = nxbox*nybox*nzbox;
ibox = 0;%counts initial box;
store = zeros(nbox,3,2);% creates a 3d matrix for the box parameters

%Assigns co-ordinatetes to discretized boxes
for iz = 1:nzbox%loop over z
    z1 = zmin+(iz-1)*dz;
    z2 = z1+dz;
    for ix = 1:nxbox%loop over x
    x1 = xmin+(ix-1)*dx;
    x2 = x1+dx;
        for iy = 1:nybox%loop over y
        y1 = ymin+(iy-1)*dy;
        y2 = y1+dy;
        ibox = ibox+1;
        store(ibox,1,:) = [x1 x2];
        store(ibox,2,:) = [y1,y2];
        store(ibox,3,:) = [z1,z2];
        end
    end
end
nxybox = nxbox*nybox;
nbox = nxbox*nybox*nzbox;
ibox = 0;%counts initial box;

%%Produce induced magnetization matrx (Mstore) for each layer of lithosphere under
%%consideration. Each layer has a maximum possible value of mag, and a min.
%%possible vaule of magnetization;
Mstore = zeros(nzbox,nxybox);
Mmax = [0.33, 0.023, 0.017, 0.007];
Mmin = [0.023, 0.017, 0.007, 0.003];
Mx = zeros(nzbox,nxybox);
nM = length(Mmax);%Number of blocks in each litho;
for im = 1:length(Mmax)%Loop over layers
    for jm = 1:nxybox%Loop over all boxes in each layer;
         Mx = (Mmax(im)-Mmin(im)).*lognrnd(-3,0.7,1,nxybox)+Mmin(im);
    end
    Mstore(im, :) = Mx;
end
m = Mstore;%magnetization;
fi = rad2deg*atan(2*tan(lat2*deg2rad));%Ambinet field inclination;
mi = fi*ones(1,nxybox);%Mag incl. of boxes same as ambient field incl.;
fd = 16;%Ambinet field decl.
md = fd*ones(1,nxybox);%Mag field decl. for boxes same as ambient field decl;
theta = 90*ones(1,nxybox);%Rotates my cor-dinates to mbox-recommended cord;

%%Anomaly due to induced mag for each layer of the litho;
tind1 = zeros(nxybox,nx0);
tind2 = zeros(nxybox,nx0);
tind3 = zeros(nxybox,nx0);
tind4 = zeros(nxybox,nx0);
for iob = 1:nx0%loop over observation
    for ibox = 1:nxybox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        %Layer 1
        %t1ind computes the ind mag field anom for all blocks at one obsv;
         t1ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m(1,ibox),theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m(1,ibox),theta(ibox(1,1)));
        %Second layer;
        t2ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m(2,ibox),theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m(2,ibox),theta(ibox(1,1)));
        %Layer2;
        t3ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m(3,ibox),theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m(3,ibox),theta(ibox(1,1)));
        %layer 3;
        t4ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m(4,ibox),theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m(4,ibox),theta(ibox(1,1)));
        %Results;
         tind1(ibox,iob) = t1ind;
         tind2(ibox,iob) = t2ind;
         tind3(ibox,iob) = t3ind;
         tind4(ibox,iob) = t4ind;
         
    end
   %Adding all box results for each observation point 
    tind1 = sum(tind1);
    tind2 = sum(tind2);
    tind3 = sum(tind3);
    tind4 = sum(tind4);
end
%
Tind = tind1+tind2+tind3+tind4;%Sums ind. magn. anomaly due to all layers;
Tind1 = repmat(Tind,nNM,1);%Replicates this to nNM calculattions of nrm anom;
T = Tnrm+Tind1;%sum of anomalies due to induced and remanence;

%Visualization of all Monte-Carlo models
for b = 1:nNM
   figure(27)
   plot(x0,T(b,:))
   hold on
end
%Adding an envelope to the Monte Carlo Models
envelope = zeros(1,nx0);
for i = 1:nx0
    envelope(i) = 2*std(T(:,i));
end

%Visualization of one model results;
figure(17)
hold on
plot(x0,Tnrm(14,:),'b')
plot(x0,Tind,'r')
T14 = Tnrm(14,:)+Tind;
plot(x0,T14,'y')
xlabel('Horizontal distance(km)','FontSize',28,'FontWeight','bold'),ylabel('Total Magnetic Anomaly(nT)','FontSize',28,'FontWeight','bold')
title('Magnetic anomalies from Siberian at 460km altitude','FontSize',28,'FontWeight','bold')
toc