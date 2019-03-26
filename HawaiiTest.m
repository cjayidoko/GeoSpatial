%Script to compute magnetic anomalies due to a hotspot;
%From Hawaii hotspot;
%using a Monte-Carlo method;
%by Idoko C.M 2017.
tic
%Basic constants;
deg2rad = pi/180;
rad2deg = 1/deg2rad;

%block of lithosphere constructed lat1 now at lat 2;
lat2 = 24;% Earth field latitude

%Magnetic anomaly due to remanence;
% Observation points
x0 = -4000:40:4000;
y0 = 0*ones(size(x0));
z0 = -460*ones(size(y0));
nx0 = length(y0);
%%Discretization parameters of a block of ocenaic mantle of dimensions: 6000 km by
%1500 km by 20 km. Hotspot lie within, dimensions: 3000 km by 500km by 20km;
xmin = -3000;xmax = 3000;dx = 300;
ymin = -750; ymax = 750;dy = 250;
zmin = 20; zmax = 40;dz = 10;%Moho = 20km, Curie depth = 40km away from swell;
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
mpackstore1 = zeros(1,nbox);%stores the magn. values used in the nrm model;
fi = rad2deg*atan(2*tan(lat2*deg2rad));%Ambient Field inclination;
theta = 90*ones(1,nbox);%Rotates to the desired direction.
tnrm = zeros(nbox,nx0);%presets values for remanent anomalies;
nmodel = 1:100;%Number of Monte-Carlo Models;
nNM = length(nmodel);
Tnrm = zeros(nNM,nx0);
for in = 1:nNM
 mpack = lognrnd(-2,0.9,1,nbox);%sets a random nrm folowing a lognormal;
 mi = [(resample(inclin(1:15),250,1)),resample(inclin(-1:-1:-15),250,1)];% random mag inclin. acroos lat2;
 mi = mi(randperm(nbox));%Magnetic inclination randonmized for each block
fd = 10;%Ambent field declin.;
%8*rand(1,nbox)+5 sets random values from +5 to 5+8;
%8*rand(1,nbox)-5 sets random values from -5 to -5+8;
md = 10*rand(1,nbox)+3.3;%2*rand(1,nbox/2)-9];%random magnetic field delin. for nboxes ;
md = md(randperm(length(md)));
for iob = 1:nx0%loop over observation
    for ibox = 1:nbox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        if x1 >= -1500 && x1 < 1500 && y1 >= -250 && y1 < 250
            m = 0;
        else
            m = mpack(ibox);
        end
        t1 = (mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m,theta(ibox(1,1)))...
           + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m,theta(ibox(1,1))));
        tnrm(ibox,iob) = t1;%stores the results of each t1;
        mpackstore1(1,ibox) = m;
    end 
    tnrm = sum(tnrm);%Adds the results for all the observations in the boxes;
end
Tnrm(in,:) = tnrm;
end

%Magnetic anomaly due to Induced Magnetization;

%%Produce induced magnetization matrx (Mstore) for each layer of mantle under
%%consideration. Each layer has a maximum possible value of magnetization (m), and a min.
%%max possible value of magnetization;
nxybox = nxbox*nybox;
Mstore = zeros(nzbox,nxybox);
Mmax = [0.025, 0.025, 0.025, 0.025];
Mmin = [0.005, 0.005, 0.005, 0.007];
Mx = zeros(nzbox,nxybox);
nM = length(Mmax);
for im = 1:nM;
    for jm = 1:nxybox;
         Mx = (Mmax(im)-Mmin(im)).*lognrnd(-2,0.9,1,nxybox)+Mmin(im);
    end
    Mstore(im, :) = Mx;
end
mpack1 = Mstore(1,:);
mpack2 = Mstore(2,:);
mpack3 = Mstore(3,:);
mpack4 = Mstore(4,:);
mpackstore11 = zeros(1,nxybox);%stor the m values used in the Ind m model;
mpackstore2 = zeros(1,nxybox);%stor the m values used in the ind m model;
mpackstore3 = zeros(1,nxybox);%stor the m values used in the ind m model;
mpackstore4 = zeros(1,nxybox);%stor the m values used in the ind m model;
fi = rad2deg*atan(2*tan(lat2*deg2rad));%Ambient Field inclination;
mi = fi*ones(1,nxybox);% random mag inclin. from 0 to 90;
fd = 10;%Ambent field declin.;
md = fd*ones(1,nxybox);%random magnetic field delin. for nboxes ;
%Initialize values for parameters;
theta = 90*ones(1,nxybox);
tind1 = zeros(nxybox,nx0);
tind2 = zeros(nxybox,nx0);
tind3 = zeros(nxybox,nx0);
%Computing Induced magnetic anomaly from the first layer;
for iob = 1:nx0%loop over observation
    for ibox = 1:nxybox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        m1 = mpack1(ibox);
        %t1ind computes the ind mag field anom for all blocks at one obsv;
        if x1 >= -1500 && x1 < 1500 && y1 >= -250 && y1 < 250
            m2 = 0;
        else
            m2 = mpack2(ibox);
        end
            %
         t1ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m1,theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m1,theta(ibox(1,1)));
         tind1(ibox,iob) = t1ind;
         mpackstore11(1,ibox) = m1;
    end
   %Adding all box results for each observation point 
    tind1 = sum(tind1);
end
%%Computing Induced magnetic anomaly from the second layer;
for iob = 1:nx0%loop over observation
    for ibox = 1:nxybox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        if x1 >= -1500 && x1 < 1500 && y1 >= -250 && y1 < 250
            m2 = 0;
        else
            m2 = mpack2(ibox);
        end
            %Second layer;
        t2ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m2,theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m2,theta(ibox(1,1)));
         tind2(ibox,iob) = t2ind;
         mpackstore2(1,ibox) = m2;
    end
   %Adding all box results for each observation point 
    tind2 = sum(tind2);

end
%%Computing Induced magnetic anomaly from the third layer;
for iob = 1:nx0%loop over observation
    for ibox = 1:nxybox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        if x1 >= -1500 && x1 < 1500 && y1 >= -250 && y1 < 250
            m3 = 0;
        else
            m3 = mpack3(ibox);
        end
        t3ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m3,theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m3,theta(ibox(1,1)));

         tind3(ibox,iob) = t3ind;
         mpackstore3(1,ibox) = m3;
  
    end
   %Adding all box results for each observation point 
    tind3 = sum(tind3);
 end
%Computing Induced magnetic anomaly from the fourth layer;
tind4 = zeros(nxybox,nx0);
for iob = 1:nx0%loop over observation
    for ibox = 1:nxybox%loop over boxes
        x1 = store(ibox,1,1);
        x2 = store(ibox,1,2);
        y1 = store(ibox,2,1);
        y2 = store(ibox,2,2);
        z1 = store(ibox,3,1);
        z2 = store(ibox,3,2);
        if x1 >= -1500 && x1 < 1500 && y1 >= -250 && y1 < 250
            m4 = 0;
        else

            m4 = mpack4(ibox);
        end
        %Layer 1
        %t1i
        t4ind = mbox(x0(iob),y0(iob),z0(iob),x1,y1,z1,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,m4,theta(ibox(1,1)))...
            + mbox(x0(iob),y0(iob),z0(iob),x1,y1,z2,x2,y2,mi(ibox(1,1)),md(ibox(1,1)),fi,fd,-m4,theta(ibox(1,1)));
        
         tind4(ibox,iob) = t4ind;
         mpackstore4(1,ibox) = m4;
    end
   %Adding all box results for each observation point 
       tind4 = sum(tind4);
end

%Combining all the induced contributions;
Tind = tind1+tind2+tind3+tind4;
%Visualizing the induced and remanent contributions;
figure(32)
hold on
plot(x0,Tnrm(7,:),'r^')
plot(x0,Tind,'b*')
xlabel('Horizontal distance(km)','FontSize',28,'FontWeight','bold'),ylabel('Total Magnetic Anomaly(nT)','FontSize',28,'FontWeight','bold')
title('Magnetic anomalies from Hawaii hotspot 460km Alt','FontSize',28,'FontWeight','bold')
%Visualizing Monte-Carlo models;
Tind = repmat(Tind,nNM,1);
T = Tnrm+Tind;
for b = 1:nNM
   figure(22)
   plot(x0,T(b,:))
   hold on
end
%Adding Envelop to the monte-Carlo models;
envelope1 = zeros(1,nx0);
envelope2 = zeros(1,nx0);
for i = 1:nx0
        envelope1(i) = 2*std(T(:,i));
        envelope2(i) = -5*std(T(:,i));
end
figure(5)
plot(x0,envelope1,'--m')
figure(6)
plot(x0,envelope2,'--m')
xlabel('Horizontal distance(km)','FontSize',28,'FontWeight','bold'),ylabel('Total Magnetic Anomaly(nT)','FontSize',28,'FontWeight','bold')
title('Magnetic anomalies from Hawaii hotspot 460km Alt','FontSize',28,'FontWeight','bold')
toc