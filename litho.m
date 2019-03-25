tic
%Litho is a code prepared to take the spherical harmonic coefficients SHC;
%obtained from SWARM satelite and forward compute the magnetic fields;
%perculiar to the lithosphere.
%This code has some adaptation from igrf code for satelite magnetic data;
%G and H coefficients should be saved in the standard;
%format as SHC file ;
%GH is the 3rd column of SHC file;
%Degree is the 1st column of SHC file
%Order is the 2nd colunmn of SHC file;

clear('gh')
%%SepGH
%This here will extract g and h coefficients and store them separately;
%in a simpler format as desired for use in magnetic field computation;
% variables gh, order, and degree were copied  out and stored separately 
gh = GH;%gh is the variable containing both g and h - 3rd column of SHC;
order = Order;%The 2nd column in the SHC file;
Ngh = length(gh);
H = zeros(Ngh,1);
G = zeros(Ngh,1);
N = zeros(Ngh,1);
j = 1:2393;
for i = 1:length(gh)
  gh1 = gh(i);
  order1 = order(i);
       if order1<0
          h = gh1;
        H(i,1) = h;
       elseif order1 == 0
         g = gh1;
         h = 1000;
         n1 = Degree(i);
        G(i,1) = g;
        H(i,1) = h;
        N(i,i) =n1;
       elseif order1>0
          g = gh1;
          n1 = Degree(i);
      G(i,1) = g;
      N(i,i) =n1;
       end 
end
G(find(G==0)) = [];%remove all zeros components;
H(find(H==0)) = [];%remove all zero values;
H(find(H==1000)) = 0;%restore all h = 1000 to zeros;
N(find(N==0)) = [];%remove all zero compnents;

%%SortGH.
%This will sort G and H into nrow and mcolums for each n;
lastNM = 0;%sets an initial value to count number of m for a given n;
Nm0 = 0;
n = 16:90;%sets all applicable degree values;
nN = length(n);%Number of degrees available;
maxM = 0:90;%orders of the max applicable degree;
nmaxM = length(maxM);%the length of the orders of max applicable degree ;
GG = 0*ones(length(n),length(maxM));%Matrix of G sorted n-rows, m-columns;
HH = 0*ones(length(n),length(maxM));%Matrix of H sorted n-rows, m-columns;
for k = 1:nN
    n1 = n(k);
    m = 0:n1;
    Nm = length(m);
    for l = 1:Nm
        G1 = G(l+lastNM);
        H1 = H(l+lastNM);
        GG(k,l) = G1;
        HH(k,l) = H1;
    end
     lastNM = Nm+lastNM;
end
Rearth_km = 6371.2;
altitude = Rearth_km+460;



%replicate altitude to have same size a latitude;
altitude = repmat(altitude,size(longitude));

numlat = numel(latitude);
numlon = numel(longitude);
numalt = numel(altitude);

% To compute cos(theta) and sin(theta).
costheta = cos((90 - latitude(:))*pi/180);
sintheta = sin((90 - latitude(:))*pi/180);
r = altitude(:);%Flattens into a vector;
cd = 1;
sd = 0;
% Special case when sin(theta) = 0.
sintheta0 = sintheta == 0;
anysintheta0 = any(sintheta0);
anysinthetanot0 = any(~sintheta0);

% Convert longitude to radians, and store in a vector form.
phi = longitude(:)*pi/180;

%set g and h equal to GG and HH respectively;
g = GG;
h = HH;
nmax = 90;%
% Precalculate cos(m*phi) and sin(m*phi) desired number of times into a
% matrix here:
cosphi = cos(bsxfun(@times, 0:nmax, phi));
sinphi = sin(bsxfun(@times, 0:nmax, phi));

%%% BEGIN MAGNETIC FIELD CALCULATION %%%
% Initialize variables used in for loop below.
Br = 0*ones(size(r));
Bt = 0*ones(size(r));
Bp = 0*ones(size(r));
lastP = 1;
lastdP_1 = 0*ones(size(costheta));
lastdP_2 = 0;

% compute dp from n = 1 to n= 15;
for n = 1 : 15
    %n = N(i);
    m = 0 : n;
    
    % Calculate legendre values. The output of the function will have each;
    % m value going down the rows;
    P = legendre(n, costheta, 'sch').';%'sch' for schmidt normalization;
    
    % Also compute the derivative of the legendre with respect to theta. It
    % is given by a recursive function of both the previous legendre values
    % as well as the previous derivatives. Functionally, it is:
    % dP(0, 0) = 0, dP(1, 1) = cos(theta)
    % dP(n, n) = sqrt(1 - 1/(2n))*(sin(theta)*dP(n-1, n-1) +
    %     cos(theta)*P(n-1, n-1))
    % dP(n, m) = (2n - 1)/sqrt(n^2 - m^2)*(cos(theta)*dP(n-1, m) -
    %     sin(theta)*P(n-1, m)) - sqrt(((n-1)^2 - m^2)/(n^2 - m^2))*
    %     dP(n-2, m)
    dP = [bsxfun(@minus, bsxfun(@times, ...
        (2*n - 1)./sqrt(n^2 - m(1:end-1).^2), ...
        bsxfun(@times, costheta, lastdP_1) - bsxfun(@times, sintheta, ...
        lastP)), bsxfun(@times, sqrt(((n - 1)^2 - m(1:end-1).^2)./...
        (n^2 - m(1:end-1).^2)), lastdP_2)), 0*ones(size(costheta))];
    if n > 1
        dP(:, end) = sqrt(1 - 1/(2*n))*...
            (sintheta.*lastdP_1(:, end) + costheta.*lastP(:, end));
        lastdP_2 = [lastdP_1 0*ones(size(costheta))];
    else
        dP(:, end) = costheta;
        lastdP_2 = lastdP_1;
    end
    lastP = P;
    lastdP_1 = dP;
    if n == 14
        lastdP14 = dP;
    elseif n == 15
        lastdP15 = dP;
    end
        
    
    
end

lastdP_1 = lastdP15;%stores this for use in upcoming computation;

% Compute lithospheric magnetic field from the 16th degree to 90th.
for n = 16 : 90
    %n = N(i);
    m = 0 : n;
    
    % Calculate legendre values. The output of the function has each m
    % value going down the rows;
    P = legendre(n, costheta, 'sch').';
    
    % Compute the derivative of the legendre with respect to theta using
    % dP(0, 0) = 0, dP(1, 1) = cos(theta)
    % dP(n, n) = sqrt(1 - 1/(2n))*(sin(theta)*dP(n-1, n-1) +
    %     cos(theta)*P(n-1, n-1))
    % dP(n, m) = (2n - 1)/sqrt(n^2 - m^2)*(cos(theta)*dP(n-1, m) -
    %     sin(theta)*P(n-1, m)) - sqrt(((n-1)^2 - m^2)/(n^2 - m^2))*
    %     dP(n-2, m)
    dP = [bsxfun(@minus, bsxfun(@times, ...
        (2*n - 1)./sqrt(n^2 - m(1:end-1).^2), ...
        bsxfun(@times, costheta, lastdP_1) - bsxfun(@times, sintheta, ...
        lastP)), bsxfun(@times, sqrt(((n - 1)^2 - m(1:end-1).^2)./...
        (n^2 - m(1:end-1).^2)), lastdP_2)), 0*ones(size(costheta))];
    if n > 1
        dP(:, end) = sqrt(1 - 1/(2*n))*...
            (sintheta.*lastdP_1(:, end) + costheta.*lastP(:, end));
        lastdP_2 = [lastdP_1 0*ones(size(costheta))];
    else
        dP(:, end) = costheta;
        lastdP_2 = lastdP_1;
    end
    lastP = P;
    lastdP_1 = dP;
    
    % Multiply coefficients by proper longitude trigonemetric term.
    gcos = bsxfun(@times, g(n-15,m+1), cosphi(:, m + 1));
    gsin = bsxfun(@times, g(n-15,m + 1), sinphi(:, m + 1));
    hcos = bsxfun(@times, h(n-15,m + 1), cosphi(:, m + 1));
    hsin = bsxfun(@times, h(n-15,m + 1), sinphi(:, m + 1));
    
    % Calculate the magnetic field components as a running sum. Find
    % explicit expressions for these in Global Earth Physics: a Handbook of
    % Physical Constants by Thomas J. Aherns (1995), pg. 49. Link:
    % http://books.google.com/books?id=aqjU_NHyre4C&lpg=PP1&dq=Global%20
    % earth%20physics%3A%20a%20handbook%20of%20physical%20constants&pg=PA49
    % #v=onepage&q&f=false
    % (except equation 6 is missing a required 1/sin(theta) and m; correct
    % equations on page 5 (equations 3a-3c) of:
    % http://hanspeterschaub.info/Papers/UnderGradStudents/
    % MagneticField.pdf)
    a_r = (Rearth_km./r).^(n + 2);
    Br = Br + a_r.*(n+1).*sum((gcos + hsin).*P, 2);
    Bt = Bt + a_r.*sum((gcos + hsin).*dP, 2);
    % Different case when sin(theta) == 0 for phi component.
    if anysinthetanot0
        Bp(~sintheta0) = Bp(~sintheta0) - 1./sintheta(~sintheta0).*...
            a_r(~sintheta0).*sum(bsxfun(@times, m, ...
            (-gsin(~sintheta0, :) + hcos(~sintheta0, :)).*...
            P(~sintheta0, :)), 2);
    end
    if anysintheta0
        Bp(sintheta0) = Bp(sintheta0) - costheta(sintheta0).*...
            a_r(sintheta0).*sum((-gsin(sintheta0, :) ...
            + hcos(sintheta0, :)).*dP(sintheta0, :), 2);
    end
    
end

% Convert from spherical to (x,y,z, ie North,East,and Vetical components;
%First we sort the radial magetic fields into theta rows and phi columns;
v = 0;%initial value of number of longitudes;
Phi = (-180:1:179)*pi/180;
Theta = (-90:1:89)*pi/180;
nphi = length(-180:1:179);
ntheta = length(-90:1:89);
Br1 = 0*ones(ntheta,1);
BR = 0*ones(ntheta,nphi);
Bp1 = 0*ones(ntheta,1);
Bt1 = 0*ones(ntheta,1);
BP = 0*ones(ntheta,nphi);
BT = 0*ones(ntheta,nphi);
for k = 1:nphi
    for l = 1:ntheta
        Br2 = Br(l+v);
        Bp2 = Bp(l+v);
        Bt2 = Bt(l+v);
        Br1(l,1) =Br2;
        Bt1(l,1) = Bt2;
        Bp1(l,1) = Bp2;
    end
    v = 180+v;
    BR(:,k) = Br1;
    BP(:,k) = Bp1;
    BT(:,k) = Bt1;
end
BZ = -BR;
BX = -BT;
BY = BP;
%Compute the total fields due to the lithosphere;
%Inc needs to be computed using the code igrf.m at 460km altitude;
Phi = Phi*180/pi;%convert back to degrees;
Theta = (-90:1:89)*pi/180;
Theta2D = repmat(Theta',1,360);
f = sqrt(BX.*BX + BY.*BY + BZ.*BZ);%TheFeildIntensity
H = sqrt(BX.*BX + BY.*BY);%The horizontal field component;
%F = BZ.*sin(Theta2D)+ H.*cos(Theta2D);%Total field anomaly;
Fa = BZ.*sin(Inc)+ H.*cos(Inc);%Total field anomaly, Inc is the inclinatn
%computed from reference field using the code igrf.m;

%%Visualization
Theta = Theta*180/pi;
figure(10)
pcolor(Phi,Theta,Fa)
shading flat
colormap jet
figure(11)
pcolor(Phi,Theta,BZ)
shading flat
colormap jet
toc
