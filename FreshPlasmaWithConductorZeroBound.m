clearvars 
close all

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated2d.gif';




e0=1;%Epsilon Naught
k=1.68637205*10^-10;%Boltzman
meterconv=3.511282*10^-14;%Convert from our units to meters
L=30*.001/meterconv;%Length
DT=1/40*30*.001/meterconv;%DeltaT
NT=2000000;%Iterations
NG=40;%Number of Grid Points
ndensity=2*10^10*meterconv^3;
totalparticles=.01%ndensity*L^3;%Total number of electrons in simulation


N=2000;%Number of Particles
me=totalparticles/N;%Mass of electron
mp=5*totalparticles/N;%mass of proton
e=totalparticles/N;%electron charge

Nin=N;

VX1=0.0;%Velocity Petubation Magnitude
VY1=0.0;%Velocity Petubation Magnitude
mode=1;%wavelengths of sin
rho_back=0;%Background charge density
dx=L/NG;%Delta X
acc=.9999999;%Number close to 1 as to make the Poisson Matrix nonsingular

dE=0;%Initialization of Energy Graph Value
dKE=0;
dUE=0;
dxmom=0;
dymom=0;
rho_sphere1=0;
rho_sphere2=0;
lostenergy=0;
esphere1=0;
esphere2=100;
changeloop=3;


%color
c=[linspace(1,1,N/2) linspace(10,10,N/2)];%Color matrix distinguishing electrons and ions

% initial loading for Particle Positions
%(Commented out some other configurations)
%xp=linspace(0,L-L/N,N)';
%xp=linspace(0,L/2,N)';
xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
%yp=L*rand(N,1)/2;
yp=L*rand(N,1);
%xp=[0;L/2]
%yp=[.05+L/2;.05+L/2]


%Mass Matrix
massmat=mp*ones(1,N);%Make all into ions
parfor n=1:N
    if n<N/2+1
        massmat(n)=me;%Change the first half into electrons
    end
end
%Charge Matrix
emat=e*ones(1,N);%Make all into ions
for n=1:N
    if n<N/2+1
        emat(n)=-e;%Change the first half into electrons
    end
end

Tempin=11700;%Initial Temperature
Temp=11700;%Temperature that will change in time
Con=massmat'/me/k/Temp;%Constant in Maxwellian Distribution
randvtot=(-2*log((1-rand(N,1)))./Con).^.5;%Initialize the speeds into a Maxwellian Distribution
randang=2*pi*rand(N,1);%Choose random direction for velocities
vxp=randvtot.*cos(randang);%put random speed and random direction together
vyp=randvtot.*sin(randang);

%Build Poisson Matrix
v=acc*ones(NG^2,1);%Vector of values close to one
v1=ones(NG^2,1);%Vector of values equal to one
%Put -4 on the diagonal and 1 on surrounding values
Poisson=-4*diag(v1)+diag(v(1:NG^2-1),1)+diag(v(1:NG^2-1),-1)+diag(v(1:NG^2-NG),-NG)+diag(v(1:NG^2-NG),NG);
Poisson=Poisson+diag(v(1:NG),NG^2-NG)+diag(v(1:NG),-NG^2+NG);

%Put ones where there are some missing
for i=1:NG
    Poisson=Poisson+full(sparse(NG*i,NG*i-NG+1,acc,NG^2,NG^2))+full(sparse(NG*i-NG+1,NG*i,acc,NG^2,NG^2));
end
%Remove ones from places where they shouldn't be
for i=1:(NG-1)
    Poisson=Poisson-full(sparse(NG*i,NG*i+1,acc,NG^2,NG^2))-full(sparse(NG*i+1,NG*i,acc,NG^2,NG^2));
end


%Preallocate memory
Ex=zeros(NG^2,1);
Ey=zeros(NG^2,1);
xmom=zeros(NT,1);
ymom=zeros(NT,1);
ForceXSphere1=zeros(NT,1);
ForceYSphere1=zeros(NT,1);
ForceXSphere2=zeros(NT,1);
ForceYSphere2=zeros(NT,1);
totEn=zeros(NT,1);
totUEn=zeros(NT,1);
totKEn=zeros(NT,1);
UdUE=zeros(NT,1);
LdUE=zeros(NT,1);
[X,Y] = meshgrid(dx:dx:L,dx:dx:L);%make grid for rho and phi

bins=10;%Bins for charge around the circle

%sphere1
sxc1=L/2;%Sphere 1 x center
syc1=2*L/5;%Sphere 1 y center
sr1=L/11;%Sphere 1 radius
vsx1=0;%X Velocity of sphere 1
vsy1=0;%Y velocity of sphere 1
ms1=1000*me*totalparticles;%Mass of sphere


sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];%Data to plot the sphere on x
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;%Data to plot sphere on y
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;%Data to plot negative sphere on y

%sphere2
sxc2=L/2;%Sphere 2 x center
syc2=3*L/5;%Sphere 2 y center
sr2=L/22;%Sphere 2 radius
vsx2=0;%X Velocity of sphere 2
vsy2=0;%Y velocity of sphere 
ms2=1000*me*totalparticles;%Mass of sphere

sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];%Data to plot the sphere on x
sy2=(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;%Data to plot sphere on y
nsy2=-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2;%Data to plot negative sphere on y

In1=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);%Particles in Sphere 1
In2=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);%Particles in sphere 2

emat(In1)=0;%Eliminate charge in sphere
emat(In2)=0;%Elimnate charge in sphere

%Capacity Matrix

%Put a particle on every grid point

ypcap=linspace(L/NG,L,NG)';
xpcap=L/NG*ones(NG,1);
xpctemp=L/NG*ones(NG,1);
ypctemp=linspace(L/NG,L,NG)';
for n=2:NG
    ypcap=[ypcap;ypctemp];
    xpcap=[xpcap;xpctemp*n];
end
%Now a particle is on each point


% Put particles on grid points
g1=cat(2,ceil(xp/dx-.5),ceil(yp/dx-.5));%put onto the lower grid points
g=[g1;g1+1];%Put onto lower and higher grid points
fraz1=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower left weight
fraz2=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower right weight
fraz3=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper left weight
fraz4=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper right weight
fraz=[fraz1 fraz2;fraz3 fraz4];%Put it all together now yall

% If it's on a nonexistant grip point, shift it to the other sides
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;



Inc1=(((xpcap-sxc1).^2+(ypcap-syc1).^2)<sr1^2);%Particles on sphere 1
Inc2=(((xpcap-sxc2).^2+(ypcap-syc2).^2)<sr2^2);%Particles on sphere 2
ematc1=zeros(NG^2,1);%Set charges to zero
ematc2=zeros(NG^2,1);%Set charges to zero
ematc1(Inc1)=1;%Charges in sphere are 1
ematc2(Inc2)=1;%Charges in sphere are one

indices1=find(ematc1);%Indexes of points in sphere
Bprime1=zeros(length(indices1));%Initialize Bprime
parfor n=1:length(indices1)
    rhotemp1=zeros(NG^2);
    rhotemp1(indices1(n))=1/dx^2;%Put charge 1 on each point in sphere
    Phitemp1=Poisson\(-rhotemp1*dx^2)/e0;
    Bprime1(:,n)=Phitemp1(indices1);
end
Cap1=inv(Bprime1);%Final Capacity Matrix

indices2=find(ematc2);%Indexes of points in sphere
Bprime2=zeros(length(indices2));%Initialize Bprime
parfor n=1:length(indices2)
    rhotemp2=zeros(NG^2);
    rhotemp2(indices2(n))=1/dx^2;%Put charge 1 on each point in sphere
    Phitemp2=Poisson\(-rhotemp2*dx^2)/e0;
    Bprime2(:,n)=Phitemp2(indices2);
end
Cap2=inv(Bprime2);%Final Capacity Matrix

%Circular Matrix representing 1 particle around Sphere 1
NCirc=100;%Number of particles for circle
angle=linspace(2*pi/NCirc,2*pi,NCirc)';%angle matrix to distribute charges
xpc=sxc1+sr1*sin(angle);
ypc=syc1+sr1*cos(angle);

charge=e/NCirc;

g1=cat(2,ceil(xpc/dx-.5),ceil(ypc/dx-.5));%put onto the lower grid points
g=[g1;g1+1];%Put onto lower and higher grid points
fraz1=((1-abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((1-abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%lower left weight
fraz2=((abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((1-abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%lower right weight
fraz3=((1-abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%upper left weight
fraz4=((abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%upper right weight
fraz=[fraz1 fraz2;fraz3 fraz4];%Put it all together now yall

% If it's on a nonexistant grip point, shift it to the other sides
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;

for n=1:NCirc
    if n==1
        weight=sparse([index2(g(1,1),g(1,2),NG) index2(g(1,1),g(NCirc+1,2),NG) index2(g(NCirc+1,1),g(1,2),NG) index2(g(NCirc+1,1),g(NCirc+1,2),NG)],[1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1)] ,NG^2,NCirc);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index2(g(n,1),g(n,2),NG) index2(g(n,1),g(NCirc+n,2),NG) index2(g(NCirc+n,1),g(n,2),NG) index2(g(NCirc+n,1),g(NCirc+n,2),NG)],[n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n)],NG^2,NCirc);%matrix stating weights of each particle on each grid point
    end
end

%Make the 1 particle sphere matrix
rho_sphere1=sum(charge.*weight(:,1:NCirc),2)/dx^2;

%Circular Matrix representing 1 particle around Sphere 2
xpc=sxc2+sr2*sin(angle);
ypc=syc2+sr2*cos(angle);

g1=cat(2,ceil(xpc/dx-.5),ceil(ypc/dx-.5));%put onto the lower grid points
g=[g1;g1+1];%Put onto lower and higher grid points
fraz1=((1-abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((1-abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%lower left weight
fraz2=((abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((1-abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%lower right weight
fraz3=((1-abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%upper left weight
fraz4=((abs(xpc/dx-g1(1:NCirc,1)+.5))'*diag((abs(ypc/dx-g1(1:NCirc,2)+.5))'))';%upper right weight
fraz=[fraz1 fraz2;fraz3 fraz4];%Put it all together now yall

% If it's on a nonexistant grip point, shift it to the other sides
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;

for n=1:NCirc
    if n==1
        weight=sparse([index2(g(1,1),g(1,2),NG) index2(g(1,1),g(NCirc+1,2),NG) index2(g(NCirc+1,1),g(1,2),NG) index2(g(NCirc+1,1),g(NCirc+1,2),NG)],[1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1)] ,NG^2,NCirc);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index2(g(n,1),g(n,2),NG) index2(g(n,1),g(NCirc+n,2),NG) index2(g(NCirc+n,1),g(n,2),NG) index2(g(NCirc+n,1),g(NCirc+n,2),NG)],[n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n)],NG^2,NCirc);%matrix stating weights of each particle on each grid point
    end
end

%Make the 1 particle sphere matrix
rho_sphere2=sum(charge.*weight(:,1:NCirc),2)/dx^2;


%color
c1=[linspace(1,1,N/2)];%Ions
c2=[linspace(10,10,N/2)];%Electrons
c=[c1 c2];%&PUT it all together

% Main computational cycle
for it=1:NT
% update xp
scatter(xp,yp,20,c,'filled') 

hold on
sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];
sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];
sy2=real((sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);
nsy2=real(-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);
plot(sx1,real(sy1),'r')
plot(sx1,real(nsy1),'r')
plot(sx2,real(sy2),'r')
plot(sx2,real(nsy2),'r')
axis([L/NG L-L/NG L/NG L-L/NG])

    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if it == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
hold off
xp=xp+vxp*DT*.5;
yp=yp+vyp*DT*.5;
syc1=syc1+vsy1*DT*.5;
syc2=syc2+vsy2*DT*.5;
sxc1=sxc1+vsx1*DT*.5;
sxc2=sxc2+vsx2*DT*.5;



%Particle Delete if in a bordering grid space or out of bounds
out=xp<L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=xp>L-L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=yp<L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

out=yp>L-L/NG;
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');

%Add Fresh Plasma
NGnew=4*NG-4;%Border Grid spaces getting replaced
Nnew=floor(NGnew*Nin/NG^2+.5);%Put in a number of particles to keep N about the same
randmat=randi([1 4],Nnew,1);%Which side will they appear?
in1=sum(randmat(:) == 1);%Put on Bottom
xp=[xp;(L-L/NG)*rand(in1,1)];
yp=[yp;(L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 2);%Put on right
xp=[xp;L-(L/NG)*rand(in1,1)];
yp=[yp;(L-L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 3);%Put on top
xp=[xp;L/NG+(L-L/NG)*rand(in1,1)];
yp=[yp;L-(L/NG)*rand(in1,1)];
in1=sum(randmat(:) == 4);%Put on left
xp=[xp;(L/NG)*rand(in1,1)];
yp=[yp;L/NG+(L-L/NG)*rand(in1,1)];

randvtot=(-2*log((1-rand(Nnew,1)))./(((mp-me)*randmat+me)/me/k/Tempin)).^.5;%Maxwellian
randang=2*pi*rand(Nnew,1);%Random direction
vxp=[vxp;randvtot.*cos(randang)];
vyp=[vyp;randvtot.*sin(randang)];

randmat=randi([0 1],Nnew,1);%Ion or Electron?

emat=[emat 2*e*randmat'-e];%Set new charges
massmat=[massmat (mp-me)*randmat'+me];%Set new masses
c=[c 9*randmat'+1];%Set correct colors

N=N+Nnew;%Add new particles to N

% Put particles on grid points
g1=cat(2,ceil(xp/dx-.5),ceil(yp/dx-.5));%put onto the lower grid points
g=[g1;g1+1];%Put onto lower and higher grid points
fraz1=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower left weight
fraz2=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))'))';%lower right weight
fraz3=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper left weight
fraz4=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))'))';%upper right weight
fraz=[fraz1 fraz2;fraz3 fraz4];%Put it all together now yall

% If it's on a nonexistant grip point, shift it to the other sides
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;

    %Put the weights into a 3d matrix, where each row in the third
    %dimension is a new particle. Adding all allements on each row should
    %equal 1 in 4 adjacent indices
for n=1:N
    if n==1
        weight=sparse([index2(g(1,1),g(1,2),NG) index2(g(1,1),g(N+1,2),NG) index2(g(N+1,1),g(1,2),NG) index2(g(N+1,1),g(N+1,2),NG)],[1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1)] ,NG^2,N);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index2(g(n,1),g(n,2),NG) index2(g(n,1),g(N+n,2),NG) index2(g(N+n,1),g(n,2),NG) index2(g(N+n,1),g(N+n,2),NG)],[n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n)],NG^2,N);%matrix stating weights of each particle on each grid point
    end
end

%Compute the charge density
rho=sum(emat.*weight(:,1:N),2)/dx^2+esphere1*rho_sphere1+esphere2*rho_sphere2;%rho of each grid point

for n=1:NG
    if n==1
        rho1=rho(1:NG);
    else
        rho1=cat(2,rho1,rho((n-1)*NG+1:n*NG));
    end
end



drho1tot=zeros(NG^2,1);
drho2tot=zeros(NG^2,1);
for n=1:changeloop
Phi=Poisson\(-rho*dx^2)/e0;

%Find the voltage of sphere 1
PhiC1=sum(Cap1'*Phi(indices1),"all")/sum(Cap1,"all");
%Change in charge
drho1=Cap1'*(PhiC1-Phi(indices1))/dx^2;
drho1tot(indices1)=drho1tot(indices1)+drho1;
rho(indices1)=rho(indices1)+drho1;%New charge

PhiC2=sum(Cap2'*Phi(indices2),"all")/sum(Cap2,"all");%Find the voltage of sphere 2
drho2=Cap2'*(PhiC2-Phi(indices2))/dx^2;%Change in charge
drho2tot(indices2)=drho2tot(indices2)+drho2;
rho(indices2)=rho(indices2)+drho2;%New charge
end



% compute Potential with Poissons equation
Phi=Poisson\(-rho*dx^2)/e0;%New Voltage

for n=1:NG
    if n==1
        Phi1=Phi(1:NG);
    else
        Phi1=cat(2,Phi1,Phi((n-1)*NG+1:n*NG));
    end
end


%Particle Delete if in a sphere
out=(((xp-sxc1).^2+(yp-syc1).^2)<sr1^2);%If in sphere
j=find(out);%Which are out
if isempty(j)==0%Add charge to sphere
    esphere1=esphere1+sum(emat(out));
    weight(:,j)=[];
end

lostenergy=lostenergy+.5*massmat(out)*(vxp(out).^2+vyp(out).^2);%Kinetic energy of particles deleted
%Delete variables
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');%Fix N

out=(((xp-sxc2).^2+(yp-syc2).^2)<sr2^2);%If in sphere
j=find(out);%Which are out
if isempty(j)==0%Add charge to sphere
    esphere2=esphere2+sum(emat(out));
    weight(:,j)=[];
end

lostenergy=lostenergy+.5*massmat(out)*(vxp(out).^2+vyp(out).^2);
%Delete Variable
xp(out)=[];
yp(out)=[];
vxp(out)=[];
vyp(out)=[];
emat(out)=[];
massmat(out)=[];
c(out)=[];
N=N-sum(out,'all');%Fix N



%Electric Field in x direciton
parfor n=1:NG^2
        in1=n+NG;%Index of right node
        in2=n-NG;%Index of left node
        if in1>NG^2
            in1=in1-NG^2;%Fix if on right boundary
        end
        if in2<1
            in2=in2+NG^2;%Fix if on left boundary
        end
        Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in x
end

for n=1:NG
    if n==1
        Ex1=Ex(1:NG);
    else
        Ex1=cat(2,Ex1,Ex((n-1)*NG+1:n*NG));
    end
end

%Electric Field in Y Direction
parfor n=1:NG^2
        in1=n+1;%Index of above node
        in2=n-1;%index of lower lode
        if mod(in1,NG)==1
            in1=in1-NG;%Fix if on the top
        end
        if mod(in2,NG)==0
            in2=in2+NG;%fix if on the bottom
        end
        Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in y
end

for n=1:NG
    if n==1
        Ey1=Ey(1:NG);
    else
        Ey1=cat(2,Ey1,Ey((n-1)*NG+1:n*NG));
    end
end

    %Update Velocity
 vxp=vxp-weight(:,1:N)'*Ex.*emat'*DT./massmat';
 vyp=vyp-weight(:,1:N)'*Ey.*emat'*DT./massmat';
    
 
 ForceXSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1/e)'*Ex;
 ForceYSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1/e)'*Ey;
 ForceXSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2/e)'*Ex;
 ForceYSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2/e)'*Ey;
 
 vsx1=vsx1;%-ForceXSphere1*DT/ms1;%Change sphere 1 x velocity
 vsy1=vsy1;%-ForceYSphere1*DT/ms1;%Change sphere 1 y velocity
 vsx2=vsx2;%-ForceXSphere2*DT/ms2;%Change sphere 2 x velocity
 vsy2=vsy2;%-ForceYSphere2*DT/ms2;%Change sphere 2 y velocity
 
 for n=1:NG
    if n==1
        drho2tot1=drho2tot(1:NG)+esphere2*rho_sphere2(1:NG)/e;
    else
        drho2tot1=cat(2,drho2tot1,drho2tot((n-1)*NG+1:n*NG)+esphere2*rho_sphere2((n-1)*NG+1:n*NG)/e);
    end
 end

 for n=1:NG
    if n==1
        drho1tot1=drho1tot(1:NG)+esphere1*rho_sphere1(1:NG)/e;
    else
        drho1tot1=cat(2,drho1tot1,drho1tot((n-1)*NG+1:n*NG)+esphere1*rho_sphere1((n-1)*NG+1:n*NG)/e);
    end
 end
 
 
%Momentum test (Check to see if momentum is constant
xmom(it)=sum(massmat*vxp);
ymom(it)=sum(massmat*vyp);

%Energy test (Check if Total energy is constant)

totEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all')+.5*dx^2*e0*sum(Ex.^2+Ey.^2,'all');


totEnWLost(it)=totEn(it)+lostenergy;

totUEn(it)=.5*dx^2*e0*sum(Ex.^2+Ey.^2,'all');


%Kinetic Energy
totKEn(it)=.5*sum(massmat*(vxp.^2+vyp.^2),'all');

%Display Iteration
it


Temp=(sum(massmat/me*(vxp.^2+vyp.^2),'all')/N)/2/k;%Calculate the temperature
LD =(k*Temp/ndensity)^.5;
if .5*LD<dx%Compare to delta x, if too small, display error
    "Error, Mesh not fine enough or not hot enough"
    LD
    dx
end
PP=(pi*me*L^2/(N*e^2))^.5;%Compute plasma period 
if .25*PP<DT%If plasma period is small display error
    "Error, large time step"
    PP 
    DT
end

%CFL Condition
C=DT/dx;
if C>1%If time step is larger than dx
    "May be unstable, C>1"
    C
end

end
figure
hold on
plot(totEn)
%plot(dKE+UdUE)
%plot(dKE+LdUE)
plot(totUEn)
plot(totKEn)
%plot(totEnWLost)
%plot(UdUE)
%plot(LdUE)
legend("Energy","Potential","Kinetic")
hold off

figure
hold on
plot(ForceYSphere2/2-ForceYSphere1/2)
plot(ForceXSphere2/2-ForceXSphere1/2)
legend("Y Force (Positive is Attractive)","X Force")
hold off

figure
surf(X,Y,Ex1)
legend("Electric Field X")

figure
surf(X,Y,Ey1)
legend("Electric Field Y")

figure
surf(X,Y,drho1tot1)
legend("changed charge sphere 2")

figure
surf(X,Y,drho2tot1)
legend("Changed Charge sphere 1")

figure
surf(X,Y,Phi1)
legend("Potential")

figure
surf(X,Y,rho1)
legend("Charge Density")


mean=mean(ForceYSphere2/2-ForceYSphere1/2)