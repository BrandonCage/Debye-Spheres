load('Temp50test7')
% 
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated2d.gif';
% 

itcont=it;

% Main computational cycle
for it=itcont:NT
% % update xp
% scatter(xp,yp,20,c,'filled') 
% 
% hold on
% sx1=[sxc1-sr1:2*sr1/100:sxc1+sr1];
% sy1=(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
% nsy1=-(sr1.^2-(sx1-sxc1).^2).^(1/2)+syc1;
% sx2=[sxc2-sr2:2*sr2/100:sxc2+sr2];
% sy2=real((sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);
% nsy2=real(-(sr2.^2-(sx2-sxc2).^2).^(1/2)+syc2);

% plot(sx1,real(sy1),'r')
% plot(sx1,real(nsy1),'r')
% plot(sx2,real(sy2),'r')
% plot(sx2,real(nsy2),'r')
% axis([L/NG L-L/NG L/NG L-L/NG])
% 
%     drawnow 
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if it == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
% hold off
xp=xp+vxp*DT;
yp=yp+vyp*DT;
syc1=syc1+vsy1*DT;
syc2=syc2+vsy2*DT;
sxc1=sxc1+vsx1*DT;
sxc2=sxc2+vsx2*DT;



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

randmat=randi([0 1],Nnew,1);%Ion or Electron?

randvtot=(-2*log((1-rand(Nnew,1)))./(((mp-me)*randmat+me)/me/k/Tempin)).^.5;%Maxwellian
randang=2*pi*rand(Nnew,1);%Random direction
vxp=[vxp;randvtot.*cos(randang)];
vyp=[vyp;randvtot.*sin(randang)];


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
    %dimension is a new particle. Adding all allements on    
    %each row should
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

Volt1(it)=PhiC1;
Volt2(it)=PhiC2;


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
            x(n)=(-Phi(in2))/(2*dx);%Electric field in x, zero out edge
        elseif in2<1
            Ex(n)=(Phi(in1))/(2*dx);%Electric field in x, zero out edge
        else
        Ex(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in x
        end
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
            Ey(n)=(-Phi(in2))/(2*dx);%Electric field in y, zere out edge
        elseif mod(in2,NG)==0
            Ey(n)=(Phi(in1))/(2*dx);%Electric field in y, zero out edge
        else
        Ey(n)=(Phi(in1)-Phi(in2))/(2*dx);%Electric field in y
        end
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
    
 
 ForceXSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ex;
 ForceYSphere1(it)=dx^2*(drho1tot+esphere1*rho_sphere1)'*Ey;
 ForceXSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ex;
 ForceYSphere2(it)=dx^2*(drho2tot+esphere2*rho_sphere2)'*Ey;
 
 vsx1=vsx1;%-ForceXSphere1*DT/ms1;%Change sphere 1 x velocity
 vsy1=vsy1;%-ForceYSphere1*DT/ms1;%Change sphere 1 y velocity
 vsx2=vsx2;%-ForceXSphere2*DT/ms2;%Change sphere 2 x velocity
 vsy2=vsy2;%-ForceYSphere2*DT/ms2;%Change sphere 2 y velocity
 
 for n=1:NG
    if n==1
        drho2tot1=drho2tot(1:NG)+esphere2*rho_sphere2(1:NG);
    else
        drho2tot1=cat(2,drho2tot1,drho2tot((n-1)*NG+1:n*NG)+esphere2*rho_sphere2((n-1)*NG+1:n*NG));
    end
 end

 for n=1:NG
    if n==1
        drho1tot1=drho1tot(1:NG)+esphere1*rho_sphere1(1:NG);
    else
        drho1tot1=cat(2,drho1tot1,drho1tot((n-1)*NG+1:n*NG)+esphere1*rho_sphere1((n-1)*NG+1:n*NG));
    end
 end
 
 for n=1:NG
    if n==1
        drho2tot2=drho2tot(1:NG);
    else
        drho2tot2=cat(2,drho2tot2,drho2tot((n-1)*NG+1:n*NG));
    end
 end

 for n=1:NG
    if n==1
        drho1tot2=drho1tot(1:NG);
    else
        drho1tot2=cat(2,drho1tot2,drho1tot((n-1)*NG+1:n*NG));
    end
 end
 
  for n=1:NG
    if n==1
        drho2tot3=esphere2*rho_sphere2(1:NG);
    else
        drho2tot3=cat(2,drho2tot3,esphere2*rho_sphere2((n-1)*NG+1:n*NG));
    end
 end

 for n=1:NG
    if n==1
        drho1tot3=esphere1*rho_sphere1(1:NG);
    else
        drho1tot3=cat(2,drho1tot3,esphere1*rho_sphere1((n-1)*NG+1:n*NG));
    end
 end
 
 Charge1(it)=esphere1;
 Charge2(it)=esphere2;
 
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

PhiC2*510999
end




% 
% meanyforce=mean(ForceYSphere2/2-ForceYSphere1/2)
% meanxforce=mean(ForceXSphere2/2-ForceXSphere1/2)
% 
% sum(rho_sphere2(1:780))
% sum(rho_sphere2(781:1600))