clear all
close all

%Initialization
h = figure;
axis tight manual %ensures that getframe() returns a consistent size
filename = 'testAnimated3d.gif';

L=2*pi;%Length
me=9.10938356*10^-31;%Mass of electron
DT=.05;%DeltaT
NT=1000;%Iterations
NG=3;%Number of Grid Points
N=500;%Number of Particles
WP=1;%Used to determine Q
QM=-1;
V0=0.2;
VXT=0.0;%Initial X velocity
VYT=0.0;%Initial Y velocity
VZT=0,.0;%Initial Z velocity
XP1=1;
YP1=1;
ZP1=1;
VX1=0.0;
VY1=0.0;
VZ1=0.0;
mode=1;
Q=WP^2/(QM*N/L);
rho_back=Q*N/L^3;
dx=L/NG;%Delta X
BoundIter=2;
acc=.9999

%color
c=linspace(1,10,N);

%initial positions
xp=linspace(0,L-L/N,N)';
xp=linspace(0,L/2,N)';
%xp=L*rand(N,1);
%yp=linspace(0,L-L/N,N)';
yp=L*rand(N,1)/2;
yp=L*rand(N,1);
%zp=zeros(N,1);
zp=L*rand(N,1)/2;
vxp=VXT*randn(N,1);
vyp=VYT*randn(N,1);
vzp=VZT*randn(N,1);



pm=[1:N]';
pm=1-2*mod(pm,2);
vxp=vxp;%+pm.*V0;
vyp=vyp;%+pm.*V0;

% Perturbation
vxp=vxp+VX1*sin(2*pi*xp/L*mode);
vyp=vyp+VY1*sin(2*pi*yp/L*mode);
vzp=vzp+VZ1*sin(2*pi*zp/L*mode);

%Building the Poisson Matrix
p=1:N;p=[p p;p p];
v=acc*ones(NG^2,1);
v1=ones(NG^2,1);
Poisson=-6*diag(v1)+diag(v(1:NG^2-1),1)+diag(v(1:NG^2-1),-1)+diag(v(1:NG^2-NG),-NG)+diag(v(1:NG^2-NG),NG);
Poisson=Poisson+diag(v(1:NG),NG^2-NG)+diag(v(1:NG),-NG^2+NG);

i=1;

for i=1:NG
    Poisson=Poisson+sparse(NG*i,NG*i-NG+1,acc,NG^2,NG^2)+sparse(NG*i-NG+1,NG*i,acc,NG^2,NG^2);
end
for i=1:(NG-1)
    Poisson=Poisson-sparse(NG*i,NG*i+1,acc,NG^2,NG^2)-sparse(NG*i+1,NG*i,acc,NG^2,NG^2);
end


i=1;
for i=1:NG
    j=2;
    if i==1
        temp=Poisson;
        for j=2:NG
            if j==2
                temp=cat(2,temp,diag(v(1:NG^2)));
            elseif j==NG
                temp=cat(2,temp,diag(v(1:NG^2)));
            else
                temp=cat(2,temp,zeros(NG^2,NG^2));
            end
            j=j+1;
        end
    else
        temp1=diag(v(1:NG^2));
        for j=2:NG
            if i==j
                temp1=cat(2,temp1,Poisson);
            elseif i==j+1
                temp1=cat(2,temp1,diag(v(1:NG^2)));
            elseif i==j-1
                temp1=cat(2,temp1,diag(v(1:NG^2)));
            elseif i==j+NG-1
                temp1=cat(2,temp1,diag(v(1:NG^2))); 
            elseif i==j-NG+1
                temp1=cat(2,temp1,diag(v(1:NG^2))); 
            else
                temp1=cat(2,temp1,zeros(NG^2,NG^2));
            end
            j=j+1;
        end
        temp=cat(1,temp,temp1);
    end
end
Poisson=temp;
clear temp;
clear temp1;
clear v;
clear acc;
clear v1;
%Poisson Built


%Iteration cycle
for it=1:NT
%plot xp,yp,zp
    plot3(xp,yp,zp,'o','Color','b','MarkerSize',5, 'MarkerFaceColor','r') 
    axis([0 L 0 L 0 L])
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
      
%update position
xp=xp+vxp*DT;
yp=yp+vyp*DT;
zp=zp+vzp*DT;

%apply boundary condition on the particle positions
n=1;
for n=1:BoundIter
out=(xp<0); 
xp(out)=xp(out)+L;
out=(xp>=L);
xp(out)=xp(out)-L;
out=(yp<0); 
yp(out)=yp(out)+L;
out=(yp>=L);
yp(out)=yp(out)-L;
out=(zp<0); 
zp(out)=zp(out)+L;
out=(zp>=L);
zp(out)=zp(out)-L;
end
%projection to grid points
g1=cat(2,ceil(xp/dx-.5),ceil(yp/dx-.5),ceil(zp/dx-.5));%put on grid points
g=[g1;g1+1];%put on grid points

%Weights
fraz1=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))')'*diag((1-abs(zp/dx-g1(1:N,3)+.5))'))';
fraz2=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))')'*diag((1-abs(zp/dx-g1(1:N,3)+.5))'))';
fraz3=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))')'*diag((1-abs(zp/dx-g1(1:N,3)+.5))'))';
fraz4=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))')'*diag((1-abs(zp/dx-g1(1:N,3)+.5))'))';
fraz5=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))')'*diag((abs(zp/dx-g1(1:N,3)+.5))'))';
fraz6=((abs(xp/dx-g1(1:N,1)+.5))'*diag((1-abs(yp/dx-g1(1:N,2)+.5))')'*diag((abs(zp/dx-g1(1:N,3)+.5))'))';
fraz7=((1-abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))')'*diag((abs(zp/dx-g1(1:N,3)+.5))'))';
fraz8=((abs(xp/dx-g1(1:N,1)+.5))'*diag((abs(yp/dx-g1(1:N,2)+.5))')'*diag((abs(zp/dx-g1(1:N,3)+.5))'))';

% apply bc on the projection
n=1;
for n=1:BoundIter
out=(g<1);
g(out)=g(out)+NG;
out=(g>NG);
g(out)=g(out)-NG;
n=n+1;
end


n=1;
for n=1:N
    if n==1
        weight=sparse([index(g(1,1),g(1,2),g(1,3)) index(g(1,1),g(1,2),g(N+1,3)) index(g(1,1),g(N+1,2),g(1,3)) index(g(1,1),g(N+1,2),g(N+1,3)) index(g(N+1,1),g(1,2),g(1,3)) index(g(N+1,1),g(1,2),g(N+1,3)) index(g(N+1,1),g(N+1,2),g(1,3)) index(g(N+1,1),g(N+1,2),g(N+1,3))],[1 1 1 1 1 1 1 1],[fraz1(1) fraz2(1) fraz3(1) fraz4(1) fraz5(1) fraz6(1) fraz7(1) fraz8(1)] ,NG^3,N);%matrix stating weights of each particle on each grid point
    else
        weight=weight+sparse([index(g(n,1),g(n,2),g(n,3)) index(g(n,1),g(n,2),g(N+n,3)) index(g(n,1),g(N+n,2),g(n,3)) index(g(n,1),g(N+n,2),g(N+n,3)) index(g(N+n,1),g(n,2),g(n,3)) index(g(N+n,1),g(n,2),g(N+n,3)) index(g(N+n,1),g(N+n,2),g(n,3)) index(g(N+n,1),g(N+n,2),g(N+n,3))],[n n n n n n n n],[fraz1(n) fraz2(n) fraz3(n) fraz4(n) fraz5(n) fraz6(n) fraz7(n) fraz8(n)] ,NG^3,N);%matrix stating weights of each particle on each grid point
    end
    n=n+1;
end


rho=-((Q/dx^2)*sum(weight,2))+rho_back;%rho of each grid point
n=1;
rho1=rho;

% computing fields
Phi=Poisson\(-rho1(1:NG^3)*dx^2);
n=1;
% for n=1:NG
%     if n==1
%         Phi1=Phi(1:NG);
%     else
%         Phi1=cat(2,Phi1,Phi((n-1)*NG+1:n*NG));
%     end
%     n=n+1;
% end


%Electric Field in x direciton
n=1;
for n=1:NG^3
    if n==1
        Ex=[(Phi(1+NG^2)-Phi(NG^3-NG^2+1))/(2*dx)];
    else
        in1=n+NG^2;
        in2=n-NG^2;
        if in1>NG^3
            in1=in1-NG^3;
        end
        if in2<1
            in2=in2+NG^3;
        end
        Ex=[Ex;(Phi(in1)-Phi(in2))/(2*dx)];
    end
    n=n+1;
end
n=1;
% for n=1:NG
%     if n==1
%         Ex1=Ex(1:NG);
%     else
%         Ex1=cat(2,Ex1,Ex((n-1)*NG+1:n*NG));
%     end
%     n=n+1;
% end


%Electric Field in Y Direction
n=1;
for n=1:NG^3
    if n==1
        Ey=[(Phi(1+NG)-Phi(NG^2-NG+1))/(2*dx)];
    else
        in1=n+NG;
        in2=n-NG;
        if mod(in2,NG^2)==0
            in2=in2+NG^2;
        elseif mod(in2,NG^2)>NG^2-NG
            in2=in2+NG^2;
        end 
        if mod(in1,NG^2)<NG+1
            if mod(in1,NG^2)~=0
                in1=in1-NG^2;
            end
        end
        Ey=[Ey;(Phi(in1)-Phi(in2))/(2*dx)];
    end
    n=n+1;
end
% n=1;
% for n=1:NG
%     if n==1
%         Ey1=Ey(1:NG);
%     else
%         Ey1=cat(2,Ey1,Ey((n-1)*NG+1:n*NG));
%     end
%     n=n+1;
% end


%Electric Field in Z Direction
n=1;
for n=1:NG^3
    if n==1
        Ez=[(Phi(2)-Phi(NG))/(2*dx)];
    else
        in1=n+1;
        in2=n-1;
        if mod(in1,NG)==1
            in1=in1-NG;
        end
        if mod(in2,NG)==0
            in2=in2+NG;
        end
        Ez=[Ez;(Phi(in1)-Phi(in2))/(2*dx)];
    end
    n=n+1;
end


% projection q->p and update of vp
n=1;
for n=1:N
    vxp(n)=vxp(n)+sum(weight(:,n)'*Ex)*QM*DT;
    n=n+1;
end
for n=1:N
    vyp(n)=vyp(n)+sum(weight(:,n)'*Ey)*QM*DT;
    n=n+1;
end
for n=1:N
    vzp(n)=vzp(n)+sum(weight(:,n)'*Ez)*QM*DT;
    n=n+1;
end
end

