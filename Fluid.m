%========================================================================
%  A very simple Navier-Stokes solver 2-D enclosed cavity. 
%The viscosity is taken to be a constant and a forward in time, 
%centered in space discretization is used. 
%The density is advected by a simple upwind scheme.
%========================================================================
%clear all
clc
%domain size and physical variables
Lx=1.0;Ly=1.0;gx=0.0;gy=00.0; rho1=1.0; rho2=1.0; m0=0.001; rro=rho1;
unorth=1;usouth=0;veast=0;vwest=0;time=0.0;bl=0;lt=0;rtt=0;rb=0;

% Numerical variables
nx=100;ny=100;dt=(1/(nx*unorth)*0.9);nstep=200000; maxit=300;maxError=0.0001;beta=1.2;
% Zero various arrys
u=zeros(nx+1,ny+2); v=zeros(nx+2,ny+1); p=zeros(nx+2,ny+2);
ut=zeros(nx+1,ny+2); vt=zeros(nx+2,ny+1); tmp1=zeros(nx+2,ny+2);
uu=zeros(nx+1,ny+1); vv=zeros(nx+1,ny+1); tmp2=zeros(nx+2,ny+2);
% Set the grid
dx=Lx/nx;dy=Ly/ny;
for i=1:nx+2; x(i)=dx*(i-1.5);end; for j=1:ny+2; y(j)=dy*(j-1.5);end;
% Set density
r=zeros(nx+2,ny+2)+rho1;

%================== START TIME LOOP======================================
for is=1:nstep,is
    uold=u;
    
% tangential velocity at boundaries
u(1:nx+1,1)=2*usouth-u(1:nx+1,2);u(1:nx+1,ny+2)=2*unorth-u(1:nx+1,ny+1);
v(1,1:ny+1)=2*vwest-v(2,1:ny+1);v(nx+2,1:ny+1)=2*veast-v(nx+1,1:ny+1);

for i=2:nx,for j=2:ny+1 % TEMPORARY u-velocity
ut(i,j)=u(i,j)+dt*(-0.25*(((u(i+1,j)+u(i,j))^2-(u(i,j)+ ...
u(i-1,j))^2)/dx+((u(i,j+1)+u(i,j))*(v(i+1,j)+ ...
v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/dy)+ ...
m0/(0.5*(r(i+1,j)+r(i,j)))*( ...
(u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ ...
 (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2 )+gx );
end,end

for i=2:nx+1,for j=2:ny % TEMPORARY v-velocity

vt(i,j)=v(i,j)+dt*(-0.25*(((u(i,j+1)+u(i,j))*(v(i+1,j)+ ...
v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/dx+ ...
((v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)/dy)+ ...
m0/(0.5*(r(i,j+1)+r(i,j)))*( ...
(v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+ ...
(v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2 )+gy );
end,end
%========================================================================
% Compute source term and the coefficient for p(i,j)
rt=r;

for i=2:nx+1,for j=2:ny+1
tmp1(i,j)= (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy );
tmp2(i,j)=1.0/( (1./dx)*( 1./(dx*(rt(i+1,j)+rt(i,j)))+...
1./(dx*(rt(i-1,j)+rt(i,j))) )+...
(1./dy)*(1./(dy*(rt(i,j+1)+rt(i,j)))+...
1./(dy*(rt(i,j-1)+rt(i,j))) ) );
end,end
for it=1:maxit % SOLVE FOR PRESSURE
oldArray=p; p(2:nx+1,1)=p(2:nx+1,2);
        p(2:nx+1,ny+2)=p(2:nx+1,ny+1);
        p(1,2:ny+1)=p(2,2:ny+1);
        p(nx+2,2:ny+1)=p(nx+1,2:ny+1);
        
for i=2:nx+1,for j=2:ny+1
p(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(...
(1./dx)*( p(i+1,j)/(dx*(rt(i+1,j)+rt(i,j)))+...
p(i-1,j)/(dx*(rt(i-1,j)+rt(i,j))) )+...
(1./dy)*( p(i,j+1)/(dy*(rt(i,j+1)+rt(i,j)))+...
p(i,j-1)/(dy*(rt(i,j-1)+rt(i,j))) ) - tmp1(i,j));
%                                     end
end,end
if max(max(abs(oldArray-p))) <maxError, break,end
end
for i=2:nx,for j=2:ny+1 % CORRECT THE u-velocity
u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
end,end
for i=2:nx+1,for j=2:ny % CORRECT THE v-velocity
v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
end,end
%=======ADVECT DENSITY using centered difference plus diffusion ==========
ro=r;

%========================================================================
time=time+dt % plot the results
uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
for i=1:nx+1,xh(i)=dx*(i-1);end; for j=1:ny+1,yh(j)=dy*(j-1);end

hold off,contour(x,y,flipud(rot90(p))),axis equal,axis([0 Lx 0 Ly]);
hold on;quiver(xh,yh,flipud(rot90(uu)),flipud(rot90(vv)),'r');
%  hold on;streamline(xh,yh,flipud(rot90(uu)),flipud(rot90(vv)),[],[]);
pause(0.00001);
% hold on, loglog(is,max(max(abs(uold-u))),'b*');drawnow();
if max(max(abs(uold-u)))<maxError break;
end

end
