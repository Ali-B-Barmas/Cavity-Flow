
clc;
clear;
close all

%%%Parameters

Re= 100;                %Reynolds Number
Pr = .7;                %Prantel Number
Gr = 10^5;              %Grashove Number
Pe = Re*Pr;             %Peccklet Number
Ri = Gr/Re^2;           %Richardson Number
r = 1;                  %Density

Nx = 44;               %Number of nodes in x-direction
L = 1;                  %Dimensionless length(Domain size)
U = 1;                  %Velocity on upper wall

deltt = 0.001;          %Time increments
max_iteration = 50000;  %Maximum iteration         
max_error = 1e-5;       %Maximum error

%%% Setting matrixes and requirements for code

Ny = Nx;
d=L/(Nx-1);
x = 0:d:L;
y = 0:d:L;
w = zeros(Nx,Ny);       %Vorticity matrix
wp = w;                 %Previous time iteration values of vorticity
S = zeros(Nx,Ny);       %Stream function
u = zeros(Nx,Ny);       %X-direction velocity
u(2:Nx-1,Ny) = U;       %Top surface velocity
v = zeros(Nx,Ny);       %Y-direction velocity
T = zeros(Nx,Ny);       %Temperature matrix
T(Nx,:) = 1;            %Boundary condition on right wall

%%%Unsteady solving for w and T

for iter = 1:max_iteration

%%%Boundary conditions

    %%%for vorticity

    w(:,Ny) = -2*S(:,Ny-1)/(d^2) - U*2/d;       % Top
    w(:,1) = -2*S(:,2) /(d^2);                  % Bottom
    w(1,:) = -2*S(2,:) /(d^2);                  % Left
    w(Nx,:) = -2*S(Nx-1,:)/(d^2);               % Right

    %%%for temperature

    T(2:Nx-1,1) = T(2:Nx-1,2);                  %Bottom
    T(2:Nx-1,Ny) = T(2:Nx-1,Ny-1);              %Top

%%%Calculating vorticity

    wp = w;
    Tp=T;

    for i = 2:Nx-1
        for j = 2:Ny-1
            w(i,j) = wp(i,j)+(-1*max(u(i,j),0).*((wp(i,j)-wp(i-1,j))/d)+...
          max(-1*u(i,j),0).*((wp(i+1,j)-wp(i,j))/d)-...
          max(v(i,j),0).*((wp(i,j)-wp(i,j-1))/d)+...
          max(-1*v(i,j),0).*((wp(i,j+1)-wp(i,j))/d)+...
          1/Re*(wp(i+1,j)+wp(i-1,j)-4*wp(i,j)+wp(i,j+1)+wp(i,j-1))/(d^2))*deltt...
          + deltt*Ri*(T(i+1,j)-T(i-1,j)/(2*d));

            T(i,j)=Tp(i,j) + (-1*u(i,j).* (Tp(i+1,j)-Tp(i-1,j))/(2*d)+...
          (-v(i,j)).*(Tp(i,j+1)-Tp(i,j-1))/(2*d)+...
          1/Pe*(Tp(i+1,j)+Tp(i-1,j)-4*Tp(i,j)+Tp(i,j+1)+Tp(i,j-1))/(d^2))*deltt;
        end
    end

%%%Calculating stream function

    for i = 2:Nx-1
        for j = 2:Ny-1
            S(i,j) = (w(i,j)*d^2 + S(i+1,j) + S(i,j+1) + S(i,j-1) + S(i-1,j))/4;
        end
    end

    %%%Calculating velocities

    for i = 2:Nx-1
        for j =2:Ny-1
            u(i,j) = (S(i,j+1)-S(i,j-1))/(2*d);
            v(i,j) = (-S(i+1,j)+S(i-1,j))/(2*d);
        end
    end

%%%Convergence check

    if iter > 10
        error = max(max(w - wp));
        if error < max_error
            break;
        end
    end
end



%%%Calculating pressure after reaching to steady state

p = ones(Nx,Ny);
err_presure = 1;
t = 0;

while err_presure > 0.00001
    
    po = p;
    for i = 2:Nx-1
        for j = 2:Ny-1
            p(i,j) = r/(4*d^2)*(2*((S(i+1,j)-2*S(i,j)+S(i-1,j))*(S(i,j+1)-2*S(i,j)+S(i,j-1))-...
                1/16*(S(i+1,j+1)-S(i-1,j+1)-S(i+1,j-1)+S(i-1,j-1))^2)-...
                Ri*d^3/2*(T(i,j+1)-T(i,j-1))+d^2/r*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)));
        end
    end

    if t > 10
        err_presure = max(max(p - po));
    end
    t = t+1;
end 

%%%Dimensionless heat transfer rate on left wall

for j = 1:Ny
    q(j,1) = -(T(1,j)-T(2,j))/d;
end

%%%Dimensionless friction force on left wall

for j = 1:Ny
    c(j,1) = 1/Re*v(2,j)/d;
end

%%%Plots

cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
figure(1); contourf(x,y,u',23,'LineColor','none');
title('U-contour'); xlabel('X'); ylabel('Y')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');

figure(2); contourf(x,y,v',23,'LineColor','none');
title('V-contour'); xlabel('X'); ylabel('Y')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');

figure(3); contourf(x,y,T',23,'LineColor','none');
title('Teta-contour'); xlabel('X'); ylabel('Y')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');

N = 1000; xstart = max(x)*rand(N,1); ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(x,y);
figure(4); h=streamline(x,y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 L 0 L]); set(h,'color','k')

figure(5); contourf(x,y,w',23,'LineColor','none');
title('Vortivity-contour'); xlabel('X'); ylabel('Y')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');

figure(6); plot(y,u(round(Nx/2),:));
title('X-direction velocity(U) on X=0.5');
xlabel('Y'); ylabel('U'); axis('square'); xlim([0 L]); grid on

figure(7); plot(x,u(:,round(Ny/2)));
title('X-direction velocity(U) on Y=0.5');
xlabel('X'); ylabel('U'); axis('square'); xlim([0 L]); grid on

figure(8); plot(y,v(round(Nx/2),:));
title('Y-direction velocity(V) on X=0.5');
xlabel('Y'); ylabel('V'); axis('square'); xlim([0 L]); grid on

figure(9); plot(x,v(:,round(Ny/2)));
title('Y-direction velocity(V) on Y=0.5');
xlabel('X'); ylabel('V'); axis('square'); xlim([0 L]); grid on

figure(10); plot(y,T(round(Nx/2),:));
title('Temperature distribution on X=0.5');
xlabel('Y'); ylabel('Teta'); axis('square'); xlim([0 L]); grid on

figure(11); plot(x,T(:,round(Ny/2)));
title('Temperature distribution on Y=0.5');
xlabel('X'); ylabel('Teta'); axis('square'); xlim([0 L]); grid on

figure(12); contourf(x,y,p',23,'LineColor','none');
title('Pressure-contour'); xlabel('X'); ylabel('Y')
axis('equal',[0 L 0 L]); colormap(cm); colorbar('westoutside');

figure(13); plot(q,y);
title('Dimensionless heat transfer rate on left wall');
xlabel('Q'); ylabel('Y'); axis('square'); ylim([0 L]); grid on

figure(14); plot(c,y);
title('Dimensionless friction force on left wall');
xlabel('C'); ylabel('Y'); axis('square'); ylim([0 L]); grid on