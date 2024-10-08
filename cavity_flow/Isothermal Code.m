%%%Isothermal code

clc;
clear;
close all

%%%Parameters

Re= 100;                %Reynolds Number
Pr = .7;                %Prantel Number
Gr = 10^4;              %Grashove Number
Pe = Re*Pr;             %Peccklet Number
Ri = Gr/Re^2;           %Richardson Number
r = 1;                  %Density

Nx = 99;               %Number of nodes in x-direction
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

%%%Plots

cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
figure(1); plot(x,v(:,round(Ny/2)));
title('Centerline Y-direction velocity');
xlabel('X'); ylabel('V'); axis('square'); xlim([0 L]); grid on

figure(2); plot(y,u(round(Ny/2),:));
title('Centerline X-direction velocity');
xlabel('Y'); ylabel('U'); axis('square'); xlim([0 L]); grid on

N = 1000; xstart = max(x)*rand(N,1); ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(x,y);
figure(4); h=streamline(x,y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 L 0 L]); set(h,'color','k')