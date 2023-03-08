close all
clear all
clc
%Midterm project _lagging coefficent 3D

%% Grid Parameters
nx = 15;
ny = 15;
nz = 5;
Ngrids = nx * ny * nz;
dx = 30;
dy = 30;
dz = 30;

% Fluid & Rock Properties
kx = 70; % mD
ky = 70; % mD
kz = 7;  % mD
phio = 0.25; % Porosity, Fraction
C_o = 1e-5;
C_r = 3e-6;
C_t = C_o + C_r;
% Fluid Properties at reference conditions
rho_o = 45;
Bo_o = 1;
% Transmissibility Geometric part
Tgx = dy*dz/dx;
Tgy = dz*dx/dy;
Tgz = dx*dy/dz;
V = dx * dy * dz; % Volume of cell
%
D_o = V*phio*C_t/Bo_o;
%Initial & Boundary condition
Pi = 3000;
Pwf = 2700;
%Ti = 170;
P0 = Pi;

% Well Specifications
Well_i = 1;
Well_j = 15;
Well_k = 3;
s = 3; % Skin Factor
rw = 0.35;
re = 0.208*dx;
a = log(re/rw)+s;
b = 2*pi*dz;
Well_order = (Well_j-1)*nx + Well_i+nx*ny*(Well_k-1);
WI=b/a;
% Time specifications
ti = 0;
dt = 0.1;
tf = 120;
tprob = [30 60 90 120]; % Solution time

% Array Allocation
x = zeros(nx,ny,nz);
y = zeros(nx,ny,nz);
z = zeros(nx,ny,nz);

%Transmissibility - flow part
Tf_E = zeros (nx,ny,nz);
Tf_W = zeros (nx,ny,nz);
Tf_S = zeros (nx,ny,nz);
Tf_N = zeros (nx,ny,nz);
Tf_B = zeros (nx,ny,nz);
Tf_T = zeros (nx,ny,nz);

%Total transmissibility
T_E = zeros (nx,ny,nz);
T_W = zeros (nx,ny,nz);
T_N = zeros (nx,ny,nz);
T_S = zeros (nx,ny,nz);
T_T = zeros (nx,ny,nz);
T_B = zeros (nx,ny,nz);
T_C = zeros (nx,ny,nz);

%Pressure and flow rate matrices
P = zeros(nx,ny,nz);              
X = ones(nx,ny,nz);               
P_save = ones(Ngrids,(tf-ti)/dt); 
Q_save = zeros((tf-ti)/dt,1);    

%Lambda = K/B-Mu
lambda_x = zeros (nx,ny,nz);
lambda_y = zeros (nx,ny,nz);
lambda_z = zeros (nx,ny,nz);


%Fluid properties
rho = zeros (nx,ny,nz);
miu = zeros (nx,ny,nz);
gamma = zeros (nx,ny,nz);
Bo = zeros(nx,ny,nz);

%Fluid properties function
rho_function = @(Xfunction) rho_o*exp(C_o*(Xfunction - P0));
Bo_function = @(Xfunction)Bo_o*exp(-C_o*(Xfunction-P0));
miu_function = @(Xfunction)6.00961538e-9*(Xfunction).^2-9.1332176e-5*(Xfunction)+1.21;
%gamma = rho_function*32.14;

%Matrix equation temrs
Tmat = zeros (Ngrids, Ngrids);
Dmat = zeros (Ngrids, Ngrids);
Dvec = zeros (Ngrids, 1);
GvecT = zeros (Ngrids, 1);
GvecB = zeros (Ngrids, 1);
Gvec = zeros (Ngrids, 1);
Pvec = zeros (Ngrids, 1);
Pnew = zeros (Ngrids, 1);
Avec = zeros(Ngrids,Ngrids);
Bvec = zeros(Ngrids,1);
Qvec = zeros(Ngrids,1);
Qmat = zeros(Ngrids,Ngrids);

%% Computation section with time step loop
% Initial condition for pressure
for k=1:nz
    P(:,:,k) = (3000+0.312*(k-3)*dz);
    z(:,:,k)=(k-1)*30;
end

for k=1:nz   
    for j=1:ny
        for i=1:nx 
            Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
            gg = Grid_order;     
            X(gg,1)=P(i,j,k);
        end
    end
end

% setting counter term
count = 1;

% Time step loop
for t=ti:dt:tf
    % calculation of fluid properties which change based on pressure
    rho = rho_function(P); 
    miu = miu_function(P);
    Bo = Bo_function(P);
    gamma = rho_function(P)*32.174;

    %calculation of lambda term
    lambda_x = kx/(Bo.*miu);
    lambda_y = ky/(Bo.*miu);
    lambda_z = kz/(Bo.*miu);

    % populate Tf matrices with harmonic averages of lambda

    for i = 2:nx
        Tf_W(i,:,:) = 2*lambda_x(i-1,:,:).*lambda_x(i,:,:)./...
            (lambda_x(i-1,:,:)+ lambda_x(i,:,:));
    end

    for i = 1:nx-1
        Tf_E(i,:,:) = 2*lambda_x(i+1,:,:).*lambda_x(i,:,:)./...
            (lambda_x(i+1,:,:)+ lambda_x(i,:,:));
    end


    for j = 1:ny-1
        Tf_N(:,j,:) = 2*lambda_y(:,j+1,:).*lambda_y(:,j,:)./...
            (lambda_y(:,j+1,:)+ lambda_y(:,j,:));
    end

    for j = 2:ny
        Tf_S(:,j,:) = 2*lambda_y(:,j-1,:).*lambda_y(:,j,:)./...
            (lambda_y(:,j-1,:)+ lambda_y(:,j,:));
    end

    for k = 2:nz
        Tf_T(:,:,k) = 2*lambda_z(:,:,k-1).*lambda_z(:,:,k)./...
            (lambda_z(:,:,k-1)+ lambda_z(:,:,k));
    end

    for k = 1:nz-1
        Tf_B(:,:,k) = 2*lambda_z(:,:,k+1).*lambda_z(:,:,k)./...
            (lambda_z(:,:,k+1)+ lambda_z(:,:,k));
    end

    %populate  transmissibility matrices 
    T_B(:,:,:) =  0.00633*Tf_B.*Tgz;
    T_T(:,:,:) =  0.00633*Tf_T.*Tgz;
    T_E(:,:,:) =  0.00633*Tf_E.*Tgx;
    T_W(:,:,:) =  0.00633*Tf_W.*Tgz;
    T_N(:,:,:) =  0.00633*Tf_N.*Tgy;
    T_S(:,:,:) =  0.00633*Tf_S.*Tgz;

    %Populating Tmat
    %central transmissibility term, which is the summation of all
    %transmissibilities
    for k=1:nz
        for j=1:ny
            for i=1:nx
                Grid_order = (j-1)*nx+i+nx*ny*(k-1);
                g = Grid_order;
                Tmat(g,g) = -(T_E(i,j,k)+T_W(i,j,k)+T_N(i,j,k)+...
                    T_S(i,j,k)+T_T(i,j,k)+T_B(i,j,k));
            end
        end
    end

    %next, the off diagonal terms east and west
    for k=1:nz
        for j=1:ny
            for i=2:nx 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h-1)=T_W(i,j,k);           
            end
        end
    end
    
    for k=1:nz
        for j=1:ny
            for i=1:(nx-1) 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h+1)=T_E(i,j,k);           
            end
        end
    end

    %next, the off diagonal terms north and south
    for k=1:nz
        for j=1:(ny-1)
            for i=1:nx 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h+nx)=T_N(i,j,k);          
            end
        end
    end

    for k=1:nz
        for j=2:ny
            for i=1:nx 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h-nx)=T_S(i,j,k);          
            end
        end
    end

    %next, the off diagonal terms top and bottom   
    for k=2:nz
        for j=1:ny
            for i=1:nx 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h-nx*ny)=T_T(i,j,k);
            end
        end
    end

    for k=1:(nz-1)
        for j=1:ny
            for i=1:nx 
                Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
                h=Grid_order;  
                Tmat(h,h+nx*ny)=T_B(i,j,k); 
            end
        end
    end

    %Now, the Gravity terms (T*rho*g) will be computed and populated

    for k=2:(nz-1)
        for j=1:ny
            for i=1:nx    
                Grid_order = (j-1)*nx+i+(k-1)*(nx*ny);
                Gvec(Grid_order) = 0.5/144*((rho(i,j,k)+rho(i,j,k-1))*...
                    T_T(i,j,k)*(z(i,j,k-1)-z(i,j,k))+(rho(i,j,k)+...
                    rho(i,j,k+1))*T_B(i,j,k)*(z(i,j,k+1)-z(i,j,k)));
            end
        end
    end

    for j=1:ny
        for i=1:nx   
            k=1;
            Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
            Gvec(Grid_order)=0.5/144*((rho(i,j,k)+rho(i,j,k+1))*...
                T_B(i,j,k)*(z(i,j,k+1)-z(i,j,k)));
        end
    end

    for j=1:ny
        for i=1:nx   
            k=5;
            Grid_order=(j-1)*nx+i+(k-1)*(nx*ny);
            Gvec(Grid_order)=0.5/144*((rho(i,j,k)+rho(i,j,k-1))...
                *T_T(i,j,k)*(z(i,j,k-1)-z(i,j,k)));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Now, the Dmatrix and D vector will be calculated
    for k=1:nz
        for j=1:(ny)
            for i=1:(nx)
                Grid_order = (j-1)*nx+i+nx*ny*(k-1);
                g = Grid_order;
                Dmat(g,g) = (D_o/dt);   % diagonal terms; off-diagonal terms are all zero
                Dvec(g,1) =  Dmat(g,g)*X(g,1); % the D vector term
            end 
        end
    end

    prod = WI*0.006333*kx/(Bo(1,15,3)*miu(1,15,3));
    Qmat(Well_order,Well_order)= prod; % all other terms are zero
    Qvec(Well_order) = -prod*Pwf; % Qvec just multiply with pwf because 
    % the future time step pressure is being calculated

    Amat =  Tmat-Dmat - Qmat;
    Bvec = -Dvec+ Qvec + Gvec;
    X = Amat\Bvec;

    P_save(:,count) = X;
    Q_save(count,1)=(X(Well_order,1)-Pwf)*WI*0.006333*300/...
        (Bo(1,15,3)*miu(1,15,3));
    count = count + 1;
    P = reshape(X,15,15,5);

    %then go to next time step

    x = 1:15;
    y = 1:15;
    [XX,YY] = meshgrid(x,y);
    zz = -6900-2*dz:dz:-6900+2*dz;
    if sum(t==tprob)==1
        figure;
        for i = 1:length(zz)
            subplot(length(zz),1,i)
            contourf(XX,YY,P(:,:,i))
            xlabel('x Dimension')
            ylabel('y Dimension')
            s = "Pressure at z = " + zz(i);
            title(s)
            c = colorbar;
            c.Label.String = 'Pressure (psi)';
        end
        s = "Pressure at t = " + t + " days";
        title(s)
    end
end
figure;
plot((0:0.1:120),Q_save)
xlabel('time (day)')
ylabel('Oil production rate (SCF/d)')
title('Oil production rate vs. Time')