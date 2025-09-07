% Steady-State 2D Conduction with convective boundary condition on bottom
% Domain parameters
Lx = 1.0; % Length in x-direction
Ly = 1.0; % Length in y-direction
nx = 20; % Number of grid points in x-direction
ny = 20; % Number of grid points in y-direction
dx = Lx / (nx - 1); % Grid spacing in x-direction
dy = Ly / (ny - 1); % Grid spacing in y-direction
k = 1; % conductivity of the solid domain (W/mK)
% Convective boundary condition parameters (for bottom boundary)
h = 20; % Convective heat transfer coefficient(W/m2K)
T_inf = 300; % Ambient temperature for convective boundary
% Boundary conditions
T_top = 500; % Fixed temperature at the top boundary
T_left = 400; % Fixed temperature at the left boundary
T_right = 400; % Fixed temperature at the right boundary
% Initial guess for temperature (uniform 300 K)
T = 300 * ones(nx, ny);
% Convergence criteria
tolerance = 1e-5;
max_iterations = 10000;
error = 1;
iter = 0;
% Steady-state iterative solver for FDM (Gauss-Seidel)
while error > tolerance && iter < max_iterations
T_old = T;
% Update the interior points using central difference
for i = 2:nx-1
for j = 2:ny-1
T(i,j) = 0.25 * (T_old(i+1,j) + T_old(i-1,j) + T_old(i,j+1) + T_old(i,j-1));
end
end
% Applying boundary conditions
% Top boundary (Dirichlet boundary condition: T = 500 K)
T(nx, :) = T_top;
% Left boundary (Dirichlet boundary condition: T = 400 K)
T(:,1) = T_left;
% Right boundary (Dirichlet boundary condition: T = 400 K)
T(:,ny) = T_right;
% Bottom boundary (convective boundary condition)
T(1, :) = (T_old(2,:) + h*dy/k * T_inf) / (1 + h*dy/k);
% Corner Cases (mating points of the Walls)
T(nx,ny) = 0.5*(T_left + T_top);
T(nx,1) = T(nx,ny);
% Calculate error for convergence
error = max(max(abs(T - T_old)));
iter = iter + 1;
end
% Visualization of the temperature distribution
x = 0:dx:Lx;
y = 0:dy:Ly;
colormap(jet);
contourf(x, y, T,50);
colorbar;
title('Temperature distribution in 2d Domain of the wall');
xlabel('X grid points');
ylabel('Y grid points');