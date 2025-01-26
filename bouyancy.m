clear all
close all

video = false; % if want to save a movie, then let video = true
if(video)
    writerObj = VideoWriter('bouy_dmp0.5_10000_.mp4','MPEG-4');
    writerObj.FrameRate = 30;
    open(writerObj);
end

%set up object

g = 9.8; %m/s^2
mass = 900; %(kg)
n = 4;
sl = 1; %side length (meters)
init_vel = 0.000001; %initial velocity
volume = sl^3; %volume of cube
sa = 6*sl^2; %surface area of cube
umin = (10^-10);
dmp = 10; %damping constant



[A, N, n_pts, CoM, areas, masses] = bouy(n, mass, sl);

A_rel = A - CoM; 
% CoM = [0, 0, -0.4];
% A(:,3) = A(:,3) - 0.4;

A_rel = A - CoM; %relative positions to CoM

%water level will be at z = 0
rho_water = 1000; %density of water (kg/m^3)
rho_obj = mass/volume;

%initial velocity
% omega = [-1, 1, 0];
omega = zeros(1, 3);

vels = [ones(n_pts,1)*0, ones(n_pts, 1)*0, ones(n_pts, 1)];
vel_cm = [0*mean(vels(:,1)), 0*mean(vels(:,2)), 0*mean(vels(:,3))];

% vels = [zeros(n_pts,1), zeros(n_pts, 1), ones(n_pts, 1)*init_vel];
% vel_cm = [mean(vels(:,1)), mean(vels(:,2)), mean(vels(:,3))];
% vel_rel = vels - vel_cm;

%moment of inertia
I = zeros(3, 3);
L = zeros(1, 3);

for k = 1:n_pts
%     L = L + masses(k,:) * cross(A_rel(k,:), vels(k,:) - vel_cm);
    I = I + masses(k,:) * (norm(A_rel(k,:))^2*eye(3) - A_rel(k,:)' * A_rel(k,:));
end
% L = sum(masses(:,:) .* cross(A_rel(:,:), vels(:,:) - vel_cm));

%set up clock/time/dt
timesteps = 1000; %timesteps
time = 20; %seconds
dt = time/timesteps;
A_save = zeros(timesteps, 3);
t_save = zeros(timesteps, 1);

%set up graph
figure(1)
fig = figure(1);

fig.WindowStyle = 'normal';
% fig.WindowState = 'fullscreen';
plt_final = plot3(A(:, 1), A(:, 2), A(:, 3), 'o');
axis([-(sl*2), (sl*2), -(sl*2), (sl*2), -(sl*2), (sl*2)])
drawnow

Fsave = zeros(timesteps, 1);
tsave = zeros(timesteps, 1);
vel_cm_save = zeros(timesteps, 1);
Lsave = zeros(timesteps, 3);
L_norm_save = zeros(timesteps, 1);
x_cm_save = zeros(timesteps, 3);
omega_norm_save = zeros(timesteps, 1);

%set up main loop
for tick = 1:timesteps
    t = tick*dt;
    
    tick;
    CoM;
    
%     tic;
%     while toc < 5
%     end
    
%     omega = (I\L')';

    A = A_rel + CoM;
    vels = cross(repmat(omega,n_pts,1),A_rel)+vel_cm;
    
    %set up pressure forces
    pr = A(:, :); 

        %pr(pr(:, 3) >= 0) = 0; 

    for k = 1:n_pts
        if pr(k, 3) >= 0
            pr(k, :) = [0, 0, 0];
        end
    end

    Fpr = pr(:,3)*g*rho_water.*areas.*(N);
    
    %damping forces
    Fdamp = -dmp*(sl^2)*vels.*pr;

    %gravitational forces
    Fg = [zeros(n_pts, 1), zeros(n_pts, 1), -masses*g];
    
    equal = sum(Fpr(:,3)) + sum(Fg(:,3));
    
    %total forces
%     F = [Fpr(:,1) - Fdamp(:,1), Fpr(:,2) - Fdamp(:,2), (Fpr(:,3)) + Fg(:,3) - Fdamp(:,3)];

    %without damping
    F = [Fpr(:,1), Fpr(:,2), (Fpr(:,3)) + Fg(:,3)] - Fdamp;
    
    totalF = [sum(F(:, 1)), sum(F(:, 2)), sum(F(:,3))]

     %calculate total angular momentum
    for k = 1:n_pts
        L = L + dt*cross(A_rel(k,:), F(k,:));

    end
    omega =  (I \ L')';


    %apply rotation to points
    P = P_omega(omega);
    omega_norm = norm(omega);
    omega_cross = cross_omega(omega);
    R = P + cos(omega_norm*dt)*(eye(3) - P)...
        + sin(omega_norm*dt)*omega_cross;
    A_rel = (R*A_rel')';

    %apply rotation to normal vectors
    P = P_omega(omega);
    omega_norm = norm(omega);
    omega_cross = cross_omega(omega);
    R = P + cos(omega_norm*dt)*(eye(3) - P)...
        + sin(omega_norm*dt)*omega_cross;
    N = (R*N')';

    %update Center of Mass
    vel_cm = vel_cm + (dt/mass).*totalF;
    CoM = CoM + dt*vel_cm;

%     L = sum(masses(:,:).* cross(A_rel(:,:), vels(:,:) - vel_cm))
%     omega =  (I \ L')'

    %plot
    figure(1);
    plot3(A(:,1), A(:,2), A(:,3), 'o')


    axis([-(sl*2), (sl*2), -(sl*2), (sl*2), -(sl*2), (sl*2)])
%     axis([-20, 20, -20, 20, -20, 20])
    if(mod(tick, 100) == 0)
        drawnow
    end

    %info for graphs
    Fsave(tick, :) = norm(totalF);
    Lsave(tick,:) = L;
    L_norm_save(tick, :) = norm(L);
    x_cm_save(tick, :) = CoM;
    tsave(tick) = t;
    omega_norm_save(tick, :) = omega_norm;
end

figure(1)

plot(tsave, x_cm_save(:,3))
hold off
title('X_{cm}(z) vs Time')


figure(2)

tiledlayout(2,1) %2 x 1 graphs on one screen

nexttile
plot(tsave, Fsave)
title('||F_{total}||_2 vs Time')


nexttile
plot(tsave, L_norm_save)
title('||L||_2 vs Time')

function P = P_omega(omega)

% output: projection matrix that projects a vector onto omega

omega_norm = norm(omega,2);

P = zeros(3,3);

for i = 1:3
    for j = 1:3
        P(i,j) = omega(i)*omega(j)/(omega_norm^2);
    end
end

end

function P = cross_omega(omega)

% take the cross product with omega 

omega_norm = norm(omega,2);

P = zeros(3,3);

P(1,2) = -omega(3);
P(1,3) = omega(2);
P(2,1) = omega(3);
P(2,3) = -omega(1);
P(3,1) = -omega(2);
P(3,2) = omega(1);
P = P/omega_norm;

end
