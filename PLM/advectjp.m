% -------------------------------------------------%
%                                                  %
% MATLAB code for Price model                      %
%   u_t + (FLUX(1+q2-2q2*u))_x = x*q2*(u^2)_xx     %
%    where u is the CDF of prices                  %
%                                                  %
% Written by Dongwook Lee and JP Rabanal           %
% -------------------------------------------------%

format long;

%clc;
close all;
clear all;


% Take inputs from users =======================
N=input('Grid resolution N (e.g. N=8) = ');
disp('');
xa=0;
xb=1;
dx = (xb-xa)/N;
ngc=2;
ibeg=ngc+1;
iend=N+ngc;

% Setup your discrete cells
x(ngc+1:N+ngc)=linspace(xa+0.5*dx,xb-0.5*dx,N);
for i=1:ngc;
    x(ngc-i+1) = x(ngc+1) - i*dx;
    x(i+N+ngc) = x(N+ngc) + i*dx;
end

% Initialize arrays to zero
u(1:N+2*ngc)=0.;
Flux(N+2*ngc+1)=0.;

% Array for u^n+1
uNew = u;
uOld = u;

q2 = input('q2 = ');
diff_factor = input('Input diff_factor (e.g., 1 (Recommended!), 10, 100, etc.) = ');
disp('');

%Set IC 
disp('[1] Uniform, [2] NE, [3] NE q2=.1');
InCond = input('Initial condition (choose any)= ');
fsol = zeros(40000,iend+2);

if (InCond == 1)
   u=x;
elseif (InCond == 2)
    for i=ibeg:iend
        if (x(i) > xa) & (x(i) < (1-q2)/(1+q2))
            u(i) = 0;
        else
            u(i) = 1-(1-q2)*(1-x(i))/(2*q2*x(i));
        end
    end
elseif (InCond == 3)
    q3 = .10;
    for i=ibeg:iend
        if (x(i) > xa) & (x(i) < (1-q3)/(1+q3))
            u(i) = 0;
        else
            u(i) = 1-(1-q2)*(1-x(i))/(2*q3*x(i));
        end
    end
end

BCtype = 3;
PDEtype = 3;
TimeAdvance = 1;

% Call BC for initial data
u=applyBC(u,ibeg,iend,ngc,BCtype);

% Final simulation time t=tmax
tmax=input('tmax= ');

% CFL as input
cfl = input('CFL (0<CFL<=1 ');

% Compute initial dt
lambda = abs(1+q2-4*q2*u(i));

dt = ComputeDt(lambda,dx,cfl);
dt_adv = dt;
kappa = q2*max(abs(x));
u_max = max(abs(u));
dt_diff = cfl*0.5*dx^2/(kappa*u_max);
dt = min(dt_adv,dt_diff);

% Slow start for maximum stability
dt = 0.1*dt;

% Maximum range for plotting
maxu=max(u);
minu=min(u);
marginDist = 0.1*(maxu-minu);

%hold on;
plot(x(ibeg:iend),u(ibeg:iend),'k','LineWidth',1.2);
axis([0 1 minu-marginDist maxu+marginDist]);

% Choose your numerical method
disp('[1]FOG (1st order), [2]PLM (2nd order)');
methodType = input('Method Type (choose any)= ');

% Choose slope limiter
if (methodType == 2)
    disp('[1]minmod   [2]mc           [3]vanLeer');
    disp('[4]centered [5]upwind (a>0) [6]downwind (a<0)');
    slopeLimiter = input('Choose slope limiter (1, 2, 3 are recommended) = ');
end

% Set your clock to zero
t=0;

% Exact solution =================================
 for i=ibeg:iend;
     if (x(i) > xa) & (x(i) < (1-q2)/(1+q2));
         uExact(i) = 0;
     else
         uExact(i) = 1-(1-q2)*(1-x(i))/(2*q2*x(i));
     end
 end

% Apply BC for exact solution
uExact=applyBC(uExact,ibeg,iend,ngc,BCtype);
%hold on
plot(x(ibeg:iend),uExact(ibeg:iend),'k','LineWidth',1.2);


nStep = 1;
% Begin the evolution loop
fsol(nStep,:) = u; 
while t<=tmax
    
    % Reset the max shock speed (or advection velocity) to be zero
    % at each time step
    sMax = 0;
    
    for i=ibeg-1:iend+1;
        % Reconstruction
        if (methodType == 1)
            % No TVD slope for first-ordr Godunov
            delta(i) = 0;
        elseif (methodType == 2)
            
            % TVD slope limiter for second-order PLM
            if (slopeLimiter == 1)
                delta(i) = minmod(u(i)-u(i-1),u(i+1)-u(i));
                
            elseif (slopeLimiter == 2)
                delta(i) = mc(u(i)-u(i-1),u(i+1)-u(i));
                
            elseif (slopeLimiter == 3)
                delta(i) = vanLeer(u(i)-u(i-1),u(i+1)-u(i));
                
            % non TVD slope limiter for second-order PLM
            elseif (slopeLimiter == 4)
                delta(i) = centered(u(i)-u(i-1),u(i+1)-u(i));
                
            elseif (slopeLimiter == 5)
                delta(i) = upwind(u(i)-u(i-1),u(i+1)-u(i));
                
            elseif (slopeLimiter == 6)
                delta(i) = downwind(u(i)-u(i-1),u(i+1)-u(i));                
            end
        end
    end
 
    for i=ibeg:iend+1;
      % Characteristic tracing at each i-1/2 interface
      s=1+q2-4*q2*0.5*(u(i-1)+u(i));
      sMax=max(sMax,abs(s));
      
      % Compute CFL number
      Ca = s*dt/dx;
        
      % Evaluate left and right Riemann states at each i-1/2
      uLeft = u(i-1)+0.5*delta(i-1)*(1.-Ca);
      uRght = u(i)  -0.5*delta(i)  *(1.+Ca);
        
      % Godunov flux to solve local Riemann problems at each i-1/2
      Flux(i) = GodunovFlux(s,uLeft,uRght,PDEtype,q2);
        
    end

    % Update
    for i=ibeg:iend;
        uNew(i) = u(i) - dt/dx*(Flux(i+1) - Flux(i));
    end

    % Apply BC
    uNew=applyBC(uNew,ibeg,iend,ngc,BCtype);
    %uOld = u;
    dt_old = dt;
    u = uNew;
    
    % Add diffusion for the Price model
        % Update
     if (TimeAdvance == 1)
         for i=ibeg:iend;
             uNew(i) = u(i) + q2*x(i)*dt/dx^2*(u(i+1)^2-2*u(i)^2+u(i-1)^2);
         end
     else %Crank-Nicolson using extrapolation
          if (nStep > 1)
              for i=1:N+2*ngc;
                  uNew_pred(i) = u(i) + 0.1*dt/dt_old*(u(i) - uOld(i));
              end
              for i=ibeg:iend;
                  uNew(i) = u(i) + 0.5*q2*x(i)*dt/dx^2 ...
                        *(  u(i+1)^2-2*u(i)^2+u(i-1)^2 ...
                        + uNew_pred(i+1)^2-2*uNew_pred(i)^2+uNew_pred(i-1)^2);
              end
          else %first time step, use explicit
              for i=ibeg:iend;
                  uNew(i) = u(i) + q2*x(i)*dt/dx^2*(u(i+1)^2-2*u(i)^2+u(i-1)^2);
              end
          end
      end
    
    % Apply BC
    uNew=applyBC(uNew,ibeg,iend,ngc,BCtype);
    uOld = u;
    u = uNew;
    dt_old = dt;
    
    % update time
    t  = t + dt;
    nStep = nStep + 1;
    
    % Call CFL for a new dt for the next update
    lambda = sMax;
    fsol(nStep,:) = u; 

    %lambda = abs(max(u));
    dt = ComputeDt(lambda,dx,cfl);
   
    dt_adv = dt;
    if (TimeAdvance == 1) % Do this if explicit
        kappa = q2*max(abs(x));
        u_max = max(abs(u));
        dt_diff = cfl*0.5*dx^2/(kappa*u_max);
        dt = min(dt_adv,diff_factor*dt_diff);
    end
    
%     % Save the updated solution
%     u = uNew;

    % Plot
    %hold on
    
    % UNCOMMENT THESE LINES IF WANT TO GENERATE MOVIE
    % (1) NUMERICAL SOLUTION
    %plot(x(ibeg:iend),uNew(ibeg:iend),'ro:','MarkerSize',5,'LineWidth',1.2);
    
    % (2) EXACT SOLUTION
    %plot(x(ibeg:iend),uExact(ibeg:iend),'k','LineWidth',1.2);
    
    % (3) AXIS
    %axis([0 1 minu-marginDist maxu+marginDist]);

    %  (3A) OR, PLOT GUARDCELLS ALSO
    %  plot(x,uNew,'ro:','MarkerSize',5,'LineWidth',1.2);
    %  axis([x(1) x(end) minu-marginDist maxu+marginDist]);
    
    % (4) MOVIE 
    %pause(0.1)
    %hold off

end
% Plot the final result
% Computed solution
tplot = round(nStep/10,0);
hold on;
plot(x(ibeg:iend),u(ibeg:iend),'ro:','MarkerSize',5,'LineWidth',1.2);
plot(x(ibeg:iend),fsol(tplot,ibeg:iend),'b--','LineWidth',1.2);
plot(x(ibeg:iend),fsol(1,ibeg:iend),'k:','LineWidth',1.2);
xlabel('x')
ylabel('F(x)')
% Exact solution
if (PDEtype == 1)
    plot(x(ibeg:iend),uExact(ibeg:iend),'k','LineWidth',1.2);
elseif (PDEtype == 2) & (ICtype == 1)
    % Let's not plot for sin wave solution to Burgers
    % as there is no exact solution in general.
    % Although we can get an asymptotic solution that looks like /|/
    plot(x(ibeg:iend),uExact(ibeg:iend),'k','LineWidth',1.2);
elseif (PDEtype == 3)
    plot(x(ibeg:iend),uExact(ibeg:iend),'k','LineWidth',1.2);
end
    
% Plot with proper scale
axis([0 1 minu-marginDist maxu+marginDist]);




