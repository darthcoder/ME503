%% Solution of the Computer Assignment 1
%  Abdul Basit Ahmad
%  Roll no: 1221ME02


%% Parameters of the problem 
h = 100;
w = 100;
t = 1;
k = 1;                  % units in W/cmC



%% Mesh Generation
N = 20;                 % The total number of divisions.
dx = w/N;               %The grid seperation along x-axis
dy = h/N;               %The grid seperation along y-axis

% Arrays of x and y coordinates of the nodes.
% Initialization:
x = zeros(N+1,1);
y = zeros(N+1,1);

% Filling up the values:
for m = 1:N+1
    x(m) = (m-1)*dx;
    y(m) = (m-1)*dy;
end

% Arrays of x and y coordinates of the control volume centers.
% Initialization:
% There will be N+2 points because there will be...
% points on the edges of the domain too.
xc = zeros(N+2,1);
yc = zeros(N+2,1); 

% Filling up the values:
xc(1) = x(1);
xc(N+2) = x(N+1);
yc(1) = y(1);
yc(N+2) = y(N+1);

for m = 2:N+1
    yc(m) = (y(m)+y(m-1))/2;
    xc(m) = (x(m)+x(m-1))/2;
end

%=========================================================================
%% ANALYTICAL SOLUTION
%=========================================================================
%% Now the boundary condition is specified.
T_anal = zeros(N+1);


% At the y=h boundary:
for m = 1:N+1
    T_anal(m,N+1) = 100*sin(pi*(m-1)*dx/w);
end

% Temperature is zero at other boundaries. 

%% Now the analytical solution is specified
for m = 2:N
    for n = 2:N
        T_anal(m,n) = 100*sinh(pi*(n-1)*dy/w)*sin(pi*(m-1)*dx/w)...
            /sinh(pi*h/w);
    end
end

%=========================================================================



%=========================================================================
%% NUMERICAL SOLUTION
%=========================================================================


iter = 0;
Told = zeros(N+2);
T = zeros(N+2);
eps = 1e-6;

%% Boundary condition:
% The only non zero b.c. is at y=h boundary.
Told(:,N+2) = 100*sin(pi*xc/w);

for iter = 0:1000


    %% Guess value of the solution:
    T = Told;

    %% For vertex control volumes:

    T(2,2) = (1/6)*(2*Told(1,2)+2*Told(2,1)+Told(3,2)+Told(2,3));
    T(N+1,2) = (1/6)*(2*Told(N+2,2)+2*Told(N+1,1)+Told(N+1,3)+Told(N,2));
    T(N+1,N+1) = (1/6)*(2*Told(N+2,N+1)+2*Told(N+1,N+2)+Told(N+1,N)...
              +Told(N,N+1));
    T(2,N+1) = (1/6)*(2*Told(1,N+1)+2*Told(2,N+2)+Told(3,N+1)+Told(2,N));

    %% For the other control volumes on the edges:

    for m = 3:N
        T(m,2) = (1/5)*(2*Told(m,1)+Told(m+1,2)+Told(m-1,2)+Told(m,3));
        T(2,m) = (1/5)*(2*Told(1,m)+Told(2,m-1)+Told(2,m+1)+Told(3,m));
        T(N+1,m) = (1/5)*(2*Told(N+2,m)+Told(N,m)+Told(N+1,m-1)+...
                    Told(N+1,m+1));
        T(m,N+1) = (1/5)*(2*Told(m,N+2)+Told(m-1,N+1)+Told(m+1,N+1)...
                    +Told(m,N));
    end


    %% For the internal control volumes:
    for m = 3:N
        for n = 3:N
            T(m,n) = (1/4)*(Told(m+1,n)+Told(m-1,n)+Told(m,n+1)...
                        +Told(m,n-1));
        end
    end
    
    difft = max(abs(T-Told),-1);
    Told = T;
    
    if difft < eps
        break
    else continue;
    end
end
%=========================================================================


%=========================================================================
%% Comparison of Numerical and Analytical Solution.
%=========================================================================

%% Transforming the Numerical solution to the physical nodes.
Tnum = zeros(N+1);
for m = 2:N
    for n = 2:N
        Tnum(m,n) = (1/4)*(T(m,n)+T(m+1,n)+T(m,n+1)+T(m+1,n+1));
    end
end

for m = 1:N+1
    Tnum(m,N+1) = 100*sin(pi*(m-1)*dx/w);
end

%% Computing the error 
error = abs(Tnum-T_anal);
%=========================================================================


%=========================================================================
%% Plots
%=========================================================================
surf(T_anal); figure(gcf)       % Analytical Solution
surf(Tnum); figure(gcf)         % Numerical Solution
surf(error); figure(gcf)        % Error
%=========================================================================
