%%% Solution of the Assignment 3.
%%% Blows up to infinity.


%%Parameters

Tin = 273+120;
Tout1 = 273+10;
Tout2 = 273+15;
% Tin = 120;
% Tout1 = 10;
% Tout2 = 15;

hin = 910;
hout1 = 8;
hout2 = 280;

% rin = 0.15;
% rout = 0.30;
rin = 15;
rout = 30;

epsilon = 1e-6;

%Constants not given in the problem taken from net, could be wrong.
% alpha = 5.6e4;
% k = 0.025;
alpha = 1;
k = 1;

%%Grid Generation

%Let there be M divisions in r direction and N divisions in theta 
%direction. 

M = 100;
N = 100;

%Total M+1 points in r and N+1 points in theta directions
deltar = (rout-rin)/M;
delta_theta = 2*pi/N;

%For deltat
deltat = 3.375e-3;
%%Initial time Temperature distribution, specified in problem.
Told = (273+30)*ones(M+1,N+1);
% Told = (30)*ones(M+1,N+1);
T = Told;
maxiter = 100000;
%The maximum number of iterations the code will run for.
%Some expressions to simplify the equations
a = (2*alpha*deltat)/(deltar^2);
b = (2*alpha*deltat*hin)/(k*deltar);
c = (alpha*hin*deltat)/(k*rin);
d = (2*alpha*deltat)/((rin^2)*(delta_theta^2));
e = (alpha*deltat/(deltar^2));
f = (alpha*deltat/(rin^2)*(delta_theta^2));

bup = (2*alpha*deltat*hout1)/(k*deltar);
bdown = (2*alpha*deltat*hout2)/(k*deltar);
cup = (alpha*hout1*deltat)/(k*rin);
cdown = (alpha*hout2*deltat)/(k*rin);
 
for t = 1:maxiter
    %Insert boundary condition
    m = 1;
    for n = 2:N
        T(m,n) = (1-a+b-c-d)*Told(m,n)+e*Told(m+1,n)+f*Told(m,n+1)...
                 +f*Told(m,n-1)+(b-c)*Tin;    
    end
    T(1,1) = (1-a+b-c-d)*Told(1,1)+e*Told(2,1)+f*Told(1,2)...
                 +f*Told(1,N+1)+(b-c)*Tin;
    T(1,N+1) = (1-a+b-c-d)*Told(1,N+1)+e*Told(2,N+1)+f*Told(1,1)...
                 +f*Told(1,N)+(b-c)*Tin;    
    
    m = M+1;
    for n = 2:N/2
        T(m,n) = a*Told(m-1,n)+(1-e-bup-cup+f)*Told(m,n)+f*Told(m,n+1)...
                 +f*Told(m,n-1)+(bup+cup)*Tout1;
    end
    T(M+1,1) = a*Told(M,1)+(1-e-bup-cup+f)*Told(M+1,1)+f*Told(M+1,2)...
                 +f*Told(M+1,N+1)+(bup+cup)*Tout1;
             
    for n = N/2+1:N
        T(m,n) = a*Told(m-1,n)+(1-e-bdown-cdown+f)*Told(m,n)...
                 +f*Told(m,n+1)+f*Told(m,n-1)+(bdown+cdown)*Tout2;
    end
    T(M+1,N+1) = a*Told(M,N+1)+(1-e-bdown-cdown+f)*Told(M+1,N+1)...
                 +f*Told(M+1,1)+f*Told(M+1,N)+(bdown+cdown)*Tout2;
    %Write governing equation
    for m = 2:M
        for n = 2:N
            r = rin+((rout-rin)/M)*m;
            g = (alpha*deltat)/(2*r*deltar);
            h = (alpha*deltat)/((r^2)*(delta_theta^2));
            T(m,n) = (e+g)*Told(m+1,n)+(1-e-h)*Told(m,n)...
                     +(e-g)*Told(m-1,n)+h*Told(m,n+1)+h*Told(m,n-1);
        end
    end
    
    %Check for convergence
    difftemp = max(-1,(norm((T-Told),'fro'))^2);
    if difftemp < epsilon
        break;
    else Told = T;
        continue
    end
end