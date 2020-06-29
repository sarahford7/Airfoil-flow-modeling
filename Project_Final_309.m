%% Numerical Final Project

% Samantha Dalfen   260525758 
% Taylor Davies     260709971
% Sarah Ford        260705444
% Mathew Kfouri     260781493

%This code uses the Gauss-Sidel method to solve for Phi in two-dimensional
%transsonic flow.
%% Set up Constants

clc 
clear
Gamma = 1.4;
R = 287.058;
T_inf = 293;
P_inf = 100;
tc = 0.08;

%Choose constants 
dx = 0.1;     %grid size
dy = 0.1;     %grid size
M_inf = 0.75;

n = 50/dx;
U_inf = M_inf*sqrt(Gamma*R*T_inf);

%% Pre computations
% This helps speed up computations by having constants based on variables
% that we choose already calculated, instead of calculating in the loop.

Ac1 = (1-M_inf^2); %First half of A equation
if Ac1 > 0 %Finding U
    U_set = 0;
else 
    U_set = 1;
end
Ac2 = -(Gamma +1)*(M_inf^2/U_inf); %Second half of A equation
% b is constant
b = 1/(dy^2);
c = b;
% dyxy = t/c(-4x+82) derivative of y
k = 1;
E = 0;
tic

%% Set up Phi Matrix

Phi = zeros(n,n); 
error = ones(2500,1)';

%% Gauss-Seidel

PhiMax = Phi(1,1);
while (error(k) > 10^-4) && (k < 20000)
    E = 0; %The value representing the error for the current iteration
for j = 2:(n-1)
    
    %The following section gets computed for i=2 as Phi(i-2,j) does not
    %exist and can not be computed using gauss-sidel
    i = 2;
    A_Prev = Ac1;
    A_Current = A_Prev;
    U_Prev = U_set;
    U = U_Prev;
    
     for i = 3:(n-1) 
         
        %Allows us to avoid storing every iteration
        PhiLast = Phi(i,j); 
        A_Prev = A_Current;
        U_Prev = U;
        
        A_Current = Ac1 + Ac2*(Phi(i+1,j)-Phi(i-1,j))/(2*dx);
            if A_Current < 0 
                U = 1;
            else 
                U = 0;
            end
            
            %Calculating coefficients
            a = ((U_Prev*A_Prev)-2*(1-U)*A_Current-2)*b;
            d = ((1-U)*A_Current-2*U_Prev*A_Prev)*b;
            e = (1-U)*A_Current*b;
            g = (U_Prev*A_Prev)*b;
            %b and c are constant and defined in the pre-computation
            %section
            
            Phi(i,j) = (-c*Phi(i,j-1) -g*Phi(i-2,j) -d*Phi(i-1,j) - e*Phi(i+1,j) - b*Phi(i,j+1))/a;
            
            E = max(E,abs(Phi(i,j)-PhiLast)); %Compare to avoid storing every error
    end
end
    %Updating boundary condition
    j = 1;
    for i = 1:n
        PhiLast = Phi(i,j);
        if  (i*dx <= 21 && i*dx >= 20) %BC on airfoil
            Phi(i,1) = Phi(i,2) - U_inf*dy*(tc*(-4*i*dx+82));
        else %BC off airfoil
            Phi(i,1) = Phi(i,2);
        end
        E = max(E,abs(Phi(i,j)-PhiLast));
    end
error(k+1) = E; %Stores the maximum error of the iteration
k = k+1;
toc
end
k
toc