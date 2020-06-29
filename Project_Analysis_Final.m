%% Analysis 1

% Samantha Dalfen   260525758 
% Taylor Davies     260709971
% Sarah Ford        260705444
% Mathew Kfouri     260781493

%This code uses the values from Project_Final_309 to form charts and graphs

%% Question 1

contour(1:n,1:n,Phi,50)
colorbar

%% Question 2
Cp = zeros(n,1);
m = (20/50)*n;
j = 1;

%Calculating Cp
for i = 20/dx:21/dx
    Phi_x = (Phi(i+1,j) -  Phi(i-1,j))/(2*dx);
    Cp(m) = 2*(Phi_x)/U_inf;
    m = m + 1;
end
m = m-1;
figure(1)
plot(20:dx:21,Cp(20/50*n:m))
title('Surface Pressure Coefficient for Mach 0.85')
xlabel('X (m)')
ylabel('-Cp')

figure(2)
semilogy(1:k,error(1:k))
title('Error Convergence Plot for Mach 0.85')
xlabel('k')
ylabel('Log of Error')

%Calculating pressure

format long
P = zeros(n,n);
AA = ((Gamma-1)/2)*(M_inf^2);
for j = 1:n
    for i = 3:n-1
    u2 = ((Phi(i+1,j) -  Phi(i-1,j))/(2*dx))^2;
    P(j,i) = P_inf*(1+AA*(1-(u2/U_inf^2)))^(Gamma/(Gamma-1));
    end
    P(j,2) = P(j,3);
    P(j,1) = P(j,2);
    P(j,n) = P(j,n-1);
end

figure(3)
contour(1:n,1:n,P,25)
xlim([(19.5/50)*n (21.5/50)*n])
ylim([1 1/dx])
xticks([19.5/dx:dx*100:21.5/dx])
xticklabels({'19.5','19.75','20','20.25','20.5','20.75','21','21.25','21.5'})
yticks([1:40*dx:1/dx])
yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
title('Pressure Contour Mach 0.75');
xlabel('X (m)');
ylabel('Y (m)');

%% Question 3

plot(20:0.1:21,Cp01_08(20/50*500:(21/50)*500))
title('Surface Pressure Coefficient for Varying Grid Size at Mach 0.80')
xlabel('X (m)')
ylabel('-Cp')
hold on
plot(20:0.05:21,Cp005_08(20/50*1000:(21/50)*1000))
plot(20:0.025:21,Cp0025_08(20/50*2000:(21/50)*2000))

legend('Grid Size 0.1','Grid Size 0.05','Grid Size 0.025')

hold off

%% Question 4

plot(20:dx:21,Cp75(20/50*n:m))
title('Surface Pressure Coefficient for Varying Mach Numbers')
xlabel('X (m)')
ylabel('-Cp')
hold on
plot(20:dx:21,Cp77(20/50*n:m))
plot(20:dx:21,Cp79(20/50*n:m))
plot(20:dx:21,Cp81(20/50*n:m))
plot(20:dx:21,Cp83(20/50*n:m))
plot(20:dx:21,Cp85(20/50*n:m))
legend('Mach 0.75','Mach 0.77','Mach 0.79','Mach 0.81','Mach 0.83','Mach 0.85');

hold off