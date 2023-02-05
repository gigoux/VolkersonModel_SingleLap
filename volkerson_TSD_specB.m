close all;
clear;
clc;

%%Parameters
E_o = 25e9; %outer adherend Young's modulus [Pa]
E_i = 25e9; %inner adherend Young's modulus [Pa]
t_o = 0.0613 * 0.0254; %outer adherend thickness [m]
t_i = 0.1274 * 0.0254; %inner adherend thickness [m]
G_a = 0.6e9; %adhesive shear modulus [Pa]
G_o = 2.6e9; %outer adherend shear modulus [Pa]
G_i = 2.6e9; %inner adherend shear modulus [Pa]
t_tot = 0.1917 * 0.0254; %total thickness (adherends + adhesive) [m]
c = (1.668 * 0.0254) / 2; %half of bonded length [m]
w_t = 1.0994 * 0.0254; %specimen width [m]
N_x = 5.2921e03 / w_t; %applied axial loading [N/m]
tau_ult = 38e6; %ultimate shear stress [Pa]


%%Volkerson Model Constants
t_a = t_tot - (t_o + t_i); %adhesive thickness [m]

lambda_sq = (G_a / t_a) * ( 1/(E_o*t_o) + 1/(E_i*t_i) );
lambda = sqrt(lambda_sq);

mu_sq = 1 + (G_a*t_o)/(4*t_a*G_o) + (G_a*t_i)/(4*t_a*G_i);
mu = sqrt(mu_sq);

C_0_OG = (G_a * N_x) / (E_i * t_a * t_i);
B_0_OG = N_x / (2*sinh(lambda*c));
A_0_OG = ( N_x - (2*C_0_OG/lambda_sq) ) / (2*cosh(lambda*c));

C_0_TSD = (G_a * N_x) / (E_i * t_a * t_i);
B_0_TSD = N_x / (2*sinh((lambda/mu)*c));
A_0_TSD = ( N_x - (2*C_0_TSD/lambda_sq) ) / (2*cosh((lambda/mu)*c));


%%Average Shear Stress
tau_ave = N_x / (2*c);


%%Calculate Adhesive Shear Stress Profile
x = linspace(-c, c, 1000);

tau_OG = lambda * ( (N_x/2 - C_0_OG/lambda_sq)*(sinh(lambda*x)/cosh(lambda*c)) + ...
    (N_x/2)*(cosh(lambda*x)/sinh(lambda*c)) );
tau_OG_max = max(tau_OG);
x_OG_max = x( find(tau_OG == tau_OG_max) );

tau_TSD = (lambda/mu)*( (N_x/2 - C_0_TSD/lambda_sq)*sinh((lambda/mu)*x)/cosh((lambda/mu)*c) + ...
        (N_x/2)*cosh((lambda/mu)*x)/sinh((lambda/mu)*c));
tau_TSD_max = max(tau_TSD);
x_TSD_max = x( find(tau_TSD == tau_TSD_max) );


%%Calculate Outer Adherend Axial Stress Profile
N_o_OG = A_0_OG*cosh(lambda*x) + B_0_OG*sinh(lambda*x) + C_0_OG/lambda_sq;
N_o_TSD = A_0_TSD*cosh((lambda/mu)*x) + B_0_TSD*sinh((lambda/mu)*x) + C_0_TSD/lambda_sq;


%%Calculate Axial Stress Failure
for ii = 1:length(x)
    test = 1;
end


%%Plot Shear Stress Profile vs Average Shear Stress
figure(1)
hold on;
grid on;
plot(x, tau_OG, 'r');
plot(x, tau_TSD, 'b');
plot(x, tau_ave*ones(length(x)), 'g');
plot(x, tau_ult*ones(length(x)), 'k');
xlabel('Bondline (from -c to c) [m]');
ylabel('Adhesive Shear Stress [Pa]');
title('Shear Stress Profile for Specimen B');
legend('Volkerson OG', 'Volkerson TSD', 'Average Shear Stress', ...
    'Ultimate Shear Stress', 'Location', 'northwest');


%%Plot Outer Adherend Axial Stress Profile
figure(2)
hold on;
grid on;
plot(x, N_o_OG, 'r')
plot(x, N_o_TSD, 'b');
xlabel('Bondline (from -c to c) [m]');
ylabel('Adherend Axial Stress [N/m]');
title('Outer Adherend Axial Stress Profile for Specimen B');
legend('Volkerson OG', 'Volkerson TSD', 'Location', 'best');