clear;close all;clc;
no_of_omega_values = 5;
omega_min = -20;
% Q is symmetric about origin
omega_max = - omega_min;
% Create omega values possible (mode set Q) and its |.|
omega_values = linspace(omega_min, omega_max, no_of_omega_values );

omega_pdf_1 = 1/no_of_omega_values * ones(1, no_of_omega_values);
omega_pdf_2 = [0.5 0.45 0.05 0.0 0.0];

figure();
stem(omega_values, omega_pdf_1,'LineWidth',4,'MarkerSize',20);
hold on;
stem(omega_values, omega_pdf_2,'LineWidth',4,'MarkerSize',20);
xlabel('\omega')
ylabel('\psi_{\omega}');
xticks(omega_values)
axis([omega_min-2 omega_max+2 0 0.6]);
set(gca,'FontSize',30);
grid on;
box on;
legend('Figure 3.a','Figure 3.b');
