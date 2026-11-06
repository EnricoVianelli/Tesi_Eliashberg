% Parameters
Tc     = 7.2;                           % K
kB     = 0.08617;                       % meV/K
c = physconst('LightSpeed') * 10^+2;    % c in cm
alpha = 3.95;                           % coupling 2*D0/KbTc in BCS
Delta0 = 0.5*alpha*kB*Tc;               % meV, BCS weak coupling

ratio_Tc = 0.8;
T = ratio_Tc*Tc;
kT = T*kB;
Delta  = Delta0*tanh(1.74*sqrt(Tc/T - 1));      % approx. BCS gap(T) for T near Tc

% Setting up the integration
f_Tz = linspace(0.01, 2, 400);                  % Freq in THz
hw = f_Tz * 4.13567;                            % Photon energies 1THz = 4.13567 meV
[s1, s2] = mattis_bardeen(hw, Delta, kT);       % Function to evaluate M-B integral conductivity

% Plotting
title_text = ['Conductivity of Pb ($T_c=$', num2str(Tc), '~K) at $T=$', num2str(T),' K'];
file_name = ['BCS_Pb_conductivity_conv_units_', num2str(T), 'K.pdf'];

figure
plot(10^12*f_Tz/c, s1, 'r', 10^12*f_Tz/c, s2, 'b', LineWidth=1.5); % Freq in cm^-1 = f(Hz) / c(in cm/s)
title(title_text, 'Interpreter', 'Latex')
ylim([-0.1 3])
xlabel('Frequency (cm^{-1})'); 
ylabel('\sigma/\sigma_n (\omega)');
legend('$\sigma_1/\sigma_n$','$\sigma_2/\sigma_n$', 'Interpreter', 'Latex', 'Fontsize', 16); grid on
set(gca, 'FontSize', 16) 
exportgraphics(gcf, file_name, 'ContentType', 'vector');