% Parameters
Tc     = 8.3;                           % K
kB     = 0.08617;                       % meV/K
alpha = 4.1;                           % coupling 2*D0/KbTc
Delta0 = 0.5*alpha*kB*Tc;               % meV, BCS weak coupling
c = physconst('LightSpeed') * 10^+2;    % c in cm
%ratio_Tc = [0.001 0.3 0.6 0.8 0.99];
%T = ratio_Tc * Tc;
T = [4.5 7];                          % Passing T in K directly
kT = T * kB;

f_THz = linspace(0.01, 2, 400);   % 1THz ~ 35cm^-1
hw = f_THz * 4.13567;
x = 10^12*f_THz/c;
nT = numel(T);
nW = numel(hw);

sigma1 = zeros(nT, nW);
sigma2 = zeros(nT, nW);

% Iterate over T for BCS prediction
for i = 1:nT
    Delta = Delta0*tanh(1.74*sqrt(Tc/T(i) - 1));
    [s1, s2] = mattis_bardeen(hw, Delta, kT(i));
    sigma1(i, :) = s1; sigma2(i, :) = s2;
end

% experimental data
sigma_dc = 2.4; % units sigma [10^5 ohm^-1 cm^-1] at T = 9 K
raw = readtable('nb_digitized_raw.csv');
s1_data = raw(raw.panel==1,:);  s2_data = raw(raw.panel==2,:);
%s1_data.sigma = s1_data.sigma/sigma_dc;           % normalizzation
%s2_data.sigma = s2_data.sigma/sigma_dc;

% colors
colors = lines(nT); 

% Plotting real part --------------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :); 
    % 1. BCS
    plot(x, sigma1(i, :)*sigma_dc, 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));

    % 2. Data
    idx_exp = abs(s1_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s1_data.nu(idx_exp), s1_data.sigma(idx_exp), '--o', ...
        'Displayname',sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Real part Conductivity of Nb ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 3])
%xlim([0.0, 100])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_1/\sigma_n (\omega)');
ylabel('\sigma_1(10^5\Omega^{-1}cm^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_Nb_real_cond_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');


% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :);

    % 1. BCS
    plot(x, sigma2(i, :)*sigma_dc, 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));

    % 2. Data
    idx_exp = abs(s2_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s2_data.nu(idx_exp), s2_data.sigma(idx_exp), '--o', ...
        'DisplayName', sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Imaginary part Conductivity of Nb ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 8.5])
%xlim([0.0, 100])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_2/\sigma_n (\omega)');
ylabel('\sigma_2(10^5\Omega^{-1}cm^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_Nb_imm_cond_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');

% ----------------------- Gaps -------------------------------

temps = unique(s1_data.T);   
temps = temps(temps < Tc); % no T > Tc because the minimum here has no sense
w_min = zeros(size(temps));
s_min  = zeros(size(temps));
range_T = linspace(0.01, Tc*0.999, 100);
Delta_BCS = 2*Delta0*tanh(1.74*sqrt(Tc./range_T - 1))./(4.13567*c*10^-12); % gap with BCS approx in cm^-1

for i = 1:numel(temps)
    idx = abs(s1_data.T - temps(i)) < 1e-6;
    nu  = s1_data.nu(idx);
    sg  = s1_data.sigma(idx);

    [nu, order] = sort(nu);          % order by freq
    sg = sg(order);
    [s_min(i), k] = min(sg);         % min value and index
    nu_min(i) = nu(k);               % freq of the min

    % parabolic interpolation to find better minimum
    if k > 1 && k < numel(nu)
        p = polyfit(nu(k-1:k+1), sg(k-1:k+1), 2);
        v = -p(2)/(2*p(1));
        if p(1) > 0 && v >= nu(k-1) && v <= nu(k+1)
            nu_min(i) = v;                    
        end
    end  
end

sigma_nu = 2*ones(numel(temps),1); % error from the paper
% Extract minimum from MB model
nu_model = zeros(size(temps));
for i = 1:numel(temps)
    Delta = Delta0*tanh(1.74*sqrt(Tc/temps(i) - 1));
    [s1m, ~] = mattis_bardeen(hw, Delta, kB*temps(i));
    [~, km] = min(s1m);
    nu_model(i) = hw(km)/(4.13567*c*10^-12); % conv in cm^-1
end

% Plot gap
figure
hold on
errorbar(temps, nu_min, sigma_nu, 'o', 'Displayname', 'Data', 'Color', 'r', 'MarkerFaceColor', 'r') % data
plot(temps, nu_model, 'v', 'Displayname', 'Mattis-Bardeen', 'Color', 'b', 'MarkerFaceColor', 'b') % MB results
plot(range_T, Delta_BCS, '-', 'Displayname', 'BCS approx', 'Color', 'k') % BCS curve approx
grid on
title('Gap of Nb')
xlim([0 8.4])
ylim([0 26])
xlabel('T(K)'); 
ylabel('2 \Delta(T) (cm^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_Gap_MB_vs_T_fitexp.pdf', 'ContentType', 'vector');