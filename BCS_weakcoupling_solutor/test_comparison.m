kB    = 0.08617;                 % meV/K
alpha = [3.53 3.95 3.65];        % 2*Delta0/(kB*Tc) 
Tc    = [1.2 7.2 9.3];           % K
tr    = 0.8;                     % T/Tc fisso
c = physconst('LightSpeed') * 10^+2;
colors = {'b', 'r', 'g'}; 
materiali = {'Al', 'Pb', 'Nb'}; 

figure; 
hold on;

% Array per salvare gli handle grafici da mostrare in legenda
h_legenda = []; 
nomi_legenda = {};

for i = 1:numel(Tc)
    Delta0 = alpha(i)*kB*Tc(i)/2;                       % meV
    T      = tr*Tc(i);  
    Delta  = Delta0*tanh(1.74*sqrt(Tc(i)/T - 1));       % gap(T) approx.

    f_THz = linspace(0.01, 2, 400);                     % freq in THz
    hw    = f_THz*4.1357;                               % 1 THz = 4.1357 meV
    [s1, s2] = mattis_bardeen(hw, Delta, kB*T);         % Main function
    % Real Part -
    p1 = plot(10^12*f_THz/c, s1, '-', 'Color', colors{i}, 'LineWidth', 1);
    % Immaginary Part --
    plot(10^12*f_THz/c, s2, '--', 'Color', colors{i}, 'LineWidth', 1);

    % Salviamo solo il plot della parte reale per identificare il materiale nella legenda
    h_legenda(end+1) = p1; 
    nomi_legenda{end+1} = sprintf('%s (Tc = %g K)', materiali{i}, Tc(i));
end

% Linee dummy (invisibili sul grafico) per spiegare lo stile della linea
h_dummy_re = plot(NaN, NaN, '-k', 'LineWidth', 1.5);
h_dummy_im = plot(NaN, NaN, '--k', 'LineWidth', 1.5);

h_legenda = [h_legenda, h_dummy_re, h_dummy_im];
nomi_legenda = [nomi_legenda, {'Real - \sigma_1/\sigma_n'}, {'Imaginary - \sigma_2/\sigma_n'}];

% Configurazione estetica
ylim([-0.1 2]);
xlabel('Frequency (cm^{-1})'); 
ylabel('\sigma/\sigma_n (\omega)'); 
title(sprintf('Conductivity for T = %g T_c', tr));
legend(h_legenda, nomi_legenda, 'Location', 'best', 'FontSize', 16); 
grid on;
set(gca, 'FontSize', 16) 
exportgraphics(gcf, 'BCS_Comparison_0.8Tc.pdf', 'ContentType', 'vector'); % meglio di print
