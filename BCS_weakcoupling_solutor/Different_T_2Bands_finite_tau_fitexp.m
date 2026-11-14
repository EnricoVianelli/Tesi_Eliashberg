%% FIT MODELLO A DUE BANDE (Zimmermann, tau finito) - MgB2, Pronin et al.
%  sigma_tot(w,T) = p*sigma_1(w;D1,g1,T) + (1-p)*sigma_2(w;D2,g2,T)
%  Parametri liberi:
%     p                    (globale, uguale per tutte le T)
%     D1(T), D2(T)  [meV]  (uno per ogni temperatura)
%     g1(T), g2(T)  [cm^-1](uno per ogni temperatura)
clear; clc; close all

%% ---------------- Impostazioni ----------------
Tc       = 32;                 % K
kB       = 0.08617;            % meV/K
cm2meV   = 0.123984;           % 1 cm^-1 = 0.123984 meV
T        = [3 15 20 27];       % K
nT       = numel(T);
sigma_dc = 0.65;               % sigma_n letto a occhio dai dati a 35 K
gs       = 100;                % gamma [cm^-1] (paper 150 +70/-50)

nStart   = 10;                 % numero di partenze casuali (per uscire da minimi locali) ~25 min
errRel   = 0.05;               % errore relativo assunto sui dati (pesi)
errFloor = 0.05;               % errore minimo assoluto = errFloor*max|dato| (per T e componente)

useGammaReg = true;            % penalizza gamma lontano dal valore normale (150 +70/-50 cm^-1)
gamma0      = 150;  gammaSig = 60;
wPen        = 10;              % peso delle penalità (ordine dei gap, monotonia, gamma)

%% ---------------- Dati sperimentali ----------------
%clean_wpd_csv('raw_sigma1_MgB2.csv', 'Real_cond_MgB2.csv', T)
%clean_wpd_csv('raw_sigma2_MgB2.csv', 'Img_cond_MgB2.csv', T)
Re = readtable('Real_cond_MgB2.csv');
Im = readtable('Img_cond_MgB2.csv');

D = struct();
D.T = T;  D.kB = kB;  D.cm2meV = cm2meV;  D.gs = gs;  D.nT = nT;
D.Tc = Tc; D.useGammaReg = useGammaReg; D.gamma0 = gamma0; D.gammaSig = gammaSig; D.wPen = wPen;
D.nuR = cell(nT,1); D.yR = D.nuR; D.eR = D.nuR;
D.nuI = cell(nT,1); D.yI = D.nuI; D.eI = D.nuI;
for i = 1:nT
    iR = abs(Re.T - T(i)) < 1e-6;
    iI = abs(Im.T - T(i)) < 1e-6;
    D.nuR{i} = Re.nu(iR);     D.yR{i} = Re.sigma(iR)/sigma_dc;   % normalizzati a sigma_n
    D.nuI{i} = Im.nu(iI);     D.yI{i} = Im.sigma(iI)/sigma_dc;
    D.eR{i}  = errRel*abs(D.yR{i}) + errFloor*max(abs(D.yR{i}));
    D.eI{i}  = errRel*abs(D.yI{i}) + errFloor*max(abs(D.yI{i}));
end

%% ---------------- Vettore parametri, guess iniziale e bounds ----------------
% x = [p, D1(1..nT), D2(1..nT), g1(1..nT)/gs, g2(1..nT)/gs]
D10 = 4.8;  D20 = 2.2;                          % guess a T=0
bcs = @(Tk) tanh(1.74*sqrt(max(Tc./Tk - 1, 0))); % andamento tipo BCS per il guess
x_init = [0.5, D10*bcs(T), D20*bcs(T), 150/gs*ones(1,nT), 150/gs*ones(1,nT)];

lb = [0,   0.05*ones(1,nT), 0.05*ones(1,nT),  30/gs*ones(1,nT),  30/gs*ones(1,nT)];
ub = [1,   8.00*ones(1,nT), 5.00*ones(1,nT), 600/gs*ones(1,nT), 600/gs*ones(1,nT)];

%% --- Fit (multi-start + lsqnonlin) ---
opts = optimoptions('lsqnonlin', 'Display','off', 'MaxIterations',400, ...
    'MaxFunctionEvaluations',1e5, 'FunctionTolerance',1e-10, 'StepTolerance',1e-10);

rng(1);
best.rn = Inf;
for k = 1:nStart
    x0 = x_init;
    if k > 1
        x0 = x_init .* (1 + 0.3*randn(size(x_init)));
        x0 = min(max(x0, lb), ub);
    end
    try
        [x, rn, res, ~, ~, ~, J] = lsqnonlin(@(x) residuals(x, D), x0, lb, ub, opts);
    catch ME
        warning('Start %d fallito: %s', k, ME.message);  continue
    end
    fprintf('start %2d: chi2 = %.4f\n', k, rn);
    if rn < best.rn
        best.rn = rn;  best.x = x;  best.res = res;  best.J = full(J);
    end
end

%% --- Risultati ---
[p, D1, D2, g1, g2] = unpack(best.x, D);
dof = numel(best.res) - numel(best.x);
% incertezze approssimate dalla Jacobiana
C   = pinv(best.J'*best.J) * best.rn/dof;
err = sqrt(abs(diag(C))).';
[ep, eD1, eD2, eg1, eg2] = unpack(err, D);  eg1 = eg1; eg2 = eg2;

fprintf('\n=== RISULTATO MIGLIORE: chi2 = %.3f, chi2/dof = %.3f ===\n', best.rn, best.rn/dof);
fprintf('p = %.3f +/- %.3f\n\n', p, ep);
fprintf('  T(K)   D1(meV)   D2(meV)   g1(cm^-1)   g2(cm^-1)\n');
for i = 1:nT
    fprintf('%6.1f  %7.3f   %7.3f   %8.1f    %8.1f\n', T(i), D1(i), D2(i), g1(i), g2(i));
end

%% --- Grafici ---
nu_fit = linspace(3, 30, 200);
colors = lines(nT);

figure; hold on                                   % Real
for i = 1:nT
    s = two_band(nu_fit, p, D1(i), D2(i), g1(i), g2(i), T(i), D)*sigma_dc;
    plot(nu_fit, real(s), '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Fit $T=%.0f$~K', T(i)));
    plot(D.nuR{i}, D.yR{i}*sigma_dc, 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('Data $T=%.0f$~K', T(i)));
end
grid on; set(gca,'FontSize',18,'XScale','log')
xlabel('Frequency (cm^{-1})'); ylabel('\sigma_1 (10^{-3} (\Omega/sq)^{-1})')
title('MgB$_2$ - Real part', 'Interpreter','latex')
legend('Interpreter','latex','Location','best', 'Fontsize', 18)

figure; hold on                                   % Imaginary
for i = 1:nT
    s = two_band(nu_fit, p, D1(i), D2(i), g1(i), g2(i), T(i), D)*sigma_dc;
    plot(nu_fit, imag(s), '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Fit $T=%.0f$~K', T(i)));
    plot(D.nuI{i}, D.yI{i}*sigma_dc, 'o', 'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
        'DisplayName', sprintf('Data $T=%.0f$~K', T(i)));
end
grid on; set(gca,'FontSize',18,'XScale','log')
xlabel('Frequency (cm^{-1})'); ylabel('\sigma_2 (10^{-3} (\Omega/sq)^{-1})')
title('MgB$_2$ - Imaginary part', 'Interpreter','latex')
legend('Interpreter','latex','Location','best', 'Fontsize', 18)

figure                                            % Parameters vs T
subplot(1,2,1); errorbar(T, D1, eD1, 'o-', 'LineWidth',1.5); hold on
errorbar(T, D2, eD2, 's-', 'LineWidth',1.5); grid on
xlabel('T (K)'); ylabel('\Delta (meV)'); legend('\Delta_1','\Delta_2'); set(gca,'FontSize',14)
subplot(1,2,2); errorbar(T, g1, eg1, 'o-', 'LineWidth',1.5); hold on
errorbar(T, g2, eg2, 's-', 'LineWidth',1.5); grid on
xlabel('T (K)'); ylabel('\gamma (cm^{-1})'); legend('\gamma_1','\gamma_2'); set(gca,'FontSize',14)

%% ================= FUNZIONI =================
function [p, D1, D2, g1, g2] = unpack(x, D)
nT = D.nT;
p  = x(1);
D1 = x(2:1+nT);
D2 = x(2+nT:1+2*nT);
g1 = x(2+2*nT:1+3*nT)*D.gs;      % cm^-1
g2 = x(2+3*nT:1+4*nT)*D.gs;      % cm^-1
end

function s = two_band(nu_cm, p, D1, D2, g1_cm, g2_cm, T, D)
% sigma/sigma_n complessa; nu in cm^-1, gap in meV, gamma in cm^-1
hw  = nu_cm(:).' * D.cm2meV;                 % meV
kT  = D.kB * T;                              % meV
s1  = zimmermann_finite_tau(hw, D1, kT, g1_cm*D.cm2meV);
s2  = zimmermann_finite_tau(hw, D2, kT, g2_cm*D.cm2meV);
s   = p*s1(:) + (1-p)*s2(:);
end

function r = residuals(x, D)
[p, D1, D2, g1, g2] = unpack(x, D);
r = [];
for i = 1:D.nT
    sR = two_band(D.nuR{i}, p, D1(i), D2(i), g1(i), g2(i), D.T(i), D);
    sI = two_band(D.nuI{i}, p, D1(i), D2(i), g1(i), g2(i), D.T(i), D);
    r  = [r; (real(sR) - D.yR{i}(:))./D.eR{i}(:); ...
             (imag(sI) - D.yI{i}(:))./D.eI{i}(:)];
end
% --- penalità (vincoli fisici) ---
pen = max(0, D2 - D1);                                    % D2 <= D1
pen = [pen, max(0, diff(D1)), max(0, diff(D2))];          % gap non crescenti con T (baias)
if D.useGammaReg
    pen = [pen, (g1 - D.gamma0)/D.gammaSig, (g2 - D.gamma0)/D.gammaSig];
end
r = [r; D.wPen*pen(:)];
end