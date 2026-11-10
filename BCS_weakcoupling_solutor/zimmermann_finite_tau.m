function [s1, s2] = zimmermann_finite_tau(hw, Delta, kT, gamma)

% gamma - scattering rate 1/tau

sigma = zeros(size(hw));
opts = {'RelTol', 1e-6, 'AbsTol', 1e-9}; 

% Regolarizzazione numerica infinitesima per evitare i NaN sulle singolarità delle radici
eps_reg = 1e-12; 

% Funzione ausiliaria per la radice quadrata complessa "safe" vicino alle singolarità
safe_sqrt = @(x) sqrt(x + 1i*eps_reg);

% Fermi-Dirac f(E)
fermi = @(E) 0.5 * (1 - tanh(E / (2 * max(kT, 1e-20)))); 

if kT > 0
    th = @(E) tanh(E / (2 * kT));
else
    th = @(E) ones(size(E)); 
end

for k = 1:numel(hw)
    w = hw(k);
    E_max = Delta + 40 * max(kT, 1e-3 * Delta);
    % Coherence Factor
    C = @(E) (E.^2 + Delta^2 + w.*E) ./ (safe_sqrt(E.^2 - Delta^2) .* safe_sqrt((E + w).^2 - Delta^2));
    % Energy Density
    N = @(E) safe_sqrt(E.^2 - Delta^2) + safe_sqrt((E + w).^2 - Delta^2);

    % --- Quasi-particle scattering regime ---
    Integrand_QP = @(E) (fermi(E) - fermi(E + w)) .* ( ...
        (C(E) - 1) ./ (N(E) - w - 1i * gamma) - (C(E) + 1) ./ (N(E) + w + 1i * gamma) );
    I_A = integral(Integrand_QP, Delta, E_max, opts{:});

    % --- Couple Breaking regime ---
    if w > 2 * Delta
        Integrand_CB = @(E) th(E + w) .* ((C(E) - 1) ./ (N(E) - w - 1i * gamma));

        I_B = integral(Integrand_CB, Delta - w, -Delta, opts{:});
        I = I_A + I_B;
    else
        I = I_A;
    end

    sigma(k) = gamma / (gamma - 1i * w) + I * 1i * gamma / w;
end

% Estrazione dei componenti secondo la convenzione fisica: sigma = s1 - i*s2
s1 = real(sigma);
s2 = imag(sigma); 
end
