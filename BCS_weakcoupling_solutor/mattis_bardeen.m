function [s1, s2] = mattis_bardeen(hw, Delta, kT)
% MATTIS_BARDEEN  sigma1/sigma_n and sigma2/sigma_n 
%   hw    : photon energy(ies) \hbar*\omega  (same units as Delta, kT; hw > 0)
%   Delta : gap at the chosen temperature
%   kT    : k_B*T (same units; kT = 0 allowed)
%   s1,s2 : sigma1/sigma_n, sigma2/sigma_n (same size as hw)

s1 = zeros(size(hw));
s2 = zeros(size(hw));

fermi = @(x) 0.5*(1 - tanh(x/(2*kT)));          % 1/(exp(x/kT)+1), overflow-safe
if kT > 0
    th = @(x) tanh(x/(2*kT));                   % = 1 - 2 f(x)
else
    th = @(x) ones(size(x));                    % T = 0 -> f(x) = 0 (arguments are always > 0 here)
end

opts = {'RelTol',1e-10,'AbsTol',1e-12};

for k = 1:numel(hw)
    w = hw(k);
    N = @(E) E.^2 + Delta^2 + w*E;

    % ---------- sigma1, term 1: thermally excited quasi-particles ----------
    I1 = 0;
    if kT > 0
        umax = acosh(1 + 40*kT/Delta);          % E = Delta*cosh(u), up to Delta+40kT
        g1 = @(u) (fermi(Delta*cosh(u)) - fermi(Delta*cosh(u) + w)) ...
            .* N(Delta*cosh(u)) ./ sqrt((Delta*cosh(u) + w).^2 - Delta^2);
        I1 = integral(g1, 0, umax, opts{:});
    end
    s1(k) = 2/w * I1;

    % ---------- sigma1, term 2: photon-induced (only for hw > 2*Delta) ----------
    if w > 2*Delta
        a = Delta - w;  b = -Delta;
        E  = @(t) (a+b)/2 + (b-a)/2*sin(t);
        g2 = @(t) th(E(t)+w) .* N(E(t)) ./ ( sqrt(Delta - E(t)) .* sqrt(E(t) + w + Delta) );
        s1(k) = s1(k) - 1/w * integral(g2, -pi/2, pi/2, opts{:});
    end

    % ---------- sigma2 ----------
    b = Delta;
    if w > 2*Delta
        a = -Delta;
        E  = @(t) (a+b)/2 + (b-a)/2*sin(t);
        g3 = @(t) th(E(t)+w) .* N(E(t)) ./ sqrt((E(t)+w).^2 - Delta^2);
    else
        a = Delta - w;
        E  = @(t) (a+b)/2 + (b-a)/2*sin(t);
        g3 = @(t) th(E(t)+w) .* N(E(t)) ./ ( sqrt(E(t) + Delta) .* sqrt(E(t) + w + Delta) );
    end
    s2(k) = 1/w * integral(g3, -pi/2, pi/2, opts{:});
end
end