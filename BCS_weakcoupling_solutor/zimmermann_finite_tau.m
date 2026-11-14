function [s1, s2] = zimmermann_finite_tau(hw, Delta, kT, gamma, Umax)
% sigma(w)in BCS, Zimmermann-Bimonte 2010 (Eqs. 8, 19, 21). Richiede Delta > 0.
if nargin < 5, Umax = 10; end
opts  = {'RelTol', 1e-8, 'AbsTol', 1e-10};
sigma = zeros(size(hw));

for k = 1:numel(hw)
    w = hw(k);

    % integrando in u: ep = Delta*sinh(u), E = Delta*cosh(u), d(eps) = E du
    F = @(u) Integrand(Delta*sinh(u), w, Delta, kT, gamma) .* (Delta*cosh(u));

    % punti di rottura in u (dove Q_- = 0)
    ub  = acosh((w + Delta)/Delta);                 % E = w + Delta
    brk = 0;
    if w - Delta > Delta, brk(end+1) = acosh((w - Delta)/Delta); end   % E = w - Delta
    brk(end+1) = ub;
    brk(end+1) = max(Umax, ub + 6);                 % coda ~ e^{-2u}

    I = 0;
    for j = 1:numel(brk)-1
        I = I + integral_sin(F, brk(j), brk(j+1), opts);
    end

    % Drude + BCS: con I = integrale su [0, infinito) da parita' in ep, il 2 si cancella.
    sigma(k) = gamma/(gamma - 1i*w) + (1i*gamma/w)*I;
end

if nargout < 2
    s1 = sigma;
else
    s1 = real(sigma);  s2 = imag(sigma);
end
end

function I = integral_sin(F, a, b, opts)
% integrale di F su [a,b] con u = (a+b)/2 + (b-a)/2*sin(t): cancella le 1/sqrt agli estremi
m = (a+b)/2;  h = (b-a)/2;
G = @(t) clean(F(m + h*sin(t)).*h.*cos(t));
I = integral(G, -pi/2, pi/2, opts{:});
end

function v = clean(v)
v(~isfinite(v)) = 0;          % protegge i nodi che cadono esattamente sull'estremo (peso ~ cos t ~ 0) (sugg. claude)
end

% ---- Integrand ----

function v = Integrand(ep, w, Delta, kT, gamma)
e2 = ep.^2;
E  = sqrt(e2 + Delta^2);
if kT > 0
    th = tanh(E ./ (2*kT));
else
    th = ones(size(E));
end
Qp = sqrt((E + w).^2 - Delta^2);
s  = E - w;
Qm = complex(zeros(size(E)));
out = abs(s) > Delta;
Qm(out)  = -sign(s(out)).*sqrt(s(out).^2 - Delta^2);   % ramo Im Q > 0 dal paper
Qm(~out) = 1i*sqrt(Delta^2 - s(~out).^2);
Ap = E.*(E + w) + Delta^2;
Am = E.*(E - w) + Delta^2;
Gp = (e2.*Qp + (Qp + 1i*gamma).*Ap) ./ (Qp .* (e2 - (Qp + 1i*gamma).^2));
Gm = (e2.*Qm + (Qm + 1i*gamma).*Am) ./ (Qm .* (e2 - (Qm + 1i*gamma).^2));
v = (th ./ E) .* (Gp + Gm);
end