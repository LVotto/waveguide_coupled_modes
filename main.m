% Reproducing Griffel's differential equation (12) for
% coupled local amplitudes.
clear;
close all;
sympref('HeavisideAtOrigin', 0);

%   Defining parameters
lambda0 = 1.5E-6; k0 = 2 * pi / lambda0;    % Wave number in free space
omega = k0 * physconst('LightSpeed');       % Harmonic dependency
permt0 = 8.8541878128E-12;                  % Free space permittivity
permb0 = 1.25663706212E-6;                  % Free space permeability
d2 = 1E-6; d4 = .3E-6;                      % Center layer widths
ns = [1, 3.3, 3.2, 3.5, 3];                 % Reffractive index
s3 = 6E-7;                                  % Waveguide separation
ws = [0, d2, 2 * s3, d4, 0];                % Layer widths
dgr = 1E-7;                                 % Grating depth
pert_x = ws(2);                             % Grating center
contrast = abs(ns(2) ^ 2 - ns(3) ^ 2); % Grating amplitude (index contrast)
%   Separating parameters for guides a and b
nsa = ns(1:3);              nsb = ns(3:5); 
% nsa = [0, ns(2), 0]; nsb = [0, ns(4), 0]; % maybe this will work??
wsa = cat(2, ws(1:2), 0);   wsb = cat(2, 0, ws(4:5));

betas_a = tmt_betas(k0, nsa, wsa, "te");
betas_b = tmt_betas(k0, nsb, wsb, "te");
ba = betas_a(end); bb = betas_b(end);       % TE_0

phi = sum(wsa) + 2 * s3;        % Where guide b's second layer starts

%   TMT coefficients
coeffs_a = double(tmt_coeffs(ba, k0, wsa, nsa, "te", true));
coeffs_b = double(tmt_coeffs(bb, k0, wsb, nsb, "te", true));

%   Determining the coupling coefficients analitically
t0a = 0; t0b = phi;
lo = -Inf; hi = Inf; delta = true;
Ktaa = omega / 4 * integrate_product(lo, hi, k0, ba, ba, ws, wsa, wsa, ...
                                     ns, nsa, nsa, coeffs_a, coeffs_a, ...
                                     t0a, t0a, delta);
Ktab = omega / 4 * integrate_product(lo, hi, k0, ba, bb, ws, wsa, wsb, ...
                                     ns, nsa, nsb, coeffs_a, coeffs_b, ...
                                     t0a, t0b, delta);
Ktba = omega / 4 * integrate_product(lo, hi, k0, bb, ba, ws, wsb, wsa, ...
                                     ns, nsb, nsa, coeffs_b, coeffs_a, ...
                                     t0b, t0a, delta);
Ktbb = omega / 4 * integrate_product(lo, hi, k0, bb, bb, ws, wsb, wsb, ...
                                     ns, nsb, nsb, coeffs_b, coeffs_b, ...
                                     t0b, t0b, delta);
                                 
lo = pert_x - dgr / 2; hi = pert_x + dgr / 2; delta = false;
Khaa = omega / 4 * integrate_product(lo, hi, k0, ba, ba, ws, wsa, wsa, ...
                                     ns, nsa, nsa, coeffs_a, coeffs_a, ...
                                     t0a, t0a, delta);
Khab = omega / 4 * integrate_product(lo, hi, k0, ba, bb, ws, wsa, wsb, ...
                                     ns, nsa, nsb, coeffs_a, coeffs_b, ...
                                     t0a, t0b, delta);
Khba = omega / 4 * integrate_product(lo, hi, k0, bb, ba, ws, wsb, wsa, ...
                                     ns, nsb, nsa, coeffs_b, coeffs_a, ...
                                     t0b, t0a, delta);
Khbb = omega / 4 * integrate_product(lo, hi, k0, bb, bb, ws, wsb, wsb, ...
                                     ns, nsb, nsb, coeffs_b, coeffs_b, ...
                                     t0b, t0b, delta);
Khaa = permt0 * Khaa; Khab = permt0 * Khab;
Khba = permt0 * Khba; Khbb = permt0 * Khbb;
                                 
c = integrate_product(-Inf, Inf, k0, ba, bb, ws, wsa, wsb,  ...
                      ns, nsa, nsb, coeffs_a, coeffs_b,     ...
                      t0a, t0b, false);
c = c * (ba + bb) / 4 / omega / permb0;

Kaa = (Ktaa - c * Ktba) / (1 - c ^ 2);
Kbb = (Ktbb - c * Ktab) / (1 - c ^ 2);
Kab = (Ktab - c * Ktbb) / (1 - c ^ 2);
Kba = (Ktba - c * Ktaa) / (1 - c ^ 2);
ga = ba + Kaa; gb = bb + Kbb;
M = [[ga, Kab]; [Kba, gb]];

kaa = (Khaa - c * Khab) / (1 - c ^ 2);
kbb = (Khbb - c * Khba) / (1 - c ^ 2);
kab = (Khba - c * Khbb) / (1 - c ^ 2);
kba = (Khab - c * Khaa) / (1 - c ^ 2);
N = [[kaa, kab]; [kba, kbb]];

L = 2 * pi / (ga - gb);
dn2 = @(z) 2 / pi * contrast * cos(2 * pi / L * z);

%   Defining the differential equation itself
% Matrix A of ODE u' = A * u
% A = 1j * (M + dn2 N)
A = @(z) 1j * (M + dn2(z) * N);

a0 = 0;            % Initial condition A(0)
zend = 10E-3;       % Greatest z value
zspan = [0 zend];  % z values to evaluate
u0 = [a0; sqrt(1 - a0 ^ 2)];  % Initial condition.
%   The ODE
[z, u] = ode45(@(z, u) (A(z) * u), zspan, u0, ...
                odeset('RelTol', 1E-9, 'AbsTol', 1E-11));

%   Results:
zscale = 1E3;
plot(z * zscale, abs(u) .^ 2, 'LineWidth', 2);
ylim([0, 1]);
grid();
xlabel('$z[mm]$', 'Interpreter', 'latex');
legend({'$|a(z)|^2$', '$|b(z)|^2$'}, 'Interpreter', 'latex');
pbaspect([zend zend / 3 zend / 3]);
hold on
x = transp(z);
u1 = abs(transp(u(:, 1)) .^ 2);
u2 = abs(transp(u(:, 2)) .^ 2);
if a0 ^ 2 > 1/2
    pl = 2;
else
    pl = 1;
end
[m, i] = max(abs(u(:, pl)) .^ 2);
line(z(i) * [1, 1] * zscale, [min(ylim), max(ylim)],...
    'color', 'red', 'HandleVisibility', 'off', 'LineWidth', 1, ...
    'LineStyle', '--');
hold off