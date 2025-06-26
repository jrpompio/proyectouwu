% Parámetros de la planta
kp = 1.304;
L = 0.12678;
tau = 0.9969;
t0 = L/tau; % Tiempo muerto normalizado

% Constantes de Kaya y Sheib IAE (Controlador clásico - Regulador)
a = 0.98089;
b = -0.76167;
c = 0.91032;
d = -1.03211;

% Cálculo de Kc
Kc_kp = a * t0^b;
Kc = Kc_kp / kp;

% Cálculo de Ti
Ti_tau = (1/c) * t0^(-d);
Ti = Ti_tau * tau;

% Mostrar resultados
fprintf('Kc = %.3f\n', Kc);
fprintf('Ti = %.3f s\n', Ti);

% Modelo de la planta (para simulación)
s = tf('s');
P = kp * exp(-L*s) / (tau*s + 1);

% Controlador PI
C = Kc * (1 + 1/(Ti*s));

% Sistema en lazo cerrado
T = feedback(C*P, 1); % Lazo cerrado con retroalimentación negativa

% Simulación de respuesta a una perturbación
t = 0:0.01:10;
u = ones(size(t)); % Perturbación escalón
[y, t] = lsim(T, u, t);
plot(t, y);
xlabel('Tiempo (s)');
ylabel('Salida');
title('Respuesta del sistema a una perturbación');
grid on;

% Cálculo del IAE
error = abs(1 - y); % Error respecto al valor deseado (1 para regulador)
IAE = trapz(t, error);
fprintf('IAE = %.3f\n', IAE);
