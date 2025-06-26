% Parámetros de la planta
kp = 1.304;
L = 0.12678;
tau = 0.9969;
t0 = L/tau; % Tiempo muerto normalizado

kc = (tau/(kp*L)) * ( (4/3) + (L/4*tau))
Ti = L*( (32 + (6*L)/tau) / (13 + (8*L)/tau) )

kcmurril = 0.984/kp * t0 ^ (-0.986) 
Timurril = tau/0.608 * t0 ^ (0.707) 