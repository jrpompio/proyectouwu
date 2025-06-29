clc; close all;

% Define el modelo continuo
num = 1.304;
den = [0.9969 1];
delay = 0.1267;
Psitb = tf(num, den, 'InputDelay', delay);

% Define el tiempo de muestreo
Ts = 0.0265; 

% Discretiza el modelo
discret = c2d(Psitb, Ts, 'zoh');

% Abre sisotool para el modelo discreto
sisotool(discret);
