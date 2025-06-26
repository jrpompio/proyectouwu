clc; close all;

% 1. Leer los datos -------------------------------------------------------
archivo = 'Flujo_delta_60a80recortado.csv';

opts = detectImportOptions(archivo, ...
    'Delimiter',';', ...
    'DecimalSeparator',',');

raw = readtable(archivo, opts);

t = raw{:,1};        % tiempo
u = raw{:,2}-60;     % entrada (señal de control)
y = raw{:,3}-53.36;  % salida (respuesta del sistema)

P123c = tf(1.32, [0.910 1], 'InputDelay', 0.138);

num = 1.304;
den = [0.9969 1];
delay = 0.1267;

Psitb = tf(num, den, 'InputDelay', delay);

%2. Crear objeto iddata --------------------------------------------------
Ts   = median(diff(t));     % tiempo de muestreo
data = iddata(y, u, Ts, 'TimeUnit','s');

% 4. Simular la salida del modelo -----------------------------------------
yP123c = lsim(P123c, u, t);
yPsitb = lsim(Psitb, u, t);

% % % 5. Graficar -------------------------------------------------------------
% figure;
% plot(t, u, 'k:', 'LineWidth', 1.5); hold on;
% plot(t, y, 'b', 'LineWidth', 1);
% plot(t, yP123c, 'r--', 'LineWidth', 1.5);
% plot(t, yPsitb, 'Color', [0.4 0.7 0.3], 'LineWidth', 1.8);
% xlabel('Tiempo [s]');
% ylabel('Amplitud normalizada');
% legend('u(t)', 'y(t)', 'y_{123c} (t)', 'y_{SITb} (t)');
% xlim([209.5 230.6]);
% % xlim([206. 230.6]);
% grid on

x 
n_activos = sum((u == 20));
yActivo = y(end - n_activos + 1 : end);
tActivo = t(end - n_activos + 1 : end);
yP123cActivo = yP123c(end - n_activos + 1 : end);
yPsitbcActivo = yPsitb(end - n_activos + 1 : end);

IAEP123c = abs(yActivo - yP123cActivo);
IAEPsitbc = abs(yActivo - yPsitbcActivo);

% IAE acumulado para cada modelo
IAE123c_acum = cumsum(IAEP123c) * Ts;
IAEsitbc_acum = cumsum(IAEPsitbc) * Ts;

% % Graficar la evolución del IAE
% figure;
% plot(tActivo, IAE123c_acum, 'r-', 'LineWidth', 1.5); hold on;
% plot(tActivo, IAEsitbc_acum, 'Color', [0.4 0.7 0.3], 'LineWidth', 1.8);
% xlabel('Tiempo [s]');
% ylabel('IAE acumulado');
% legend('IAE 123c', 'IAE SITb');
% title('Evolución del IAE en el tramo activo');
% grid on;

discret = c2d(Psitb, Ts, 'zoh');
sisotool(discret);
fprintf('IAE total para P123c:   %.4f\n', IAE123c_acum(end))
fprintf('IAE total para Psitbc:  %.4f\n', IAEsitbc_acum(end))
