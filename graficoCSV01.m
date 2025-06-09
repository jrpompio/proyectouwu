clc; close all;

%% 1. Leer los datos -------------------------------------------------------
archivo = 'Flujo_delta_60a80.csv';      % <- tu archivo CSV

opts = detectImportOptions(archivo, ...
    'Delimiter',';', ...
    'DecimalSeparator',',');

raw = readtable(archivo, opts);

t = raw{:,1};                % tiempo
u = raw{:,2};                % entrada (señal de control)
y = raw{:,3};                % salida (respuesta del sistema)


%% 2. Crear objeto iddata --------------------------------------------------
Ts   = median(diff(t));      % tiempo de muestreo
data = iddata(y, u, Ts, 'TimeUnit','s');

%% 3. Dividir en estimación y validación ----------------------------------
porcEst = 0.70;                      % 70 %
N       = length(data.y);            % total de muestras
Nest    = floor(porcEst*N);          % número para estimar

ze = data(1:Nest);                   % *estimación*   (primer 70 %)
zv = data(Nest+1:end);               % *validación*   (último 30 %)

% --- Si prefieres mezclar las muestras al azar, usa: ---------------------
% rng default                         % para reproducibilidad
% idx = randperm(N);
% ze  = data(idx(1:Nest));
% zv  = data(idx(Nest+1:end));
% zv  = sort(zv);                     % vuelve a ordenar por tiempo
% -------------------------------------------------------------------------

%% 4. Calibrar modelos -----------------------------------------------------
% 4.1 ARX discreto  (na, nb, nk) = (2, 2, 1)
na = 2; nb = 2; nk = 1;
sysARX = arx(ze, [na nb nk]);

% 4.2 Función de transferencia continua (2 polos, 1 cero, sin retardo)
np = 2; nz = 1; ioDelay = 0;
sysTF = tfest(ze, np, nz, ioDelay, 'Ts', 0);  % 'Ts',0 = modelo continuo

%% 5. Validación cruzada ---------------------------------------------------
figure
compare(zv, sysARX, sysTF)
legend('Datos de validación','ARX','TF','Location','best')

% Ajuste (porcentaje que entrega compare)
[~, fitARX] = compare(zv, sysARX);
[~, fitTF ] = compare(zv, sysTF );

fprintf('\n%% de ajuste (NRMSE inverso)\n');
fprintf('  ARX : %.2f %%\n', fitARX);
fprintf('  TF  : %.2f %%\n', fitTF);

%% 6. Visualizar la dinámica identificada ----------------------------------
figure, bode(sysTF),  grid on, title('Bode – Modelo TF continuo')
figure, step(sysTF, max(t)-min(t)), grid on, title('Respuesta al escalón – Modelo TF')
