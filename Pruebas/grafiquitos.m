% grafiquitos.m  –  versión con readtable
archivos = { '2-8 murril.lvm', '2-8 mendez.lvm' };

for k = 1:numel(archivos)
    T = readtable(archivos{k}, ...
                  'FileType', 'text', ...      % << clave
                  'Delimiter', '\t', ...
                  'ReadVariableNames', false, ...
                  'DecimalSeparator', ',');

    tiempo = T{:,1};
    datos  = T{:,2:end};

    figure;
    plot(tiempo, datos, 'LineWidth', 2);
    legend({'Referencia','Flujo','Señal de control'});
    xlabel('tiempo');
    ylabel('magnitud');
    grid on;
end
