function parcial2d

% Parámetros iniciales
Dt = 0.16755;
A1 = 1/4;

% Leer datos del archivo de texto
filename = 'registro-13-nov.txt';
data = load(filename);
tn = (0:length(data)-1) * Dt;
gn = data;

% Gráfica de la función discreta gn en función de tn
figure;
plot(tn, gn, 'o-');
xlabel('tn');
ylabel('g_n');
title('Gráfica de la función discreta g_n en función de t_n');

% Encontrar g_max y t_gmax en el rango t_n desde 0 hasta t(N/3)
N = length(tn)
tn_limit = tn(ceil(N/3));
indices = find(tn <= tn_limit);
[g_max, idx_max] = max(gn(indices));
tg_max = tn(indices(idx_max));
fprintf('g_max = %.2f, tg_max = %.2f\n', g_max, tg_max);

% Calcular la TDF de g_n y graficar módulo
G_tdf = fft(gn);
f = (0:N-1) / (N * Dt); % Frecuencias asociadas
mod_G_tdf = abs(G_tdf);

figure;
plot(f, mod_G_tdf, 'o-');
xlabel('Frecuencia (k)');
ylabel('|G{tdf}(k)|');
title('Módulo de la Transformada de Fourier de g_n');

% Encontrar G_tdf_MAX
G_tdf_MAX = max(mod_G_tdf);
fprintf('G_tdf_MAX = %.2f\n', G_tdf_MAX);

% Paso 4: Determinación automática de kc
threshold = 0.05 * G_tdf_MAX; % Umbral del 5% de G_tdf_MAX
valid_indices = find(mod_G_tdf > threshold); % ?�ndices con módulos mayores al umbral
kc = max(valid_indices(valid_indices <= 10)); % kc debe estar en los primeros 10 valores
G_tdf_kc = mod_G_tdf(kc); % Modulo asociado a kc
fprintf('kc = %d, G_tdf_kc = %.2f\n', kc, G_tdf_kc);

% Paso 5: Aproximación de Mínimos Cuadrados de g_n usando kc frecuencias
Dw = 2 * pi / (N * Dt); % Frecuencia angular
Pa = zeros(size(tn)); % Inicializar la aproximación

% Construcción de la matriz de base trigonométrica
Phi = ones(N, 2*kc + 1); % Matriz Phi con la base de funciones
for k = 1:kc
    Phi(:, 2*k) = cos(k * Dw * tn);
    Phi(:, 2*k + 1) = sin(k * Dw * tn);
end

% Resolución por mínimos cuadrados para obtener los coeficientes
alpha = Phi \ gn; % Resuelve Phi * alpha = gn

% Construir la función aproximada Pa(tn)
Pa = Phi * alpha;

% Graficar la aproximación Pa
figure;
plot(tn, Pa, 'x-', 'DisplayName', 'P_a');
xlabel('t_n');
ylabel('P_a');
title('Aproximación de Mínimos Cuadrados Pa(t_n)');

% Encontrar Pa_max y t_Pa_max en el rango t_n desde 0 hasta t(N/3)
[Pa_max, idx_Pa_max] = max(Pa(indices));
tPa_max = tn(indices(idx_Pa_max));
fprintf('Pa_max = %.2f, tPa_max = %.2f\n', Pa_max, tPa_max);

% Paso 6: Generar la función discreta de respuesta a impulso h_n
p = kc * Dw;
hn = A1 * exp(-p * tn);

% Paso 7: Convolución entre h_n y g_n para obtener y_c(t_n)
yc_full = conv(hn, gn) * Dt; % Escala la convolución por Dt para la correcta amplitud
yc = yc_full(1:N); % Recorta para mantener el tamaño original de N puntos

% Encontrar yc_max y t_yc_max en el rango t_n desde 0 hasta t(N/3)
[yc_max, idx_yc_max] = max(yc(indices));
tyc_max = tn(indices(idx_yc_max));
fprintf('yc_max = %.2f, tyc_max = %.2f\n', yc_max, tyc_max);

% Paso 8: Gráfica simultánea de g_n, Pa y y_c
figure;
hold on;
grid on;
plot(tn, gn, '-r', 'DisplayName', 'g_n');
plot(tn, Pa, '-b', 'DisplayName', 'P_a');
plot(tn, yc, '-g', 'DisplayName', 'y_c');
xlabel('t_n');
ylabel('Amplitud');
title('Comparación de g_n, P_a y y_c');
legend;
hold off;
disp('')

endfunction
