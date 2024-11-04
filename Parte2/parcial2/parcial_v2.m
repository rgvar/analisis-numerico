function procesar_senal
    % Limpiar entorno
    clear; clc;

    % Parámetros iniciales
    Dt = 0.16755;  % Intervalo de muestreo en segundos
    A1 = 0.25;     % Escala para la respuesta al impulso

    % Cargar datos y definir tiempo discreto
    data = load('registro_4_nov.txt');
    N = length(data);
    t = (0:N-1) * Dt;  % Vector de tiempo
    g = data;          % Señal discreta

    % Graficar señal discreta
    figure;
    plot(t, g, 'o-');
    xlabel('Tiempo (s)');
    ylabel('Amplitud g[n]');
    title('Señal Discreta g[n] vs Tiempo');

    % Calcular máximo de la señal en el primer tercio del tiempo
    limite_tiempo = t(ceil(N / 3));
    indices_tercio = find(t <= limite_tiempo);
    [g_max, idx_g_max] = max(g(indices_tercio));
    t_g_max = t(indices_tercio(idx_g_max));
    fprintf('g_max = %.2f, t_g_max = %.2f\n', g_max, t_g_max);

    % Calcular la Transformada Discreta de Fourier (TDF) y su magnitud
    G_fft = fft(g);
    frec = (0:N-1) / (N * Dt);   % Vector de frecuencias
    G_modulo = abs(G_fft);

    % Graficar módulo de la TDF
    figure;
    plot(frec, G_modulo, 'o-');
    xlabel('Frecuencia (Hz)');
    ylabel('|G(f)|');
    title('Módulo de la Transformada de Fourier de g[n]');

    % Determinar frecuencia de corte kc (5% del máximo)
    G_max = max(G_modulo);
    umbral = 0.05 * G_max;
    frec_significativas = find(G_modulo > umbral);
    kc = max(frec_significativas(frec_significativas <= 10));
    G_kc = G_modulo(kc);
    fprintf('kc = %d, |G(kc)| = %.2f\n', kc, G_kc);

    % Aproximación de g usando kc frecuencias (método de mínimos cuadrados)
    Dw = 2 * pi / (N * Dt);    % Frecuencia angular fundamental
    Phi = ones(N, 2*kc + 1);   % Base trigonométrica

    % Construcción de la matriz base
    for k = 1:kc
        Phi(:, 2*k) = cos(k * Dw * t);
        Phi(:, 2*k + 1) = sin(k * Dw * t);
    end

    % Calcular coeficientes de mínimos cuadrados y aproximación
    coeficientes = Phi \ g;
    g_aprox = Phi * coeficientes;

    % Graficar la señal aproximada
    figure;
    plot(t, g_aprox, 'x-', 'DisplayName', 'g_{aprox}');
    xlabel('Tiempo (s)');
    ylabel('Aproximación de g');
    title('Aproximación de g[n] usando Mínimos Cuadrados');

    % Máximo de la aproximación en el primer tercio del tiempo
    [g_aprox_max, idx_g_aprox_max] = max(g_aprox(indices_tercio));
    t_g_aprox_max = t(indices_tercio(idx_g_aprox_max));
    fprintf('g_aprox_max = %.2f, t_g_aprox_max = %.2f\n', g_aprox_max, t_g_aprox_max);

    % Generar función de respuesta al impulso h[n]
    p = kc * Dw;
    h = A1 * exp(-p * t);

    % Convolución entre h y g para obtener y_c[n]
    y_c_full = conv(h, g) * Dt;
    y_c = y_c_full(1:N);  % Recortar para mantener tamaño original

    % Máximo de y_c en el primer tercio del tiempo
    [y_c_max, idx_y_c_max] = max(y_c(indices_tercio));
    t_y_c_max = t(indices_tercio(idx_y_c_max));
    fprintf('y_c_max = %.2f, t_y_c_max = %.2f\n', y_c_max, t_y_c_max);

    % Graficar comparativa de g, g_aprox y y_c
    figure;
    hold on;
    plot(t, g, 'r-', 'DisplayName', 'g[n]');
    plot(t, g_aprox, 'b-', 'DisplayName', 'g_{aprox}');
    plot(t, y_c, 'g-', 'DisplayName', 'y_c');
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    title('Comparación de g[n], g_{aprox} y y_c');
    legend;
    hold off;

end

