

function aprox_min_2_polinomio

    x = [ 1:10 ]';
    y = [ 22.1 ; 22.4 ; 21.9 ; 22.5 ; 22.2 ; 22.8 ; 22.7 ; 22.3 ; 22.9 ; 22.8 ];

    % queremos ajustar al polinomio y = a0 + a1 * x + a2 * x2

    phi = [ ones(1, 10) ; x' ; (x.^2)' ]';

    % phi' * phi * a = phi' * y

    A = phi' * phi;
    b = phi' * y;

    ec = linsolve(A,b);


    % Mostrar los coeficientes
    disp('Coeficientes del polinomio:');
    disp(ec);

    % Generar los valores ajustados
    x_fit = linspace(min(x), max(x), 100)';
    y_fit = ec(1) + ec(2) * x_fit + ec(3) * x_fit.^2;

    % Graficar los datos originales y la curva ajustada
    figure;
    plot(x, y, 'x', 'DisplayName', 'Datos originales');
    hold on;
    plot(x_fit, y_fit, '-', 'DisplayName', 'Ajuste polinomial');
    legend('show');
    title('Ajuste de datos a un polinomio de segundo grado');
    hold off;

endfunction
