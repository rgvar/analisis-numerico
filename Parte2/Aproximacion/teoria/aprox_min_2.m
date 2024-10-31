

function aprox_min_2

    %{
        se quiere ajustar los datos (x,y) a un modelo lineal:
            y(i) = a0 + a1 * x(i)
        o   y(i) = phi * a
    %}
    # datos

    x = [ 0.5 ; 1.5 ; 2.5 ; 3.5 ];
    y = [ 1.3 ; 4.2 ; 5.7 ; 8.2 ];

    phi = [ 1 , 1 , 1 , 1 ; x' ]';

    a = [ a0=0 ; a1=0 ];
    %{
        Queremos encontrar los valores a0 y a1 que minimicen la dif entre
        los valores observados y(i) y los valores del modelo a0 + a1 * x(i)
    %}
    %{
        phi'*phi * a = phi' * y
                ==
        A11 * a0 + A12 * a1 = phi'1 * y1 (b1)
        A21 * a0 + A22 * a1 = phi'2 * y2 (b2)
    %}

    A = phi'*phi;
    b = phi'*y;

    sol = linsolve(A,b)





endfunction
