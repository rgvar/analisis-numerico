

function aprox_vector_r2

    y = [ 3.2 , 6.8 ]; % vector dado
    phi=[ 1 , 2 ]; % vector base
    %{
        queremos proyectar el vector dado y, perteneciente a un plano,
        sobre una línea L definida por el vector base phi (ϕ)
        para obtener un vector p (resultado de la proyección)
    %}

    % r = y - p = y - a * phi

    % a = <y,phi> / <phi,phi>
    a = dot(y,phi)/ dot(phi,phi);

    % p = a * phi;
    p = a * phi;

    % r = y - p
    r = y - p;

    disp(['a = ', num2str(a)]);
    disp(['p = [ ', num2str(p), ' ]']);
    disp(['r = [ ', num2str(r), ' ]']);

    % verificar ortogonalidad
    if ( abs(dot(r,phi)) < 10^-10 )
        disp("r es ortogonal a ϕ => es correcta la proyección");
    endif

endfunction
