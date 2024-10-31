

function aprox_vector_r3

    y = [ 2.3 , 4.5 , 7.3 ];
    phi = [ 1 , 2 , 3];

    % r = y - a * phi
    % <r,phi> = 0

    % a = <r,phi> / <phi,phi>
    a = dot(y,phi) / dot(phi,phi);

    % p = a * phi
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
