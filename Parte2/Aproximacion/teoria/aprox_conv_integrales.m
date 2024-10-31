

function aprox_conv_integrales

    format long
    I_exacto = 0.3183;

    % Im
    % x = [ 5 ; 10 ; 20 ; 40 ; 80 ];
    % y = [ 0.2656876 ; 0.2926551 ; 0.3056462 ; 0.312019 ; 0.31517 ];

    % It
    % x = [ 5 ; 10 ; 20 ; 40 ; 80 ];
    % y = [ 0.31568758 ; 0.31765512 ; 0.31814624 ; 0.31826898 ; 0.3183 ];


    n = length(x);

    phi = [ ones(1,n) ; (1 ./ (x.^2))' ]';

    A = phi'*phi;
    b = phi'*y;

    sol = linsolve(A,b);

    I_conv = sol(1);

    proporcion = (I_conv / I_exacto)

    # disp(['It = ', num2str(proporcion),'% Iex']);





endfunction


0.31568758 ; 0.31765512 ; 0.31814624

