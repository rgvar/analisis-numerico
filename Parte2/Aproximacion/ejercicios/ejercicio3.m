function ejercicio3
    clc,clear
    graficar = 1;

    # g(t) = 10cos(0.5t) + 2cos(5t)
    format long
    # parametros
    N = 100;
    tf = 3*(2*pi/0.5); t0 = 0;
    Dt = (tf-t0)/N
    m = (N/2)-1

    w0 = 2*pi/N
    dw = w0/Dt

    ## discretizo la función
    f = @(t) 10*cos(0.5*t) + 2*cos(5*t);

    t = zeros(N,1);
    g = zeros(N,1);
    t(1,1) = 0;
    g(1,1) = f(0);
    for j=2:N
        t(j) = j * Dt;
        g(j) = f(t(j));
    endfor

    ## armo la base
    FI = zeros(N, m*2+1);
    for j=1:N
        for k=1:m
            FI(j,k*2) = cos(k*w0*(j-1));
            FI(j,k*2+1) = sin(k*w0*(j-1));
        endfor
    endfor

    # A alfa = b
    A = diag(diag(FI'*FI));
    b = FI'* g;
    alfa = linsolve(A,b);

    P = FI * alfa;
    r = g - P;

    if (graficar == 1)
        figure(1)
        plot(t,g,'r*-',t,P,'b*-',t,r,'r*-')
        grid on
        title('g(t),P(t),r(t)')
    endif


endfunction
