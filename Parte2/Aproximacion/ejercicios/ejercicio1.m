
function ejercicio1

    clc,clear

    graficos = 1;

    # función
    N=128;
    tf=10;
    t0=0;
    Dt=(tf-t0)/N;
    m=(N/2)-1;
    ncol=m*2+1;

    # períodos
    Tg=10;
    Tm=5;
    Tp=10;

    # frecuencias
    w0=2*pi/N;
    dw=w0/Dt;

    t=zeros(N,1);
    g=zeros(N,1);

    # construcción función discreta
    for j=1:N
        t(j) = (j-1)*Dt;
        g(j) = sierra(t(j),Tm);
    endfor

    # construcción base
    FI = zeros(N, ncol);
    FI(:,1) = 1;
    for j=1:N
        for k=1:m
            FI(j,k*2) = cos(k * w0 * (j-1));
            FI(j,k*2+1) = sin(k * w0 * (j-1));
        endfor
    endfor

    A = diag(diag(FI'*FI));
    b = FI' * g;
    alfa = linsolve(A,b);

    P = FI * alfa;
    r = g - P;

    # gráficos
    if (graficos == 1)

        figure(1)
        plot(t,g,'ro-')
        grid on
        title('g(t)')

        figure(2)
        plot(t, P, 'bo-')
        grid on
        title('P(t)')

        figure(3)
        plot(t,g, 'r-',t,P,'b-')
        grid on
        title('comparación')

    endif



endfunction

function y = pulso(t,m)
    if (t < m)
        y = 1;
    else
        y = 0;
    endif
endfunction

function y = sierra(t,m)
    if (t < m)
        y = t;
    else
        y = 0;
    endif
endfunction
