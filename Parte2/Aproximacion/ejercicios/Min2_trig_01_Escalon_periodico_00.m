function Min2_trig_01_Escalon_periodico_00
    clc, clear
    % se usa como g(t) una escalon periodica
    disp('el periodo de funcion escalon es: '),Tg=2*5
    disp('la 1era duración  del escalon es: '),Td=5
    disp('las Amplitudes    de g(t) son:'),A1=1, A2=0
    %
    disp ('el rango de t para g(t) es desde t0 hasta tf')
    t0=0,      tf= Tg     % da perfecto 2*pi*10
    disp('se discretiza en un N° de intervalos N= '), N=4,  %4 8 32
    disp(' el incremento de t resulta Dt= '), Dt=(tf-t0)/N
    disp('el periodo  es Tp=(tf-t0)= '),Tp=N*Dt
    disp('la frecuencia  Fundamental para Min 2,  es w0=2*pi/N'),w0=2*pi/N
    disp('el incrmento de frecuencia para Min 2,  es Dw=w0/Dt='),dw=w0/Dt
    disp('el multiplo máximo para las frecunecias es m=N/2= '),m_max=round(N/2)
    disp('se eligen para senos y cosenos desde 1 a m= '), m=-1+N/2
    disp('en la Base para Min 2, se toman 2m+1 elementos Nb'),Nb=2*m+1
    %
    compara_g_P=1      % si es 1 Arma la version continua de P(t)
    graficar=1;
    %
    t=zeros(N,1);
    g=zeros(N,1);
    % Armado de función g(t) discreta
    for j=1:N
        t(j)=t0+(j-1)*Dt;
        if t(j)< Td
            g(j)=A1;
        else
            g(j)=A2;
        end
    endfor


    if graficar==1
        figure (1)
        plot (t,g,'-or')
        grid on
        title ('g(t)')
    endif

    %%%%%%%%% INIICO DE MIN 2
    % Inicialización para Min2
    P=zeros(N,1);
    r=zeros(N,1);
    FI=zeros(N,2*m+1);
    alfa=zeros(2*m+1,1);
    kw0=zeros(1,m+1);
    c  =zeros(1,m+1);
    c_fft=zeros(1,m+1);
    %Elementos de la BASE para Min2 trigonometricos
    FI(:,1)=1;
    kw0(1)=0;
    kw(1)=0;
    for k=1:m
        kw0(k+1)= k*w0;
        kw(k+1)=k*dw;
        for j=1:N
            FI(j,2*k)  =cos(k*w0*(j-1));
            FI(j,2*k+1)=sin(k*w0*(j-1));
        end
    end
    FI
    FI'*FI
    D=diag(diag(FI'*FI));
    b=FI'*g;

    alfa = D \ b
    disp(' ')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%:')
    disp('los coeficientes b de MIn 2 son:')
    b
    disp('los coeficientes alfa de MIn 2 son:')
    alfa
    disp(' ')
    % FIN DE MIN 2
    %
    % COMPARA P de Min 2 con g(t)
    disp('COMPARA P de Min 2 con g(t):')
    P=FI*alfa;
    r=g-P;
    disp('la norma del residuo es'),Norma_r=norm(r,2)
    disp('la norma de g(t) es'),Norma_g=norm(g,2)

    if graficar==1
        figure (2)
        plot(t, P,'-xb',t,g,'or')
        grid on
        title ('g(t); P(t) discretos')

        figure (3)
        plot(t, P,'-xb',t,g,'or',t,r,'*r')
        grid on
        title ('g(t); P(t); r(t) discretos')
    endif

    disp(' ')
    % ARMA Rta en Frecuencia
    disp('Rta en Frecuencias')
    c(1)=abs(alfa(1,1));
    fase(1)=0;
    for k=2:m+1
        j=(k-1)*2;
        c(k)=(alfa(j,1)^2+alfa(j+1,1)^2)^0.5;
        fase(k)=atan(alfa(j+1,1)/alfa(j,1));
    end
    disp('las frecuencias para cada kw0 son:')
    kw0
    disp('las frecuencias para cada kw son:')
    kw
    disp('la amplitud de c para cada kw0 son:')
    c
    disp('la fase para cada kw0 son:')
    fase
    disp('')
    if graficar==1
        figure (4)
        subplot(2,1,1)
        stem(kw,c,'b')
        grid on
        title ('Amplitud para cada Frecuencia')
        subplot(2,1,2)
        plot(kw,fase,'-ob')
        grid on
        title ('Fase para cada Frecuencia')
    endif

    % FIN de  Rta en Frecuencia

    % Comparacion entre version discreta de g(t) y continua de P(t)
    if (compara_g_P ==1 && graficar == 1)
        disp('Esta armando la version continua de P(t)')
        tcom=0:1/100:2*Tp;
        dim_com=length(tcom)
        for j=1:dim_com
            P_tc(j)=alfa(1,1)*1;
            for k=1:m
                P_tc(j)=P_tc(j)+alfa(2*k)*cos(k*dw*(tcom(j)-t0));
                P_tc(j)=P_tc(j)+alfa(2*k+1)*sin(k*dw*(tcom(j)-t0));
            end
        end
        figure (10)
        plot(tcom,P_tc, t,g,'or')
        grid on
        title ('g(t) discreta y Aproximación P(t) continua')
    end
    % FIn de COMPARACION entre P de Min 2 con g(t)


endfunction

