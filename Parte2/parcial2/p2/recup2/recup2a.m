function recup2

    clc, clear

    #CARGADO DE ARCHIVO
    g=load("datos.txt", "-ascii");
    N=length(g)

    ##DECLARACION DE VARIABLES Y VECTORES
    Dim=N+1;
    Dt=0.083776;
    t0=0;
    tf=Dt*N;
    t=zeros(1,Dim);
    gNorma=0;

    ##LLENADO VECTOR t
    for i=1:Dim-1
      t(i) = (i-1) * Dt;
      gNorma=gNorma+(g(i)*g(i));
    endfor

    ##IMPRESION DE gNorma
    gNorma2=sqrt(gNorma);
    disp(['gNorma2 = ',num2str(gNorma2,'%.2f')]);

    ##IMPRESION FUNCION DISCRETA
    figure(1)
    title("Figura dato")
    plot(t(1:N),g,'-r')
    grid on

  ##TRANSFORMADA DISCRETA DE FOURIER(TDF)
    gTDF=fft(g,N);

    ##CALCULA E IMPRIME dw
    dw=2*pi/(Dt*N);
    disp(['dw = ',num2str(dw,'%.4f')]);

    ##CALCULAMOS modG
    modG=0;
    for i=1:N
      modG=modG+(abs(gTDF(i))*abs(gTDF(i)));
    endfor

    ##IMPRIME modG2
    modG2=sqrt(modG);
    disp(['modG2 = ',num2str(modG2,'%.4f')]);

    figure(2)
    title("modulo G")
    stem(abs(gTDF(1:110)))
    grid on

  ##CONVOLUCION
    kc=8;
    aUno=1/3;
    w=kc*dw;
    p=0.3*w;
    h=zeros(Dim,1);

     for i=1:N
        h(i,1)=aUno*((e^(-p*t(i)))/w)*sin(w*t(i));
     endfor

    gConv=Dt*conv(h,g);

    ##IMPRIMIMOS LO PEDIDO
    figure(3)
    title( "h")
    plot(t(1:N),h(1:N,1),'-b')
    grid on

    figure(4)
    title( "Convolucion")
    plot(t(1:N),gConv(1:N,1),'-b')
    grid on

    figure(5)
    title( "Convolucion y h")
    plot(t(1:N),g(1:N),'-r',t(1:N),gConv(1:N),'-b')
    grid on


endfunction
