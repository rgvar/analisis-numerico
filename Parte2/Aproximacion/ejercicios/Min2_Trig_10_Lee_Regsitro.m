function Min2_trig_10_Lee_Resgistro
  clc, clear
  %%%%% se lee una g(t) desde archivo y se guarad en g como un vector
  %%  save( "Registro_Ej_4_reducido.txt", "-ascii", "g")
  
  gread=load( "Registro_240901.txt", "-ascii");   %tiene  16001 regsitros
  %gread=load( "EDO_Ord1_entrada.txt", "-ascii");  %tiene  90001 regsitros
  %gread=load( "EDO_Ord1_salida.txt", "-ascii");    %tiene  90001 regsitros
  
  N=length(gread)      % conviene detener y comprobar N
  N=16000
  disp('se discretizo en un N° de intervalos N= '), N
  disp('el regsitro tiene frec de muestreo fs= '), fs=16000
  disp('con intervalos  de t iguales a Dt: '), Dt=1/fs
  disp('el rango de t para g(t) es desde t0: '), t0=0
  disp('hasta el tiempo final   tf=Tp: '), Tp=t0+N*Dt, tf=Tp
  disp('se considera un periodo  Tp=(tf-t0)= N*Dt')
  disp('la Amplitud    de g(t) es:'),A1=1
  
 %%%%%%% Armado de función g(t) discreta
  for j=1:N
    t(j)=t0+(j-1)*Dt;
    g(j,1)=gread(j);
  endfor
 
  figure (1)
  plot (t,g,'r')
  grid on
  title ('g(t)')
  %
   %%%%%%%%% INICIO DE MIN 2
  disp('la frecuencia  Fundamental para Min 2,  es w0=2*pi/N'),w0=2*pi/N
  disp('el incrmento de frecuencia para Min 2,  es Dw=w0/Dt='),dw=w0/Dt
  disp('el multiplo máximo para las frecunecias es m=N/2= '),m_max=round(N/2)
  disp('se eligen para senos y cosenos desde 1 a m= '), m=-1+N/2, m=1500
  disp('en la Base para Min 2, se toman 2m+1 elementos Nb'),Nb=2*m+1
  %
  %
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
     kw(k+1)=k;
     for j=1:N
      FI(j,2*k)  =cos(k*w0*(j-1));
      FI(j,2*k+1)=sin(k*w0*(j-1));
     end
  end
  %FI
  %FI'*FI
  D=diag(diag(FI'*FI));
  b=FI'*g;
  for j=1:2*m+1
    alfa(j,1)=b(j)/D(j,j);
  end
  disp('los coeficientes b de MIn 2 son:')
  b 
  disp('los coeficientes alfa de MIn 2 son:')
  alfa
  %disp('la norma de alfa: Coef de MIn 2 es:'), norm(alfa,2)
  % FIN DE MIN 2
%
% COMPARA P de Min 2 con g(t) en versiones discretas
% 
  P=FI*alfa;
  r=g-P;
  disp('la norma del residuo es'),Norma_r=norm(r,2)
  disp('la norma de g(t) es'),Norma_g=norm(g,2)
   
  figure (2)
  plot(t, P,'-xb',t,g,'r')
  grid on
  title ('g(t); P(t) discretos')
  
  figure (3)
  plot(t, P,'-xb',t,g,'or',t,r,'*r')
  grid on
  title ('g(t); P(t); r(t) discretos')
% FIn de COMPARACION entre P de Min 2 con g(t)
%  
%%%%%%%%% ARMA RESPUESTA EN Frecuencia
  %disp('las frecuencias kw0 son:')
  %kw0
  c(1)=alfa(1,1);
  fase(1)=0;
  for k=2:m+1
    j=(k-1)*2;
    c(k)=(alfa(j,1)^2+alfa(j+1,1)^2)^0.5;
    fase(k)=atan(alfa(j+1,1)/alfa(j,1));
  end
  disp('la amplitud de c para cada kw0 son:')
  c
  disp('la Norma de la amplitud c(kw0) es:'),norm(c,2)
  disp('la fase para cada kw0 son:')
  fase
  disp('')
  
  figure (4)
  stem(kw, abs(c),'b') %bar(kw, abs(c),'-xb')
  grid on
  title ('Amplitud para cada k-ésima frecuencia')
  
  figure (5)
  plot(kw, fase,'b') %bar(kw, abs(c),'-xb')
  grid on
  title ('Fase para cada k-ésima frecuencia')%
  
  %%%%%%%%%%%%%%%%%%%%%% AMPLITUDES tipo FFT: FI'*g
  c_fft(1)=b(1,1);
  fase(1)=0;
  for k=2:m+1
    j=(k-1)*2;
    c_fft(k)=(b(j,1)^2+b(j+1,1)^2)^0.5;
    fase(k)=atan(b(j+1,1)/b(j,1));
  end
##  disp('la amplitud de c_fft para cada kw0 son:')
##  c_fft
 
  figure (6)
    stem(kw,abs(c_fft),'b')  %bar(kw,abs(c_fft),'xb')
    grid on
    title ('Amplitudes desde FI T*g')
 
endfunction
