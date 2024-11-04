function parcial2Burky_corregido
  # función g con datos brindados
  gn = load("registro-13-nov.txt", "-ascii");

  # datos iniciales
  Dt=0.16755;
  t0=0;
  A1=1/4;

  N=length(gn); # para saber cantidad de datos, son 1024

  # vector tn (abscisas)
  for k=1:N
    tn(k)=t0+(k-1)*Dt;
  endfor

  ### GENERACIÓN DE LA FUNCIÓN DISCRETA DATO (ejercicio 1) ###
  figure(1)
  plot(tn,gn,'r')
  title('Generación de la función discreta dato')
  grid on
  title('figura 1')

  # valores máximos de tn y gn en el rango de tn

  tg_max=0;
  g_max=0;

  for i=1:(N/3)
    if(gn(i)>g_max)
      g_max = gn(i);
      tg_max = tn(i);
    endif
  endfor

  tg_max
  g_max

  ### TRANSFORMADA DISCRETA DE FOURIER (ejercicio 2) ###
  G_tdf=fft(gn,N);#Sacamos la transformada discreta de furier

  # para ver el módulo de cada una
  for k=1:N
    G_tdf_mod(k)=abs(G_tdf(k)); #para todos los valores obetenmos el valor absoluto
  endfor

  G_tdf_MAX = max(G_tdf_mod) # módulo de G_tdf punto maximo

  G_tdf_kc=0;
for k = 1:10
    if G_tdf_mod(k) > (0.05 * G_tdf_MAX)
        G_tdf_kc = G_tdf_mod(k);
        kc = k;
    endif
endfor

  G_tdf_kc
  kc

 figure(2)
  stem(abs(G_tdf))
  title('Transformada discreta de Fourier (módulo)')
  grid on
   title('figura 2')

### APROXIMACIÓN (ejercicio 3) ###
  m=kc;
  dim = 2*m+1;
  Tp = Dt*N; % no lo usé
  Dw = (2*pi)/Tp; % no lo usé
  w0 = (2*pi)/N;
  O = zeros(N,dim); # es la matriz phi
  for i=1:N
    O(i,1)=1;
    for j=1:m
      O(i,2*j)=cos(j*w0*(i-1));
      O(i,2*j+1)=sin(j*w0*(i-1));
    endfor
  endfor

  C = O'*gn; # es la matriz de coeficientes

  # vector a
  a=zeros(dim,1);
  for i=1:dim
    if (i==1)
      a(i,1) = C(1,1)/N;
   else
      a(i,1) = C(i,1)/(N/2);
   endif
  endfor

  Pa=O*a; # Pa es la matriz fi (acá la llamo O), multiplicada por el vector alfa

  tPa_max=0;
  Pa_max=0;

  for i=1:(N/3)
    if(Pa(i)>Pa_max)
      Pa_max = Pa(i);
      tPa_max = tn(i);
    endif
  endfor

  tPa_max
  Pa_max
  figure(3)
  plot(tn,Pa,'b')
  title('Aproximación con Pa')
  grid on
   title('figura 3')


  ### CONVOLUCIÓN (ejercicio 4) ###
  Tp = Dt*N;
  Dw = (2*pi)/Tp;
  p=kc*Dw;

  for i=1:N
    hn(i)=A1*exp((-p)*tn(i)); # función h para convolución
  endfor
  yc=Dt*conv(hn,gn); # queda en yc la convolución entre hn y gn

  figure(4)
  plot(tn,yc(1:N),'g')
  title('Convolución')
  grid on
   title('figura 4')

  tyc_max=0;
  yc_max=0;

  for i=1:(N/3)
    if(yc(i)>yc_max)
      yc_max = yc(i);
      tyc_max = tn(i);
    endif
  endfor

  tyc_max
  yc_max

  ### COMPARAR (ejercicio 5) ###
  figure(5)
  plot(tn,gn,'r',tn,Pa,'b',tn,yc(1:N),'g')
  title('Comparación final')
  grid on
   title('figura 5')

endfunction
