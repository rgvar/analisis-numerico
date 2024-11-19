function globalLucas

  clc,clear

  ##CARGAMOS FUNCION
  gr=load("global.txt","-ascii");

  ##VARIABLES
  N=2688;
  rho=1.17;
  dt=0.02;
  t0=0;
  tf=(N-1)*dt;

  ##EJERCICIO 1
  ##EDO
  t2(1)=0;
  dx1(1)=0;
  dx2(1)=0;
  dx3(1)=0;
  h3(1)=0;

  ##SISTEMA EDO
  for k=1:N
    t2(k+1)=t2(k)+dt;
    g2(k)=1;
    dx1(k+1)=dx1(k)+dt*((-rho)*dx1(k)+0*dx2(k)+0*dx3(k)+1*g2(k));
    dx2(k+1)=dx2(k)+dt*(1*dx1(k)+(-rho)*dx2(k)+0*dx3(k)+0*g2(k));
    dx3(k+1)=dx3(k)+dt*(0*dx1(k)+1*dx2(k)+(-rho)*dx3(k)+0*g2(k));
    h3(k+1)=(0*dx1(k)+1*dx2(k)+(-rho)*dx3(k));
  endfor

  figure(1)
  plot(t2,dx3,'black')
  title("dx3 Euler")


  figure(2)
  plot(t2,h3,'red')
  title("h3")

  ##EJERCICIO 2

  ##VECTOR T
  for i=1:N
    t(i)= t0+ (i-1)*dt;
  endfor

  ##PLOT FIGURE 1
  figure (3)
  plot(t,gr)
  title("Funcion discreta dada por el registro")

  ##EJERCICIO 3
  H3=fft(h3);

  figure(4)
  stem(abs(H3(1:70)))

  ##EJERCICIO 4
  Gr=fft(gr);

  figure(5)
  stem(abs(Gr(1:N/2)))

  ##EJERCICIO 5
  for k=1:N
    Ycf(k)=H3(k)*Gr(k);
  endfor

  figure(6)
  stem(abs(Ycf(1:N/20)))

  ##EJERCICIO 6
  ycf=dt*ifft(Ycf);

  figure(7)
  plot(t,ycf,'g')

  ##EJERCICIO 7
  cn=dt*conv(h3,gr);

  figure(8)
  plot(t(1:N-1),cn(1:N-1),'g')


endfunction
