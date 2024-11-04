function recuperatorio2

%%===============DATOS=====================
g=load('registro-13-nov.txt','-ascii');
N= length(g)
dt=0.083776;
t=zeros(1,N);

for j=1:N
  t(j)=dt*(j-1);
end
figure (1)
plot (t,g,'r')
grid on
legend("Funcion g(tn)")
%==========NORMA CUADRATICA=================
Norma_2_g = norm(g)
%===========TDF====================
Tp=dt*N;
dw=(2*pi)/Tp

G_tdf= fft(g,N);

mod_G=abs(G_tdf);

Norma_2_G= norm(mod_G)
 %Crear un vector de frecuencias
frecuencias = dw * (0:(N-1));

  % Graficar el módulo de la TDF
  figure(2)
  #grafico reducido
  stem(abs(G_tdf(1:N/9)))
  title("Transformada de Fourier (modulo)")
%=============CONVOLUCION============
A1=1/3;
kc=8;
w=kc*dw;
a=0.3;
p=a*w;
for i=1:N
    h(i)= A1 * exp(-p*t(i)) * sin(w * t(i))/w;
end

yc=dt*conv(h,g);
figure(3)
plot(t,yc(1:N),'-b')
grid on
legend("Convolucion")

figure(4)
plot(t,h(1:N),'b')
grid on
legend("h")

figure (5)
plot(t,g,'r',t,yc(1:N),'b')
grid on
legend ("G(t), Yc")
end

