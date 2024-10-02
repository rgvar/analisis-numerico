

function euler_atras
    # implícito
	t0=0;
	y0=4;
	dt=0.01;
	NDt=1000;

	f = @(t,y) (t/2 - y/2);

	t=zeros(NDt,1);
	y=zeros(NDt,1);

	t(1)=t0;
	y(1)=y0;
	for m=1:NDt - 1
		k = dt * f(t(m), y(m));
		t(m + 1) = t(m) + dt;
		y(m + 1) = y(m) + k;
		% fprintf("m:%d ; t=%f ; y=%f ; k1=%f \n",m,t(m),y(m),k);
	end

	figure(1)
	plot(t,y(:,1),'b')
	grid on

endfunction
