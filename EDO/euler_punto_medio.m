

function euler_punto_medio

    t0=0;
	y0=4;
	dt=0.01;
	NDt=1000;

	f = @(t,y) (t/2 - y/2);

	t=zeros(NDt,1);
	y=zeros(NDt,1);

	t(1)=t0;
	y(1)=y0;
	for k=2:NDt - 1
		t(k+1) = t(k) + dt;
		y(k+1) = y(k-1) + 2*dt*f(t(k),y(k));
	end

	figure(1)
	plot(t,y(:,1),'b')
	grid on



endfunction
