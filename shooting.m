%Setting parameter values, initial conditions and algorithm parametrization
delta = 0.0625;
theta = .67;
beta = 4/(1-theta+4*(1-delta));
z0=1;
h=.31;
z=1;
k0=((beta^(-1)-1+delta)/((1-theta)*(z0*h)^theta))^(-1/theta); %initial capital level
max_iter=2000; %maximum number of iterations
conv=0.01; %convergence criterion (in percent of steady-state consumption)

%Steady state
k_ss = ((beta^(-1)-1+delta)/((1-theta)*(z*h)^theta))^(-1/theta)
c_ss = (z*h)^theta*((beta^(-1)-1+delta)/((1-theta)*(z*h)^theta))^(-(1-theta)/theta)-delta*((beta^(-1)-1+delta)/((1-theta)*(z*h)^theta))^(-1/theta);

%Setting interval for c0 
if k0<k_ss
	c_minmax=[0 c_ss];
elseif k0>k_ss
	c_minmax=[c_ss (1-delta)*k0+k0^(1-theta)*(z0*h)^theta];
else
	c0=c_ss;
	disp(['The economy starts from the steady state. c0=' num2str(c0)]);
end

for i=1:max_iter
	%Guess for c0 
	c0=0.5*(c_minmax(1,2)+c_minmax(1,1));
	%Calculating k and c at t=1
    kt=(1-delta)*k0+k0^(1-theta)*(z*h)^theta-c0;
    ct=(beta*((1-theta)*kt^(-theta)*(z*h)^theta+1-delta))*c0;
	%Initializing vectors in which we will store k and c trajectories
    k=[k0;kt];
	  c=[c0;ct];
    %Sub-loop generating trajectories of k and c for t>1 until they "bend"
    while (k0<k_ss & kt>=k(end-1) & ct>=c(end-1)) | ...
	  (k0>k_ss & kt<=k(end-1) & ct<=c(end-1))
        kt=(1-delta)*k(end)+k(end)^(1-theta)*(z*h)^theta-c(end);
		    ct=(beta*((1-theta)*(kt^(-theta))*(z*h)^theta+1-delta))*c(end);
        k=[k;kt];
        c=[c;ct];
    end
    %Checking if we are close enough to the steady state
    if 100*abs(c(end)/c_ss-1)<conv
        break;
    end
    %Verifying if guess of c0 is too big or too small and appropriate correction of c_minmax
	if c(end)>c(end-1)
		c_minmax(1,2)=c0;
	else
		c_minmax(1,1)=c0;
    end
end
%Checking if the number of iterations was sufficient
if i==max_iter
    disp('Maximum number of iterations has been reached. Convergence criterion is not met. Increase max_iter od double-check your code.');
else
    disp(['Solution: c0=' num2str(c0)]);
end
