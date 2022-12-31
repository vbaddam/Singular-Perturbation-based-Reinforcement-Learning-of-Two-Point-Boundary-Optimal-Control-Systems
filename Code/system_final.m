function dX = system_final(t,X,K,ifLearned,expl_noise_freq, test, e)


x = X(1);

if ~ifLearned   % See if learning is stopped
	%u = sum(sin(expl_noise_freq*t));
    u = sum(sin(expl_noise_freq*t));
else
	u = -K*x;    % Exploitation
end

if ~test
    dx = act_sys(x,u, t, e);
    dxx = kron(x',x')';
    dux = kron(x',u')';
    dX  = [dx;dxx;dux];

else
    dX = act_sys(x,u, t, e);
    
end 

end


function dx = act_sys(x,u, t, e)
%% Actual dynamics of the desel engine. 
%  This is the system you can customize.

A = -(1+0.2*t);
%A = -1.2;
B = 1+0.2*t;

dx = 1/e*(A*x+B*u);
end

function dy = ep_sys(x, u, e)
    
    A = -1.2;
    
    B = 1.2;
    
    dy = -1/e*(-A*x-B*u);
end