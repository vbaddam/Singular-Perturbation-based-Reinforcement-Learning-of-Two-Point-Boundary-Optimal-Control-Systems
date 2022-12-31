yn = 1;
un = 1;
e = 0.05;

Q = 1;
R = 1;

Kinit  = -10;  
Kf = Kinit;
N  = 80;           
MaxIteration = 12;
T  = 0.01;          

y0 = 0; 

expl_noise_freq = (rand(un,100)-.5)*100; 

Dyy=[];
Iyy=[];
Iyu=[];

X=[y0;kron(y0',y0')';kron(y0,zeros(un,1))]';



ifLearned = 0;
test = 0;

x_save=[];
t_save=[];

ran = 0;
v = 0;

for i=1:N
%while ran~=2

    [t,X] = ode45(@(t,x)system_final(t,x,Kf,ifLearned,expl_noise_freq, test, e),[(i-1)*T,i*T],X(end,:));

    Dyy=[Dyy;kron(X(end,1:yn),X(end,1:yn))-kron(X(1,1:yn),X(1,1:yn))];
    Iyy=[Iyy;X(end,yn+1:yn+yn^2)-X(1,yn+1:yn+yn^2)];
    Iyu=[Iyu;X(end,yn+yn^2+1:end)-X(1,yn+yn^2+1:end)];
    
    ran = rank([Iyy Iyu]);

    x_save=[x_save;X];
    t_save=[t_save;t];
end


ifLearned = 1; 
test = 1;

P_old = zeros(yn);
%P = P_old;
P = eye(yn)*10;    
it = 0;            
p_save = [];      
k_save = [];       

P_old = zeros(yn);
P = eye(yn)*10;

it = 0;


while  norm(P-P_old)>1e-8
    
    it = it+1;                        
    P_old = P;
    QK = Q+Kf'*R*Kf;                   
    

    
    Theta = [Dyy,-Iyy*kron(eye(yn),Kf')-Iyu]; 
    
    
    Xi = -Iyy*QK(:);                 
    pp = pinv(Theta)*Xi;

    


    P = reshape(pp(1:yn*yn), [yn, yn]); 
    P = (P + P')/2;
    
   
    
    BPv = pp(end-(yn*un-1):end);
    Kf = inv(R)*reshape(BPv,un,yn)/2;
    %K = reshape(pp(5:6), 1, 2);
    
    p_save = [p_save, P];   
    k_save = [k_save, Kf];
    
    disp(['K_', num2str(it), '=']);
    disp(Kf);
    

end

% disp(size(p_save))
% 
% figure(1)
% subplot(211)
% plot(0:length(p_save)-1,p_save,'o',0:length(p_save)-1,p_save,'Linewidth',2);
% 
% legend('||P_k-P^*||')
% xlabel('Number of iterations')
% 
% subplot(212)
% plot(0:length(k_save)-1,k_save,'^',0:length(k_save)-1,k_save,'Linewidth',2)
% 
% legend('||K_k-K^*||')
% xlabel('Number of iterations')
% 
% disp(t(end));
% 
% e = 0.1;
% % 
tspan1 = linspace(1, 0, 100);
[tt1,xfin]=ode45(@(t,x)system_final(t,x,Kf,ifLearned,expl_noise_freq, test, e),tspan1,1);
% % % % % 
% % % % % % Keep track of the post-learning trajectories
% % % % % t_final = [t_save;tt];
% % % % % x_final = [x_save;xx];
% % % % % 
% % % % % x_ep1 = x_final + 0.01;
% % % % % x_ep2 = x_final + 0.001;
% % % % % 
% % % % % figure(2)
%plot(tt1, xfin, 'Linewidth',2);
% % 

