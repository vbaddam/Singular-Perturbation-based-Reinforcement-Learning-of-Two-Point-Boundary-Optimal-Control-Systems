yn = 1;
un = 1;
e = 1;

Q = 1;
R = 1;

Kinit  = 0;
Ki = Kinit;
N  = 80;           
MaxIteration = 12;
%T  = 0.000000001; 
T = 0.01;

y0 = 1; 

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

%tw = linspace(0, 0.01, 30);

for i=1:29
%while ran~=2

    

    [t,X] = ode45(@(t,x)system_ini(t,x,Ki,ifLearned,expl_noise_freq, test, e, Q, R),[(i-1)*T,i*T],X(end,:));

    Dyy=[Dyy;kron(X(end,1:yn),X(end,1:yn))-kron(X(1,1:yn),X(1,1:yn))];
    Iyy=[Iyy;X(end,yn+1:yn+yn^2)-X(1,yn+1:yn+yn^2)];
    Iyu=[Iyu;X(end,yn+yn^2+1:end)-X(1,yn+yn^2+1:end)];
    
    ran = rank([Iyy Iyu]);

    x_save=[x_save;X];
    t_save=[t_save;t];
end
% 

ifLearned = 1; 
test = 1;

P_old = rand(yn);
%P = P_old;
P = eye(yn)*10;    
it = 0;            
p_save = [];      
k_save = [];       

P_old = rand(yn);
P = eye(yn)*10;

it = 0;


while norm(P-P_old)>1e-12
    
    it = it+1;                        
    P_old = P;
    QK = Q+Ki'*R*Ki;                   
    

    
    Theta = [Dyy,-Iyy*kron(eye(yn),Ki')-Iyu]; 
    
    
    Xi = -Iyy*QK(:);                 
    pp = pinv(Theta)*Xi;

    


    P = reshape(pp(1:yn*yn), [yn, yn]); 
    P = (P + P')/2;
    
   
    
    BPv = pp(end-(yn*un-1):end);
    Ki = inv(R)*reshape(BPv,un,yn)/2;
    %K = reshape(pp(5:6), 1, 2);
    
    p_save = [p_save, P];   
    k_save = [k_save, Ki];
    
    disp(['K_', num2str(it), '=']);
    disp(Ki);
    

end
