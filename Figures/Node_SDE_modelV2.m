function [All_X,All_eta, All_v]=Node_SDE_modelV2(V1,V2,X0,eta0,Parameters)

%% Parameters setting
alpha=Parameters.alpha;
lambda=Parameters.lambda;

n = Parameters.n; 
rho1 = Parameters.rho1; 
sigma = Parameters.sigma; 
tau = Parameters.tau;
theta = Parameters.theta; 
dt = Parameters.dt; 
total_t = Parameters.total_t;

t_seq = round(total_t/dt); 
gene = length(X0); 

All_X=zeros(t_seq,gene);
All_eta=zeros(t_seq,gene);

%% model solution

All_X(1,:)=X0;
All_eta(1,:)=eta0;
All_v(1,:,:)=V1;

ts=1/3*t_seq;
te=2/3*t_seq;
flag=randi([1, 4]);

for t=2:t_seq
    % stage-analysis
    if t<ts
        V=V1;
    elseif t>te
        V=V2;
    else
        tp=(t-ts)/(te-ts);
        V=Network_Evolution(V1,V2,flag,tp);
    end
 
    %% weight computing    
    temp_w=zeros(1,gene);
    for i=1:length(X0)
        tw=sum(alpha(i,:).*V(i,:).*X0);
        temp_w(i)=tw;
    end

    %% perturbation
    Z=randn(1,gene);

    %% model solution

    eta=eta0-1/tau*eta0*dt+sqrt(2/tau)*sigma*sqrt(dt)*Z;

    X=(X0+lambda.*F(rho1,temp_w,theta,n)*dt)./(1+exp(eta0 - sigma^2/2)*dt);

    %% Ê±¼äµü´ú
    
    All_X(t,:)=X;
    
    All_eta(t,:)=eta;

    All_v(t,:,:)=V;

    X0=X;
    
    eta0=eta;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
function f=F(rho1,w,theta,n)

f=rho1+(1-rho1)*(sqrt(w)./theta);

end

%% 
function new_v=Network_Evolution(V1,V2,flag,tp)

switch flag
    case 1
        new_v=V1+(V2-V1)*tp;
    case 2
        h=(tanh(5*(tp-1/2))+1)/2;
        new_v=V1+(V2-V1)*h;
    case 3
        a=3;
        h=(exp(a*tp)-1)/(exp(a)-1);
        new_v=V1+(V2-V1)*h;
    case 4
        a=3;
        h=1-exp(-a*tp);
        new_v=V1+(V2-V1)*h;
end

end
