function [Rho,Beta,model_free,NLL] = ML_fitting (u)

% Num_Trials_per_Block = 35; % less than 56
%
% Initial value for parameters and matrix
% Total_Num_Trials = Num_Trials_per_Block*Num_Blocks;
data_choice = u(:,1);
data_mag = u (:,2);
data_prob = u(:,3);
sure_mag = u(:,5) ;


P_star_gamble = mean((data_prob .* data_mag - sure_mag)> 0);
% P_star_gamble model the P_gamble for neutral agent i.e. the frequency of
% choosing risky gamble by neutral agent.
P_gamble= mean(data_choice);
model_free = P_gamble;
phi0 = [1;1];
%phi0 = [1;0];
% phi0(1):"rho" for exponential function
% phi0(2): beta for softmax Function
utility_fun = @(phi,data_prob,data_mag,sure_mag)...
    data_prob.*(data_mag.^ phi(1))- sure_mag;
softmax_fun = @(phi,data_prob,data_mag,sure_mag) ...
    1./(1 + exp(-phi(2).*utility_fun(phi,data_prob,data_mag,sure_mag)));
LL_fun = @(phi,data_prob,data_mag,sure_mag)...
    (log(softmax_fun(phi,data_prob,data_mag,sure_mag)) .* data_choice)...
    + (log(1- softmax_fun(phi,data_prob,data_mag,sure_mag)) .* (1-data_choice)) ;
MyNLL = @(phi) -sum((LL_fun(phi,data_prob,data_mag,sure_mag)));
%% Test the Model
test_utility_fun = feval(utility_fun, phi0, data_prob,data_mag,sure_mag); %test utility_fun
test_softmax_fun = feval(softmax_fun, phi0, data_prob,data_mag,sure_mag); % test softmax_fun
test_LL_fun = feval (LL_fun , phi0,data_prob,data_mag,sure_mag); % test LL fun
test_myNLL = feval (MyNLL,phi0);
%%  Parameters Estimation with maximum likelihood estimation (MLE)
%         opts = optimset('fminsearch');
%         [PhiHatML , NLL] = fminsearch(MyNLL, phi0, opts);
A = []; b = []; Aeq = []; beq = []; lb = [0.80 0]; ub = [1.2 1000];nonlcon=[];
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[PhiHatML,NLL,exitflag,output,lambda,grad,hessian]=...
    fmincon(MyNLL,phi0,A,b,Aeq,beq,lb,ub,nonlcon,options);
Rho = round(PhiHatML(1),2);
Beta = PhiHatML(2);

