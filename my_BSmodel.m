function price=my_BSmodel(S0,K,r,q,sigma,T,NumT,OptionType,ExerciseType)
% Inputs
% S0: the initial stock price
% K: the strike price
% r: the risk-free interest rate
% q: the yield on the underlying asset
% sigma: volatility of the market
% T: time to maturity in years
% NumT: number of step in one year (must be a multiple of 12)
% OptionType: 'p' for put and 'c' for call
% ExerciseType: 'a' for American and 'e' for European
% Outputs
% price: the price of the option
%-----------------------------------------------------------------------
% Example:
% To value an European call option with 
% T=4 months, K=21, S0=20, r=0.1, sigma=0.3, yield=0
% use the command my_BSmodel(20,21,0.1,0,0.3,4/12,3600,'c','e')
%-----------------------------------------------------------------------
NT=NumT*T;  %total number of time interval
dt=T./NT;
u=exp(sigma.*sqrt(dt));  %the up movement
d=1./u;  %the down movement
a=exp(r*dt);
p=(a-d)./(u-d);  %the probability of an up movement
f=zeros(NT+1,NT+1);  %initialization of the values
%option prices at maturity
for i=0:NT  %fixed at the (NT+1)th column
    if (OptionType=='p')
        f(NT+1-i,NT+1)=max(K-S0.*(u.^i).*(d.^(NT-i)),0);
    else
        f(NT+1-i,NT+1)=max(S0.*(u.^i).*(d.^(NT-i))-K,0);
    end
end
%backwards induction process
for j=NT:-1:1
    for i=0:(j-1)
        B_Value=exp(-r.*dt).*(p.*f(j-i,j+1)+(1-p).*f(j-i+1,j+1));
        %American put option
        if ((ExerciseType=='a')&&(OptionType=='p'))
            E_Payoff=max(K-S0.*(u.^i).*(d.^(j-i-1)),0);
            f(j-i,j)=max(B_Value,E_Payoff);
        %American call option
        elseif ((ExerciseType=='a')&&(OptionType=='c'))
            E_Payoff=max(S0.*(u.^i).*(d.^(j-i-1))-K,0);
            f(j-i,j)=max(B_Value,E_Payoff);
        %European put and call options
        else
            f(j-i,j)=max(B_Value,0);
        end
    end
end
price=f(1,1);
end