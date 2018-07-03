function [ys,check] = ABK_2016_steadystate(ys,exo)
% function [ys,check] = modeloOpen1_steadystate(ys,exo)
% computes the steady state for the NK_baseline.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values fpr the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 of not (allows to impos restriction on parameters)

global M_ 

% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end
% initialize indicator
check = 0;

%% Enter model equations here

%options=optimset(); % set options for numerical solver
%if gammma1<0 % parameter violates restriction; Preventing this cannot be implemented via prior restriction as it is a composite of different parameters and the valid prior region has unknown form
%    check=1; %set failure indicator
%    return; %return without updating steady states
%end

%% ============================================
%% II. Stationary equilibrium
%% ============================================
A = 1;
rer  = 1;
Rf = 1.005;

% II.[1] Direct solution
pi = 0;
Q = 1;
SDF = Beta;
R = 1/Beta;
mC= (eta-1)/eta;
tauDf = 0;
tauK  = omegatauK;

% II.[2] Banks
sf = 1-Beta*Rf;
[s,x,fopt]=HelpSS_ABK(Beta,Sigma,xi,sf,theta,Gamma,tauK,tauDf);          % solve optimal contract

phi = (Beta-Sigma)/(Sigma*(s+sf*x) + xi*(1+s));
Theta = theta*(1+Gamma/2*x^2);
psi = Theta*phi;
Z = R*(1+s) - lambda;
Omega = Beta*(1-Sigma+Sigma*psi);
mu = Omega*(Z+lambda - (1-tauK)*R);
mudf = Omega*((1-tauDf)*R-Rf);
muf = mudf/mu;
tauN = tauK*phi + tauDf*phi*x;     
v =   Omega*(1+tauN)*R;

% II.[3] Families and producers
Kh = s/chi;
% chiKh = chi*Kh^2/2;
w =  (mC*A/(Z^alphaK*rer^alphaM))^(1/(1-alphaK-alphaM));
K = (w^(1+zeta)/zeta0)^(1/zeta)*Z^(-1)*alphaK/(1-alphaK-alphaM);
Kb = K-Kh;
M = alphaM*Z*K/(alphaK*rer);
Y = (Z/(alphaK*mC))*K;
GDP  = Y - rer*M;
B  = RR*GDP;
Yx = ( alphaM*Z/(alphaK*rer) - (x*Kb+(x-1)*B)*(1-Rf)/K )*K;
Yf = Yx/(rer^varphi);
C_K = Z/(alphaK*mC) - (1-lambda) - Yx/K - chi*Kh^2/(2*K);
C   = C_K*K;
I   = (1-lambda)*K;
Df  = (Kb+B)*x/rer;
i   = R-1;

N = (Kb+B)/phi;
D = Kb + B- N-rer*Df;
F = B*rer;
L = (w/zeta0)^(1/zeta);
Rb = Z+lambda;

Phi = 0;
Phiprim=0;
aux1 = SDF*(Z+lambda*Q)/(Q+chi*Kh);
aux2 = SDF*R;
aux3 = mudf;
aux4 = mu;
aux5 = v;

instutil = log((1-lambda_cost)*C-zeta0*L^(1+zeta)/(1+zeta));
U        = (1-Beta)^(-1)*instutil;

%% Logaritmos
GDP = log(GDP);
Q = log(Q);
C=log(C);
L=log(L);
N=log(N);
D=log(D);
Df=log(Df);
rer = log(rer);
M = log(M);
Yx = log(Yx);
F  = log(F);
B  = log(B);
I  = log(I);
K  = log(K);
Kh = log(Kh);
Kb = log(Kb);
Y  = log(Y);
Yf = log(Yf);


%% end own model equations

for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'M_.params(' num2str(iter) ') = ' M_.param_names(iter,:) ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
