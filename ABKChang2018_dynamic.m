function [residual, g1, g2, g3] = ABKChang2018_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(40, 1);
SDF__ = params(1)*(exp(y(12))-params(4)*exp((1+params(5))*y(20))/(1+params(5)))/(exp(y(52))-params(4)*exp((1+params(5))*y(56))/(1+params(5)));
Phi__ = params(11)*(exp(y(26)-(steady_state(15)))-1)^2/2;
Phiprim__ = params(11)*(exp(y(26)-(steady_state(15)))-1);
T20 = params(1)*(exp(y(12))-params(4)*exp((1+params(5))*y(20))/(1+params(5)));
T28 = exp(y(52))-params(4)*exp((1+params(5))*y(56))/(1+params(5));
T54 = exp(y(15))+params(3)*exp(y(17));
T73 = y(19)^(1-params(8)-params(9))*y(14)^params(8);
T86 = (exp(y(5))/params(8))^params(8);
T87 = y(25)*T86;
T91 = (exp(y(16))/params(9))^params(9);
T92 = T87*T91;
T94 = exp(y(20))/(1-params(8)-params(9));
T95 = T94^(1-params(8)-params(9));
T120 = exp(y(12))+(1+Phi__)*exp(y(26))+exp(y(13))+params(3)*exp(y(17))^2/2;
T126 = 1-params(6)*y(21)^2/2;
T217 = SDF__*(1-params(12)+params(12)*y(60));
T226 = y(55)*(1-y(44))-exp(y(59)-y(24))*y(33);
T246 = sqrt(1+2*y(35)^2/params(14));
T247 = (-1)+T246;
T398 = exp(y(12))*(1-params(27))-params(4)*exp(y(20))^(1+params(5))/(1+params(5));
T407 = params(1)*exp(y(12))/T28;
T432 = (-(T20*exp(y(52))))/(T28*T28);
T482 = (-((-(exp(y(15)+y(28))*exp(y(24)+y(32))))/((exp(y(15)+y(28))+exp(y(47)))*(exp(y(15)+y(28))+exp(y(47))))));
T527 = params(1)*(-(params(4)*(1+params(5))*exp((1+params(5))*y(20))/(1+params(5))))/T28;
T565 = (-(T20*(-(params(4)*(1+params(5))*exp((1+params(5))*y(56))/(1+params(5))))))/(T28*T28);
T618 = (-(exp(y(24)+y(32))/(exp(y(15)+y(28))+exp(y(47)))));
T643 = params(11)*exp(y(26)-(steady_state(15)))*2*(exp(y(26)-(steady_state(15)))-1)/2;
lhs =1;
rhs =SDF__*(y(53)+params(2)*exp(y(54)))/T54;
residual(1)= lhs-rhs;
lhs =1;
rhs =SDF__*y(55);
residual(2)= lhs-rhs;
lhs =y(19);
rhs =params(4)*exp(params(5)*y(20));
residual(3)= lhs-rhs;
lhs =y(22);
rhs =T73*exp(params(9)*y(24))/y(25);
residual(4)= lhs-rhs;
lhs =exp(y(23));
rhs =T92*T95;
residual(5)= lhs-rhs;
lhs =exp(y(24)+y(16));
rhs =exp(y(5))*y(14)*params(9)/params(8);
residual(6)= lhs-rhs;
lhs =y(19)*exp(y(20));
rhs =exp(y(5))*y(14)*(1-params(8)-params(9))/params(8);
residual(7)= lhs-rhs;
lhs =exp(y(23));
rhs =T120/T126;
residual(8)= lhs-rhs;
lhs =exp(y(13));
rhs =exp(y(24)*params(10))*exp(y(27));
residual(9)= lhs-rhs;
lhs =y(21)*(1+y(21));
rhs =1/params(6)*(1+y(22)*params(7)-params(7))+SDF__*y(57)*(1+y(57))*exp(y(58)-y(23));
residual(10)= lhs-rhs;
lhs =exp(y(26));
rhs =exp(y(29))-params(2)*exp(y(5));
residual(11)= lhs-rhs;
lhs =exp(y(15));
rhs =1+Phi__+exp(y(26)-(steady_state(15)))*Phiprim__;
residual(12)= lhs-rhs;
lhs =exp(y(29));
rhs =exp(y(17))+exp(y(28));
residual(13)= lhs-rhs;
lhs =exp(y(15)+y(28))+exp(y(47));
rhs =exp(y(30))+exp(y(31))+exp(y(24)+y(32));
residual(14)= lhs-rhs;
lhs =exp(y(30));
rhs =(params(12)+params(13))*y(48)*(exp(y(1)+y(4))+exp(y(11)))-params(12)*y(18)*exp(y(6))-params(12)*y(8)*exp(y(24)+y(7));
residual(15)= lhs-rhs;
lhs =y(34);
rhs =exp(y(24)+y(32))/(exp(y(15)+y(28))+exp(y(47)));
residual(16)= lhs-rhs;
lhs =y(36);
rhs =T217*T226;
residual(17)= lhs-rhs;
lhs =y(37);
rhs =T217*(y(61)-y(55)*(1+y(43)));
residual(18)= lhs-rhs;
lhs =y(35);
rhs =y(36)/y(37);
residual(19)= lhs-rhs;
lhs =y(34);
rhs =1/y(35)*T247;
residual(20)= lhs-rhs;
lhs =y(38);
rhs =y(55)*T217*(1+y(42));
residual(21)= lhs-rhs;
lhs =y(39);
rhs =(exp(y(15)+y(28))+exp(y(47)))/exp(y(30));
residual(22)= lhs-rhs;
lhs =y(39);
rhs =y(38)/(y(41)-(y(37)+y(34)*y(36)));
residual(23)= lhs-rhs;
lhs =y(41);
rhs =params(15)*(1+params(14)*y(34)^2/2);
residual(24)= lhs-rhs;
lhs =y(40);
rhs =y(39)*y(41);
residual(25)= lhs-rhs;
lhs =exp(y(32));
rhs =(-(exp(y(13)-y(24))-exp(y(16))))+y(8)*exp(y(7))+exp(y(46))-y(8)*exp(y(10));
residual(26)= lhs-rhs;
lhs =y(18);
rhs =(1+y(9))/(1+y(21));
residual(27)= lhs-rhs;
lhs =y(48);
rhs =(y(14)+params(2)*exp(y(15)))/exp(y(1));
residual(28)= lhs-rhs;
lhs =y(42);
rhs =y(43)*y(39)+y(34)*y(44)*y(39);
residual(29)= lhs-rhs;
lhs =y(45);
rhs =(steady_state(34))+y(21)*(1-params(16))*params(17)+params(16)*(y(9)-(steady_state(34)))+params(29)*x(it_, 1);
residual(30)= lhs-rhs;
lhs =y(33);
rhs =(1-params(18))*(steady_state(22))+y(8)*params(18)+params(32)*x(it_, 2);
residual(31)= lhs-rhs;
lhs =y(27);
rhs =(1-params(19))*(steady_state(16))+params(19)*y(3)+params(31)*x(it_, 3);
residual(32)= lhs-rhs;
lhs =log(y(25));
rhs =(1-params(20))*log((steady_state(14)))+params(20)*log(y(2))+params(30)*x(it_, 4);
residual(33)= lhs-rhs;
lhs =exp(y(47));
rhs =exp(y(24)+y(46));
residual(34)= lhs-rhs;
lhs =y(47);
rhs =(1-params(24))*(steady_state(36))+y(11)*params(24)-params(25)*(y(24)-(steady_state(13)))+params(28)*x(it_, 5);
residual(35)= lhs-rhs;
lhs =y(44);
rhs =params(21)*(log(y(28))-log((steady_state(17))));
residual(36)= lhs-rhs;
lhs =y(43);
rhs =params(22);
residual(37)= lhs-rhs;
lhs =exp(y(49));
rhs =exp(y(23))-exp(y(24)+y(16));
residual(38)= lhs-rhs;
lhs =y(50);
rhs =log(T398);
residual(39)= lhs-rhs;
lhs =y(51);
rhs =y(50)+params(1)*y(62);
residual(40)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(40, 67);

  %
  % Jacobian matrix
  %

  g1(1,12)=(-((y(53)+params(2)*exp(y(54)))*T407/T54));
  g1(1,52)=(-((y(53)+params(2)*exp(y(54)))*T432/T54));
  g1(1,53)=(-(SDF__/T54));
  g1(1,15)=(-((-(SDF__*(y(53)+params(2)*exp(y(54)))*exp(y(15))))/(T54*T54)));
  g1(1,54)=(-(SDF__*params(2)*exp(y(54))/T54));
  g1(1,17)=(-((-(SDF__*(y(53)+params(2)*exp(y(54)))*params(3)*exp(y(17))))/(T54*T54)));
  g1(1,20)=(-((y(53)+params(2)*exp(y(54)))*T527/T54));
  g1(1,56)=(-((y(53)+params(2)*exp(y(54)))*T565/T54));
  g1(2,12)=(-(y(55)*T407));
  g1(2,52)=(-(y(55)*T432));
  g1(2,55)=(-SDF__);
  g1(2,20)=(-(y(55)*T527));
  g1(2,56)=(-(y(55)*T565));
  g1(3,19)=1;
  g1(3,20)=(-(params(4)*params(5)*exp(params(5)*y(20))));
  g1(4,14)=(-(exp(params(9)*y(24))*y(19)^(1-params(8)-params(9))*getPowerDeriv(y(14),params(8),1)/y(25)));
  g1(4,19)=(-(exp(params(9)*y(24))*y(14)^params(8)*getPowerDeriv(y(19),1-params(8)-params(9),1)/y(25)));
  g1(4,22)=1;
  g1(4,24)=(-(T73*params(9)*exp(params(9)*y(24))/y(25)));
  g1(4,25)=(-((-(T73*exp(params(9)*y(24))))/(y(25)*y(25))));
  g1(5,16)=(-(T95*T87*exp(y(16))/params(9)*getPowerDeriv(exp(y(16))/params(9),params(9),1)));
  g1(5,20)=(-(T92*T94*getPowerDeriv(T94,1-params(8)-params(9),1)));
  g1(5,23)=exp(y(23));
  g1(5,25)=(-(T95*T86*T91));
  g1(5,5)=(-(T95*T91*y(25)*exp(y(5))/params(8)*getPowerDeriv(exp(y(5))/params(8),params(8),1)));
  g1(6,14)=(-(exp(y(5))*params(9)/params(8)));
  g1(6,16)=exp(y(24)+y(16));
  g1(6,24)=exp(y(24)+y(16));
  g1(6,5)=(-(exp(y(5))*y(14)*params(9)/params(8)));
  g1(7,14)=(-(exp(y(5))*(1-params(8)-params(9))/params(8)));
  g1(7,19)=exp(y(20));
  g1(7,20)=y(19)*exp(y(20));
  g1(7,5)=(-(exp(y(5))*y(14)*(1-params(8)-params(9))/params(8)));
  g1(8,12)=(-(exp(y(12))/T126));
  g1(8,13)=(-(exp(y(13))/T126));
  g1(8,17)=(-(params(3)*exp(y(17))*2*exp(y(17))/2/T126));
  g1(8,21)=(-((-(T120*(-(params(6)*2*y(21)/2))))/(T126*T126)));
  g1(8,23)=exp(y(23));
  g1(8,26)=(-(((1+Phi__)*exp(y(26))+exp(y(26))*T643)/T126));
  g1(9,13)=exp(y(13));
  g1(9,24)=(-(exp(y(27))*params(10)*exp(y(24)*params(10))));
  g1(9,27)=(-(exp(y(24)*params(10))*exp(y(27))));
  g1(10,12)=(-(exp(y(58)-y(23))*(1+y(57))*y(57)*T407));
  g1(10,52)=(-(exp(y(58)-y(23))*(1+y(57))*y(57)*T432));
  g1(10,20)=(-(exp(y(58)-y(23))*(1+y(57))*y(57)*T527));
  g1(10,56)=(-(exp(y(58)-y(23))*(1+y(57))*y(57)*T565));
  g1(10,21)=y(21)+1+y(21);
  g1(10,57)=(-(exp(y(58)-y(23))*(SDF__*y(57)+SDF__*(1+y(57)))));
  g1(10,22)=(-(1/params(6)*params(7)));
  g1(10,23)=(-(SDF__*y(57)*(1+y(57))*(-exp(y(58)-y(23)))));
  g1(10,58)=(-(SDF__*y(57)*(1+y(57))*exp(y(58)-y(23))));
  g1(11,26)=exp(y(26));
  g1(11,5)=params(2)*exp(y(5));
  g1(11,29)=(-exp(y(29)));
  g1(12,15)=exp(y(15));
  g1(12,26)=(-(T643+exp(y(26)-(steady_state(15)))*Phiprim__+exp(y(26)-(steady_state(15)))*params(11)*exp(y(26)-(steady_state(15)))));
  g1(13,17)=(-exp(y(17)));
  g1(13,28)=(-exp(y(28)));
  g1(13,29)=exp(y(29));
  g1(14,15)=exp(y(15)+y(28));
  g1(14,24)=(-exp(y(24)+y(32)));
  g1(14,28)=exp(y(15)+y(28));
  g1(14,30)=(-exp(y(30)));
  g1(14,31)=(-exp(y(31)));
  g1(14,32)=(-exp(y(24)+y(32)));
  g1(14,47)=exp(y(47));
  g1(15,1)=(-((params(12)+params(13))*y(48)*exp(y(1)+y(4))));
  g1(15,18)=params(12)*exp(y(6));
  g1(15,24)=params(12)*y(8)*exp(y(24)+y(7));
  g1(15,4)=(-((params(12)+params(13))*y(48)*exp(y(1)+y(4))));
  g1(15,30)=exp(y(30));
  g1(15,6)=params(12)*y(18)*exp(y(6));
  g1(15,7)=params(12)*y(8)*exp(y(24)+y(7));
  g1(15,8)=params(12)*exp(y(24)+y(7));
  g1(15,11)=(-((params(12)+params(13))*y(48)*exp(y(11))));
  g1(15,48)=(-((params(12)+params(13))*(exp(y(1)+y(4))+exp(y(11)))));
  g1(16,15)=T482;
  g1(16,24)=T618;
  g1(16,28)=T482;
  g1(16,32)=T618;
  g1(16,34)=1;
  g1(16,47)=(-((-(exp(y(47))*exp(y(24)+y(32))))/((exp(y(15)+y(28))+exp(y(47)))*(exp(y(15)+y(28))+exp(y(47))))));
  g1(17,12)=(-(T226*(1-params(12)+params(12)*y(60))*T407));
  g1(17,52)=(-(T226*(1-params(12)+params(12)*y(60))*T432));
  g1(17,55)=(-(T217*(1-y(44))));
  g1(17,20)=(-(T226*(1-params(12)+params(12)*y(60))*T527));
  g1(17,56)=(-(T226*(1-params(12)+params(12)*y(60))*T565));
  g1(17,24)=(-(T217*(-(y(33)*(-exp(y(59)-y(24)))))));
  g1(17,59)=(-(T217*(-(exp(y(59)-y(24))*y(33)))));
  g1(17,33)=(-(T217*(-exp(y(59)-y(24)))));
  g1(17,36)=1;
  g1(17,60)=(-(T226*SDF__*params(12)));
  g1(17,44)=(-(T217*(-y(55))));
  g1(18,12)=(-((y(61)-y(55)*(1+y(43)))*(1-params(12)+params(12)*y(60))*T407));
  g1(18,52)=(-((y(61)-y(55)*(1+y(43)))*(1-params(12)+params(12)*y(60))*T432));
  g1(18,55)=(-(T217*(-(1+y(43)))));
  g1(18,20)=(-((y(61)-y(55)*(1+y(43)))*(1-params(12)+params(12)*y(60))*T527));
  g1(18,56)=(-((y(61)-y(55)*(1+y(43)))*(1-params(12)+params(12)*y(60))*T565));
  g1(18,37)=1;
  g1(18,60)=(-((y(61)-y(55)*(1+y(43)))*SDF__*params(12)));
  g1(18,43)=(-(T217*(-y(55))));
  g1(18,61)=(-T217);
  g1(19,35)=1;
  g1(19,36)=(-(1/y(37)));
  g1(19,37)=(-((-y(36))/(y(37)*y(37))));
  g1(20,34)=1;
  g1(20,35)=(-(T247*(-1)/(y(35)*y(35))+1/y(35)*2*2*y(35)/params(14)/(T246+T246)));
  g1(21,12)=(-(y(55)*(1+y(42))*(1-params(12)+params(12)*y(60))*T407));
  g1(21,52)=(-(y(55)*(1+y(42))*(1-params(12)+params(12)*y(60))*T432));
  g1(21,55)=(-(T217*(1+y(42))));
  g1(21,20)=(-(y(55)*(1+y(42))*(1-params(12)+params(12)*y(60))*T527));
  g1(21,56)=(-(y(55)*(1+y(42))*(1-params(12)+params(12)*y(60))*T565));
  g1(21,38)=1;
  g1(21,60)=(-(y(55)*(1+y(42))*SDF__*params(12)));
  g1(21,42)=(-(y(55)*T217));
  g1(22,15)=(-(exp(y(15)+y(28))/exp(y(30))));
  g1(22,28)=(-(exp(y(15)+y(28))/exp(y(30))));
  g1(22,30)=(-((-((exp(y(15)+y(28))+exp(y(47)))*exp(y(30))))/(exp(y(30))*exp(y(30)))));
  g1(22,39)=1;
  g1(22,47)=(-(exp(y(47))/exp(y(30))));
  g1(23,34)=(-((-(y(38)*(-y(36))))/((y(41)-(y(37)+y(34)*y(36)))*(y(41)-(y(37)+y(34)*y(36))))));
  g1(23,36)=(-((-(y(38)*(-y(34))))/((y(41)-(y(37)+y(34)*y(36)))*(y(41)-(y(37)+y(34)*y(36))))));
  g1(23,37)=(-(y(38)/((y(41)-(y(37)+y(34)*y(36)))*(y(41)-(y(37)+y(34)*y(36))))));
  g1(23,38)=(-(1/(y(41)-(y(37)+y(34)*y(36)))));
  g1(23,39)=1;
  g1(23,41)=(-((-y(38))/((y(41)-(y(37)+y(34)*y(36)))*(y(41)-(y(37)+y(34)*y(36))))));
  g1(24,34)=(-(params(15)*params(14)*2*y(34)/2));
  g1(24,41)=1;
  g1(25,39)=(-y(41));
  g1(25,40)=1;
  g1(25,41)=(-y(39));
  g1(26,13)=exp(y(13)-y(24));
  g1(26,16)=(-exp(y(16)));
  g1(26,24)=(-exp(y(13)-y(24)));
  g1(26,7)=(-(y(8)*exp(y(7))));
  g1(26,32)=exp(y(32));
  g1(26,8)=(-(exp(y(7))-exp(y(10))));
  g1(26,10)=y(8)*exp(y(10));
  g1(26,46)=(-exp(y(46)));
  g1(27,18)=1;
  g1(27,21)=(-((-(1+y(9)))/((1+y(21))*(1+y(21)))));
  g1(27,9)=(-(1/(1+y(21))));
  g1(28,14)=(-(1/exp(y(1))));
  g1(28,1)=(-((-((y(14)+params(2)*exp(y(15)))*exp(y(1))))/(exp(y(1))*exp(y(1)))));
  g1(28,15)=(-(params(2)*exp(y(15))/exp(y(1))));
  g1(28,48)=1;
  g1(29,34)=(-(y(44)*y(39)));
  g1(29,39)=(-(y(43)+y(34)*y(44)));
  g1(29,42)=1;
  g1(29,43)=(-y(39));
  g1(29,44)=(-(y(34)*y(39)));
  g1(30,21)=(-((1-params(16))*params(17)));
  g1(30,9)=(-params(16));
  g1(30,45)=1;
  g1(30,63)=(-params(29));
  g1(31,8)=(-params(18));
  g1(31,33)=1;
  g1(31,64)=(-params(32));
  g1(32,3)=(-params(19));
  g1(32,27)=1;
  g1(32,65)=(-params(31));
  g1(33,2)=(-(params(20)*1/y(2)));
  g1(33,25)=1/y(25);
  g1(33,66)=(-params(30));
  g1(34,24)=(-exp(y(24)+y(46)));
  g1(34,46)=(-exp(y(24)+y(46)));
  g1(34,47)=exp(y(47));
  g1(35,24)=params(25);
  g1(35,11)=(-params(24));
  g1(35,47)=1;
  g1(35,67)=(-params(28));
  g1(36,28)=(-(params(21)*1/y(28)));
  g1(36,44)=1;
  g1(37,43)=1;
  g1(38,16)=exp(y(24)+y(16));
  g1(38,23)=(-exp(y(23)));
  g1(38,24)=exp(y(24)+y(16));
  g1(38,49)=exp(y(49));
  g1(39,12)=(-(exp(y(12))*(1-params(27))/T398));
  g1(39,20)=(-((-(params(4)*exp(y(20))*getPowerDeriv(exp(y(20)),1+params(5),1)/(1+params(5))))/T398));
  g1(39,50)=1;
  g1(40,50)=(-1);
  g1(40,51)=1;
  g1(40,62)=(-params(1));
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],40,4489);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],40,300763);
end
end
