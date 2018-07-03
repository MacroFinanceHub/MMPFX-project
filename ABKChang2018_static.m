function [residual, g1, g2] = ABKChang2018_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 40, 1);

%
% Model equations
%

SDF__ = params(1)*(exp(y(1))-params(4)*exp((1+params(5))*y(9))/(1+params(5)))/(exp(y(1))-params(4)*exp((1+params(5))*y(9))/(1+params(5)));
Phi__ = params(11)*(exp(y(15)-(y(15)))-1)^2/2;
Phiprim__ = params(11)*(exp(y(15)-(y(15)))-1);
T19 = exp(y(1))-params(4)*exp((1+params(5))*y(9))/(1+params(5));
T62 = y(8)^(1-params(8)-params(9))*y(3)^params(8);
T75 = (exp(y(18))/params(8))^params(8);
T76 = y(14)*T75;
T80 = (exp(y(5))/params(9))^params(9);
T81 = T76*T80;
T83 = exp(y(9))/(1-params(8)-params(9));
T84 = T83^(1-params(8)-params(9));
T109 = exp(y(1))+(1+Phi__)*exp(y(15))+exp(y(2))+params(3)*exp(y(6))^2/2;
T115 = 1-params(6)*y(10)^2/2;
T152 = exp(y(4)+y(17))+exp(y(36));
T208 = sqrt(1+2*y(24)^2/params(14));
T209 = (-1)+T208;
T349 = exp(y(1))*(1-params(27))-params(4)*exp(y(9))^(1+params(5))/(1+params(5));
T361 = (T19*params(1)*exp(y(1))-exp(y(1))*params(1)*T19)/(T19*T19);
T411 = (-((-(exp(y(4)+y(17))*exp(y(13)+y(21))))/(T152*T152)));
T460 = (T19*params(1)*(-(params(4)*(1+params(5))*exp((1+params(5))*y(9))/(1+params(5))))-params(1)*T19*(-(params(4)*(1+params(5))*exp((1+params(5))*y(9))/(1+params(5)))))/(T19*T19);
T521 = (-(exp(y(13)+y(21))/T152));
lhs =1;
rhs =SDF__*(y(3)+params(2)*exp(y(4)))/(exp(y(4))+params(3)*exp(y(6)));
residual(1)= lhs-rhs;
lhs =1;
rhs =SDF__*y(7);
residual(2)= lhs-rhs;
lhs =y(8);
rhs =params(4)*exp(params(5)*y(9));
residual(3)= lhs-rhs;
lhs =y(11);
rhs =T62*exp(params(9)*y(13))/y(14);
residual(4)= lhs-rhs;
lhs =exp(y(12));
rhs =T81*T84;
residual(5)= lhs-rhs;
lhs =exp(y(13)+y(5));
rhs =exp(y(18))*y(3)*params(9)/params(8);
residual(6)= lhs-rhs;
lhs =y(8)*exp(y(9));
rhs =exp(y(18))*y(3)*(1-params(8)-params(9))/params(8);
residual(7)= lhs-rhs;
lhs =exp(y(12));
rhs =T109/T115;
residual(8)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(13)*params(10))*exp(y(16));
residual(9)= lhs-rhs;
lhs =y(10)*(1+y(10));
rhs =1/params(6)*(1+y(11)*params(7)-params(7))+(1+y(10))*SDF__*y(10);
residual(10)= lhs-rhs;
lhs =exp(y(15));
rhs =exp(y(18))-params(2)*exp(y(18));
residual(11)= lhs-rhs;
lhs =exp(y(4));
rhs =1+Phi__+exp(y(15)-(y(15)))*Phiprim__;
residual(12)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(6))+exp(y(17));
residual(13)= lhs-rhs;
lhs =T152;
rhs =exp(y(19))+exp(y(20))+exp(y(13)+y(21));
residual(14)= lhs-rhs;
lhs =exp(y(19));
rhs =T152*(params(12)+params(13))*y(37)-exp(y(20))*y(7)*params(12)-exp(y(13)+y(21))*params(12)*y(22);
residual(15)= lhs-rhs;
lhs =y(23);
rhs =exp(y(13)+y(21))/T152;
residual(16)= lhs-rhs;
lhs =y(25);
rhs =SDF__*(1-params(12)+params(12)*y(29))*(y(7)*(1-y(33))-y(22));
residual(17)= lhs-rhs;
lhs =y(26);
rhs =SDF__*(1-params(12)+params(12)*y(29))*(y(37)-y(7)*(1+y(32)));
residual(18)= lhs-rhs;
lhs =y(24);
rhs =y(25)/y(26);
residual(19)= lhs-rhs;
lhs =y(23);
rhs =1/y(24)*T209;
residual(20)= lhs-rhs;
lhs =y(27);
rhs =y(7)*SDF__*(1-params(12)+params(12)*y(29))*(1+y(31));
residual(21)= lhs-rhs;
lhs =y(28);
rhs =T152/exp(y(19));
residual(22)= lhs-rhs;
lhs =y(28);
rhs =y(27)/(y(30)-(y(26)+y(23)*y(25)));
residual(23)= lhs-rhs;
lhs =y(30);
rhs =params(15)*(1+params(14)*y(23)^2/2);
residual(24)= lhs-rhs;
lhs =y(29);
rhs =y(28)*y(30);
residual(25)= lhs-rhs;
lhs =exp(y(21));
rhs =(-(exp(y(2)-y(13))-exp(y(5))))+y(22)*exp(y(21))+exp(y(35))-y(22)*exp(y(35));
residual(26)= lhs-rhs;
lhs =y(7);
rhs =(1+y(34))/(1+y(10));
residual(27)= lhs-rhs;
lhs =y(37);
rhs =(y(3)+params(2)*exp(y(4)))/exp(y(4));
residual(28)= lhs-rhs;
lhs =y(31);
rhs =y(32)*y(28)+y(23)*y(33)*y(28);
residual(29)= lhs-rhs;
lhs =y(34);
rhs =(y(34))+y(10)*(1-params(16))*params(17)+params(16)*(y(34)-(y(34)))+params(29)*x(1);
residual(30)= lhs-rhs;
lhs =y(22);
rhs =(1-params(18))*(y(22))+y(22)*params(18)+params(32)*x(2);
residual(31)= lhs-rhs;
lhs =y(16);
rhs =(1-params(19))*(y(16))+y(16)*params(19)+params(31)*x(3);
residual(32)= lhs-rhs;
lhs =log(y(14));
rhs =(1-params(20))*log((y(14)))+log(y(14))*params(20)+params(30)*x(4);
residual(33)= lhs-rhs;
lhs =exp(y(36));
rhs =exp(y(13)+y(35));
residual(34)= lhs-rhs;
lhs =y(36);
rhs =(1-params(24))*(y(36))+y(36)*params(24)-params(25)*(y(13)-(y(13)))+params(28)*x(5);
residual(35)= lhs-rhs;
lhs =y(33);
rhs =params(21)*(log(y(17))-log((y(17))));
residual(36)= lhs-rhs;
lhs =y(32);
rhs =params(22);
residual(37)= lhs-rhs;
lhs =exp(y(38));
rhs =exp(y(12))-exp(y(13)+y(5));
residual(38)= lhs-rhs;
lhs =y(39);
rhs =log(T349);
residual(39)= lhs-rhs;
lhs =y(40);
rhs =y(39)+params(1)*y(40);
residual(40)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(40, 40);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-((y(3)+params(2)*exp(y(4)))*T361/(exp(y(4))+params(3)*exp(y(6)))));
  g1(1,3)=(-(SDF__/(exp(y(4))+params(3)*exp(y(6)))));
  g1(1,4)=(-(((exp(y(4))+params(3)*exp(y(6)))*SDF__*params(2)*exp(y(4))-exp(y(4))*SDF__*(y(3)+params(2)*exp(y(4))))/((exp(y(4))+params(3)*exp(y(6)))*(exp(y(4))+params(3)*exp(y(6))))));
  g1(1,6)=(-((-(SDF__*(y(3)+params(2)*exp(y(4)))*params(3)*exp(y(6))))/((exp(y(4))+params(3)*exp(y(6)))*(exp(y(4))+params(3)*exp(y(6))))));
  g1(1,9)=(-((y(3)+params(2)*exp(y(4)))*T460/(exp(y(4))+params(3)*exp(y(6)))));
  g1(2,1)=(-(y(7)*T361));
  g1(2,7)=(-SDF__);
  g1(2,9)=(-(y(7)*T460));
  g1(3,8)=1;
  g1(3,9)=(-(params(4)*params(5)*exp(params(5)*y(9))));
  g1(4,3)=(-(exp(params(9)*y(13))*y(8)^(1-params(8)-params(9))*getPowerDeriv(y(3),params(8),1)/y(14)));
  g1(4,8)=(-(exp(params(9)*y(13))*y(3)^params(8)*getPowerDeriv(y(8),1-params(8)-params(9),1)/y(14)));
  g1(4,11)=1;
  g1(4,13)=(-(T62*params(9)*exp(params(9)*y(13))/y(14)));
  g1(4,14)=(-((-(T62*exp(params(9)*y(13))))/(y(14)*y(14))));
  g1(5,5)=(-(T84*T76*exp(y(5))/params(9)*getPowerDeriv(exp(y(5))/params(9),params(9),1)));
  g1(5,9)=(-(T81*T83*getPowerDeriv(T83,1-params(8)-params(9),1)));
  g1(5,12)=exp(y(12));
  g1(5,14)=(-(T84*T75*T80));
  g1(5,18)=(-(T84*T80*y(14)*exp(y(18))/params(8)*getPowerDeriv(exp(y(18))/params(8),params(8),1)));
  g1(6,3)=(-(exp(y(18))*params(9)/params(8)));
  g1(6,5)=exp(y(13)+y(5));
  g1(6,13)=exp(y(13)+y(5));
  g1(6,18)=(-(exp(y(18))*y(3)*params(9)/params(8)));
  g1(7,3)=(-(exp(y(18))*(1-params(8)-params(9))/params(8)));
  g1(7,8)=exp(y(9));
  g1(7,9)=y(8)*exp(y(9));
  g1(7,18)=(-(exp(y(18))*y(3)*(1-params(8)-params(9))/params(8)));
  g1(8,1)=(-(exp(y(1))/T115));
  g1(8,2)=(-(exp(y(2))/T115));
  g1(8,6)=(-(params(3)*exp(y(6))*2*exp(y(6))/2/T115));
  g1(8,10)=(-((-(T109*(-(params(6)*2*y(10)/2))))/(T115*T115)));
  g1(8,12)=exp(y(12));
  g1(8,15)=(-((1+Phi__)*exp(y(15))/T115));
  g1(9,2)=exp(y(2));
  g1(9,13)=(-(exp(y(16))*params(10)*exp(y(13)*params(10))));
  g1(9,16)=(-(exp(y(13)*params(10))*exp(y(16))));
  g1(10,1)=(-((1+y(10))*y(10)*T361));
  g1(10,9)=(-((1+y(10))*y(10)*T460));
  g1(10,10)=y(10)+1+y(10)-(SDF__*y(10)+SDF__*(1+y(10)));
  g1(10,11)=(-(1/params(6)*params(7)));
  g1(11,15)=exp(y(15));
  g1(11,18)=(-(exp(y(18))-params(2)*exp(y(18))));
  g1(12,4)=exp(y(4));
  g1(13,6)=(-exp(y(6)));
  g1(13,17)=(-exp(y(17)));
  g1(13,18)=exp(y(18));
  g1(14,4)=exp(y(4)+y(17));
  g1(14,13)=(-exp(y(13)+y(21)));
  g1(14,17)=exp(y(4)+y(17));
  g1(14,19)=(-exp(y(19)));
  g1(14,20)=(-exp(y(20)));
  g1(14,21)=(-exp(y(13)+y(21)));
  g1(14,36)=exp(y(36));
  g1(15,4)=(-(exp(y(4)+y(17))*(params(12)+params(13))*y(37)));
  g1(15,7)=exp(y(20))*params(12);
  g1(15,13)=exp(y(13)+y(21))*params(12)*y(22);
  g1(15,17)=(-(exp(y(4)+y(17))*(params(12)+params(13))*y(37)));
  g1(15,19)=exp(y(19));
  g1(15,20)=exp(y(20))*y(7)*params(12);
  g1(15,21)=exp(y(13)+y(21))*params(12)*y(22);
  g1(15,22)=exp(y(13)+y(21))*params(12);
  g1(15,36)=(-(exp(y(36))*(params(12)+params(13))*y(37)));
  g1(15,37)=(-(T152*(params(12)+params(13))));
  g1(16,4)=T411;
  g1(16,13)=T521;
  g1(16,17)=T411;
  g1(16,21)=T521;
  g1(16,23)=1;
  g1(16,36)=(-((-(exp(y(36))*exp(y(13)+y(21))))/(T152*T152)));
  g1(17,1)=(-((y(7)*(1-y(33))-y(22))*(1-params(12)+params(12)*y(29))*T361));
  g1(17,7)=(-(SDF__*(1-params(12)+params(12)*y(29))*(1-y(33))));
  g1(17,9)=(-((y(7)*(1-y(33))-y(22))*(1-params(12)+params(12)*y(29))*T460));
  g1(17,22)=SDF__*(1-params(12)+params(12)*y(29));
  g1(17,25)=1;
  g1(17,29)=(-((y(7)*(1-y(33))-y(22))*SDF__*params(12)));
  g1(17,33)=(-(SDF__*(1-params(12)+params(12)*y(29))*(-y(7))));
  g1(18,1)=(-((y(37)-y(7)*(1+y(32)))*(1-params(12)+params(12)*y(29))*T361));
  g1(18,7)=(-(SDF__*(1-params(12)+params(12)*y(29))*(-(1+y(32)))));
  g1(18,9)=(-((y(37)-y(7)*(1+y(32)))*(1-params(12)+params(12)*y(29))*T460));
  g1(18,26)=1;
  g1(18,29)=(-((y(37)-y(7)*(1+y(32)))*SDF__*params(12)));
  g1(18,32)=(-(SDF__*(1-params(12)+params(12)*y(29))*(-y(7))));
  g1(18,37)=(-(SDF__*(1-params(12)+params(12)*y(29))));
  g1(19,24)=1;
  g1(19,25)=(-(1/y(26)));
  g1(19,26)=(-((-y(25))/(y(26)*y(26))));
  g1(20,23)=1;
  g1(20,24)=(-(T209*(-1)/(y(24)*y(24))+1/y(24)*2*2*y(24)/params(14)/(T208+T208)));
  g1(21,1)=(-(y(7)*(1+y(31))*(1-params(12)+params(12)*y(29))*T361));
  g1(21,7)=(-(SDF__*(1-params(12)+params(12)*y(29))*(1+y(31))));
  g1(21,9)=(-(y(7)*(1+y(31))*(1-params(12)+params(12)*y(29))*T460));
  g1(21,27)=1;
  g1(21,29)=(-(y(7)*(1+y(31))*SDF__*params(12)));
  g1(21,31)=(-(y(7)*SDF__*(1-params(12)+params(12)*y(29))));
  g1(22,4)=(-(exp(y(4)+y(17))/exp(y(19))));
  g1(22,17)=(-(exp(y(4)+y(17))/exp(y(19))));
  g1(22,19)=(-((-(T152*exp(y(19))))/(exp(y(19))*exp(y(19)))));
  g1(22,28)=1;
  g1(22,36)=(-(exp(y(36))/exp(y(19))));
  g1(23,23)=(-((-(y(27)*(-y(25))))/((y(30)-(y(26)+y(23)*y(25)))*(y(30)-(y(26)+y(23)*y(25))))));
  g1(23,25)=(-((-(y(27)*(-y(23))))/((y(30)-(y(26)+y(23)*y(25)))*(y(30)-(y(26)+y(23)*y(25))))));
  g1(23,26)=(-(y(27)/((y(30)-(y(26)+y(23)*y(25)))*(y(30)-(y(26)+y(23)*y(25))))));
  g1(23,27)=(-(1/(y(30)-(y(26)+y(23)*y(25)))));
  g1(23,28)=1;
  g1(23,30)=(-((-y(27))/((y(30)-(y(26)+y(23)*y(25)))*(y(30)-(y(26)+y(23)*y(25))))));
  g1(24,23)=(-(params(15)*params(14)*2*y(23)/2));
  g1(24,30)=1;
  g1(25,28)=(-y(30));
  g1(25,29)=1;
  g1(25,30)=(-y(28));
  g1(26,2)=exp(y(2)-y(13));
  g1(26,5)=(-exp(y(5)));
  g1(26,13)=(-exp(y(2)-y(13)));
  g1(26,21)=exp(y(21))-y(22)*exp(y(21));
  g1(26,22)=(-(exp(y(21))-exp(y(35))));
  g1(26,35)=(-(exp(y(35))-y(22)*exp(y(35))));
  g1(27,7)=1;
  g1(27,10)=(-((-(1+y(34)))/((1+y(10))*(1+y(10)))));
  g1(27,34)=(-(1/(1+y(10))));
  g1(28,3)=(-(1/exp(y(4))));
  g1(28,4)=(-((exp(y(4))*params(2)*exp(y(4))-exp(y(4))*(y(3)+params(2)*exp(y(4))))/(exp(y(4))*exp(y(4)))));
  g1(28,37)=1;
  g1(29,23)=(-(y(33)*y(28)));
  g1(29,28)=(-(y(32)+y(23)*y(33)));
  g1(29,31)=1;
  g1(29,32)=(-y(28));
  g1(29,33)=(-(y(23)*y(28)));
  g1(30,10)=(-((1-params(16))*params(17)));
  g1(31,22)=1-(params(18)+1-params(18));
  g1(32,16)=1-(params(19)+1-params(19));
  g1(33,14)=1/y(14)-((1-params(20))*1/(y(14))+params(20)*1/y(14));
  g1(34,13)=(-exp(y(13)+y(35)));
  g1(34,35)=(-exp(y(13)+y(35)));
  g1(34,36)=exp(y(36));
  g1(35,36)=1-(params(24)+1-params(24));
  g1(36,17)=(-(params(21)*(1/y(17)-1/(y(17)))));
  g1(36,33)=1;
  g1(37,32)=1;
  g1(38,5)=exp(y(13)+y(5));
  g1(38,12)=(-exp(y(12)));
  g1(38,13)=exp(y(13)+y(5));
  g1(38,38)=exp(y(38));
  g1(39,1)=(-(exp(y(1))*(1-params(27))/T349));
  g1(39,9)=(-((-(params(4)*exp(y(9))*getPowerDeriv(exp(y(9)),1+params(5),1)/(1+params(5))))/T349));
  g1(39,39)=1;
  g1(40,39)=(-1);
  g1(40,40)=1-params(1);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],40,1600);
end
end
