function [s,x,fopt]=HelpSS_ABK(beta,sigma,xi,sf,theta,gamma,tauK,tauD)
%% Ecuentra de s y x para el modelo Aoki, Benigno y Kiyotaki (2016)
% ============
% Argumentos
% ============
% Como insumo necesita el valor de algunos parametros
% [1] beta: tasa de descuento subjetiva de familias
% [2] sigma: Probabilidad del banco a sobrevivir dentro de un período
% [3] xi: 'Seed' para los bancos que entran al mercado (porcentaje)
% [4] sf: Ventajas en costos de endeudarse externamente (1-beta*Rf)
% [5] theta: Porcentaje que el banco puede 'robar' si sólo se endeudaría domésticamente.
% [6] gamma: Controla el incremento en la habilidad de 'robar' del banco.
% [7] tauK: Impuestos al retorno del capital bancario.
% [8] tauD: Impuesto al endeudamiento externo.
% ============
% Output
% ============
% [1] s: Exceso de retorno del capital respecto al costo de deposito nacional
% [2] x: Porción del capital bancario financiado con depósitos externos.
% =======================================================
% By Alex Carrasco, 2018
% =======================================================

    function yy=findSS(X,beta,sigma,xi,sf,theta,gamma,tauK,tauD)
        s=X(1);
        x=X(2);
        eq1 =  (1-sigma)*(beta*(s+sf*x) + xi*(1+s))*(sigma*(s+sf*x) + xi*(1+s)) - ...
                        theta*(1+gamma*x^2/2)*(beta-sigma)*( sigma*(1-beta)*(s+sf*x)  + (1-sigma)*xi*(1+s) );
        eq2 =  (s-tauK)/(sf - tauD)*(-1 + sqrt(1+2/gamma*( (sf-tauD)/(s-tauK))^2));
%         eq2 =  (s-tauK)/(sf - tauD)*(-1 + sqrt(1+( (sf-tauD)/(s-tauK))^2));
        yy(1) = eq1;
        yy(2) = x-eq2;
    end

fobj = @(X) findSS(X,beta,sigma,xi,sf,theta,gamma,tauK,tauD);      % Declaring the function

X0 =[0.001 ,0.2 ];                                                       % Initial Guess

optims=optimoptions('fsolve','Display','off','MaxIter',1e15,'MaxFunEvals',2e15,'TolFun',1e-15,'TolX',1e-15);
[XOpt,fopt,~]=fsolve(fobj, X0, optims);                        % Solving non linear equations


s=XOpt(1);
x=XOpt(2);
end