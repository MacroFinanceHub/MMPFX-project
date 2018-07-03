% script for replicate analysis 
% Subject: Macro IV
% By Alex Carrasco, 2018

%% [I] Settings
clear all; clc; close all;
plotprop;

% Graphs
path_g       = [cd '\Graphs'];
imprime      = @(x) print( gcf, '-depsc2', [path_g filesep x]);
imprpdf      = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis   = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20, 'Box', 'On', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 22, 'Fontangle', 'normal');
label_x      = @(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20,'interpreter','latex');
label_y      = @(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 16,'interpreter','latex');
label_z      = @(x) zlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 16,'interpreter','latex');

global Wbase Wbase_unc oo_ M_ options_
flag_welfare = 0;
flag_irf     = 1-flag_welfare;

taylor       = [1.25 1.5 2];
param_ve     = [0 1/2 1 2 3];
param_omeg   = [0 0.05 0.1 0.2];

%% [II] First order analysis
if flag_irf 
dynare ABKChang2018 -Dcompute_optimal_policy=0 -Dcompute_2ndIRF=0  noclearall;
oo_order1_low   = oo_;

set_param_value('RR',0.4);
[oo_.dr,~,M_,~,oo_] = resol(0,M_,options_,oo_);
info = stoch_simul(var_list_); %get decision rules and moments
oo_order1_high  = oo_;

set_param_value('v_e',2);
[oo_.dr,~,M_,~,oo_] = resol(0,M_,options_,oo_);
info = stoch_simul(var_list_); %get decision rules and moments
oo_order1_FXrule_ = oo_;

set_param_value('v_e',0);
set_param_value('omegatauDf',0.1);
[oo_.dr,~,M_,~,oo_] = resol(0,M_,options_,oo_);
info = stoch_simul(var_list_); %get decision rules and moments
oo_order1_MPrule_ = oo_;

set_param_value('v_e',2);
[oo_.dr,~,M_,~,oo_] = resol(0,M_,options_,oo_);
info = stoch_simul(var_list_); %get decision rules and moments
oo_order1_FXMPrule_ = oo_;

choque = {'res_b' 'res_Rf'};
vars   = {'GDP' 'C' 'Yx' 'M' 'Df' 'rer' 'Q' 'N' 'pi' 'i' 'Rf' 'A' 'Yf' 'tauDf' 'F' 'K' 'Kb' 'U' 'D' 'B'};

nn= 25; 
rng=1:nn;
cc=1;

for jj=1:numel(choque)
    for ii=1:numel(vars)
        nom = [vars{ii} '_' choque{jj}];     
        figure(cc);
        hold on;           
        irf1 = oo_order1_low.irfs.(nom)(rng)';
        irf2 = oo_order1_high.irfs.(nom)(rng)';
        irf3 = oo_order1_FXrule_.irfs.(nom)(rng)';
        irf4 = oo_order1_MPrule_.irfs.(nom)(rng)';
        irf5 = oo_order1_FXMPrule_.irfs.(nom)(rng)';
        
        h=plot(rng,[irf2 irf3 irf4 irf5]*100); 
        set(h,PropName,PropValues(1:4,:)); 
        
        formataxis(gca);
        label_x('periods');
        plot(rng,zeros(1,numel(rng)),':r','linewidth',1);
        axis tight;   
        
        if ii==1;
            legend('Baseline','FX rule','MPrud rule','FX+MPrud rule');
            formatlegend('southeast');        
        end
        imprime(nom);
        imprpdf(nom);
        cc=cc+1;
    end
end

delete([path_g '\*.eps']);
close all;
end 

%% [III] Second Order Analysis
if flag_welfare
    W_u = nan(numel(param_ve),numel(param_omeg),numel(taylor));
    W_c = W_u;
    dynare ABKChang2018 -Dcompute_optimal_policy=0 -Dcompute_2ndIRF=1 noclearall;
    optims  = optimoptions('fsolve','Display','iter','MaxFunEvals',1e5,'TolFun',1e-10);
    options_.nocorr=1;
    options_.noprint=1;
    options_.verbosity=0;
%     options_.nofunctions=1;
    
    for k = 1:numel(taylor)
        for i=1:numel(param_ve)                       
            for j=1:numel(param_omeg)
                % getting the benchmark welfare  
                set_param_value('omegapi',taylor(k));
                set_param_value('v_e',param_ve(i));                                              
                set_param_value('omegatauDf',param_omeg(j));                       
                set_param_value('lambda_cost',0);
                
                info = stoch_simul([]); %get decision rules and moments
                
                [oo_.dr,~,M_,options_,oo_] = resol(0,M_,options_,oo_); %get decision rules
                initial_condition_states = repmat(oo_.dr.ys,1,M_.maximum_lag); %get steady state as initial condition
                shock_matrix = zeros(1,M_.exo_nbr); %create shock matrix with number of time periods in rows
                y_sim = simult_(initial_condition_states,oo_.dr,shock_matrix,options_.order); %simulate one period to get value               
                
                Wbase     = y_sim(strmatch('U',M_.endo_names,'exact'),2);
                Wbase_unc = oo_.mean(strmatch('U',M_.endo_names,'exact'));            
                
                set_param_value('omegapi',1.5);
                set_param_value('v_e',0);
                set_param_value('omegatauDf',0);                     
                [oo_.dr,~,M_,options_,oo_] = resol(0,M_,options_,oo_); %get decision rules

                lambda_cost_unc  = csolve('get_consumption_equivalent_unconditional_welfare',0,[],1e-8,1000);
                lambda_cost_cond = csolve('get_consumption_equivalent_conditional_welfare',lambda_cost_unc,[],1e-8,1000);

                W_u(i,j,k)=lambda_cost_unc;
                W_c(i,j,k)=lambda_cost_cond;    
            end
            fprintf('$\\mathbf{%1.2f}$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ \\\\ \n', [param_ve(i) -W_c(i,:,k)*100]);
        end
    end
end

