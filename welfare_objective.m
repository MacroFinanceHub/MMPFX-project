function [outvalue]=welfare_objective(x_opt,x_opt_name)
% function [outvalue]=welfare_objective(x_opt,x_opt_name)

global oo_ options_ M_

%% set parameter for use in Dynare
for ii=1:size(x_opt_name,1)
    set_param_value(x_opt_name{ii,1},x_opt(ii));
end

if any(x_opt<cell2mat(x_opt_name(:,2))) || any(x_opt>cell2mat(x_opt_name(:,3))) %make sure parameters are inside their bounds
    outvalue=10e6+sum([x_opt].^2); %penalty function
    return
end

var_list_ = char('U');
info = stoch_simul(var_list_); %get decision rules and moments
if info(1) %filter out error code
    outvalue=1e5+sum([x_opt].^2);
    return;
end
outvalue=-oo_.mean(strmatch('U',var_list_,'exact')); %extract Welfare
% oo_.mean
end