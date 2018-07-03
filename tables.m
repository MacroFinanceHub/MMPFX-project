% Aoki, Benigno y Kiyotaki (2016) + Chang (2018)
% By Alex Carrasco 2018.
% Rio 2018

%% I. Parameters
load ABKChang2018_results.mat;

% III.[2] Reportando Calibracion (table_BLCal.tex)
fid = fopen('baseline_calibration.tex' ,'wt');
fprintf(fid,'\\begin{tabular}{c | c}\n');
fprintf(fid,'\\textsf{Parametro} & \\textsf{Valor}\n');
fprintf(fid,'\\\\ \n');
fprintf(fid,'\\hline \n');
for ii=1:numel(M_.params)
    fprintf(fid,['\\textsf{%s} &  %2.2f'], M_.param_names(ii,:), M_.params(ii));
    fprintf(fid,'\\\\ \n');
end
fprintf(fid, '\\end{tabular}  \n');
fclose(fid);
