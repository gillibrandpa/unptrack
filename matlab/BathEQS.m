function results = BathEQS(gridFileName, EQS, time, depth)

%
% function to analyse results from ptrack model for bath treatment EQS
% compliance. This function assumes results are output hourly from the
% model.
%
% Usage: results = BathEQS(EQS, EQS_time, depth);
%
% where:    EQS is the Environmental Quality Standard to be met (ug/L)
%           time is the time (seconds) of the model output files
%           depth is the depth range (m) over which the concentrations are
%           averaged.
%
%           results = [time Cmax ExceedanceArea];
%

% Check time is a column vector
[nrow,ncol] = size(time);
if nrow == 1 && ncol > nrow
    time = time';
end
ntime = length(time);

% Initialise output arrays
ExceedanceArea = zeros(size(time));
Cmax = zeros(size(time));

% Loop through output files
for it = 1:ntime
    itime = time(it);

    % Read data 
    [pconc,xg,yg,area] = plot_particle_conc(gridFileName,itime,depth,1,'n');
    
    % Convert the concentrations from kg/m3 to ug/L
    pconc = pconc * 1e6;
   
    % Get maximum concentration
    Cmax(it,1) = max(max(pconc));
    
    % Calculate the exceedance area
    id = pconc > EQS;
    ExceedanceArea(it,1) = sum(area(id)) / 1e6;
       
    % Display output
    disp(['Time = ',num2str(itime/3600),' hours: Max Conc = ',num2str(Cmax(it),'%10.5f'), ...
        ' ug/L; Area > EQS = ',num2str(ExceedanceArea(it),'%10.5f'),' km2']); 
    
end

% Write output to file
fileout = ['BathEQS_Results_0-',num2str(abs(depth(1))),'m.dat'];
fid = fopen(fileout,'w');
results =[time Cmax ExceedanceArea];
fprintf(fid,'%s\n','Time (h), Maximum Concentration (ug/L), Area > EQS (km^2)');
fprintf(fid,'%d, %10.5f, %10.5f\n',results');
fclose(fid);

end
    