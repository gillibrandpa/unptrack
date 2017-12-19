function [pconc, xc, yc, area] = plot_particle_conc(gridFileName,pTime,Depth,pStage,units,varargin)

%
% Plot horizontal distributions of particle concentrations.
%
% Usage:    [xg,yg,pconc] = plot_particle_conc(gridFileName,pTime,Depth,pStage,plotYN,units);
%
%           pTime = time in hours of particle output file
%                   e.g. ptime = 1; file = ptrack_000001.dat
%           Depth = the limits of the layer to plot e.g. depth = [10 20] plots particle concentrations between 10 and 20 m.
%                   if depth has a single value, z, the surface layer from 0 - z m is calculated.
%                   if depth = 999, the depth-integrated concentration over the whole water column is calculated.
%           pStage: 0 = plot all stages, 1 = plot nauplii, 2 = plot copepodids only.
%           units: the units e.g. 'ug/L' for the plots (data is returned
%                   in kg/m3).
% and
%           varargin can include:
%                   plotYN: switch 'y' or 'n' to plot results [default is 'n'];
%                   axlim: axis limits [xmin xmax ymin ymax] defining the area to be plotted
%

% start time
tstring = num2str(pTime(1));

% load particle data
fileno = num2str(pTime(1) + 100000000);
fileno(1) = '0';
filename = ['pdnsty-',fileno,'.dat'];
disp(['Loading ',filename]);

% Load data
pxyzc = load(filename);

if size(pxyzc,1) == 0
    disp([filename, ' is empty. Stopping']);
    return
end

% load grid
xyz = load(gridFileName);
np = xyz(1,1);
ne = xyz(1,2);
nz = xyz(1,3);
xg = xyz(2:np+1,1);
yg = xyz(2:np+1,2);
h = xyz(2:np+1,3);
h(find(h == 0)) = NaN;
ele = xyz(np+2:end,:);
for ie = 1:ne
    xc(ie,1) = sum(xg(ele(ie,:)))/3;
    yc(ie,1) = sum(yg(ele(ie,:)))/3;
end

% calculate element areas
for k = 1: ne
    xn = xg(ele(k,:));
    yn = yg(ele(k,:));
    f1(1) = xn(1) - xn(2);
    f1(2) = yn(1) - yn(2);
    f2(1) = xn(1) - xn(3);
    f2(2) = yn(1) - yn(3);
    area(k,1) = abs(f1(1) * f2(2) - f1(2) * f2(1)) / 2.0;
end

% Derive particle elements and depth
if size(pxyzc,2) == 4,
    zc = unique(pxyzc(:,3));
    zp = pxyzc(:,3);
    ndep = length(unique(pxyzc(:,3)));
    pe = repmat([1:ne], ndep,1);
    pe = reshape(pe,ne*ndep,1);
elseif size(pxyzc,2) == 5,
    pe = unique(pxyzc(:,1));
    zp = pxyzc(:,4);
end

% Check the layer range
if size(Depth) == 1;
    Depth = [-Depth 0];
else
    Depth = [-Depth(2) -Depth(1)];
end
zrange = Depth(2) - Depth(1);

% Initialise the concentration vector
pconc = zeros(ne,1);

% Calculate depth-averaged values for specified range
for ie = 1:length(pe)
    element = pe(ie);
    id1 = pxyzc(:,1) == element;
    id2 = zp > Depth(1) & zp <= Depth(2);
    id = id1 & id2;
    if sum(id) > 0
        pconc(element,1) = sum(pxyzc(id,end)) / zrange;
    end
end

if nargin > 5
    plotYN = varargin{1};
else
    plotYN = 'n';
end
if strcmpi(plotYN,'n'); return; end

% Plot particle concentrations

% set figure properties
fh = findall(0,'type','figure');
if isempty(fh)
    fh = 1;
else
    fh = max(fh.Number);
end
figure(fh);
clf(fh);
set(gcf,'position',get(0,'Screensize'));
set(gcf,'papertype','a4');
orient portrait;
colormap('default');
axes;
hold on;

% Number of frames/images to plot
nFrame = length(pTime);

% If more than one frame, make movie
if nFrame > 1;
    OBJ = VideoWriter('pConcPlayBack.avi');
    dt = pTime(2) - pTime(1);
    OBJ.FrameRate = 7200 / dt;
    open(OBJ);
end

% Loop through files, plotting particles and capturing image
for ic = 1:nFrame
   
    % Delete contours from previous step
    if ic > 1
        delete(phan);
        %delete(hbar);
        % load next data file
        fileno = num2str(pTime(ic) + 100000000);
        fileno(1) = '0';
        filename = ['pdnsty-',fileno,'.dat'];
        disp(['Loading ',filename])
        pxyzc = load(filename);
        pe = unique(pxyzc(:,1));
        zp = pxyzc(:,4);
        pconc = zeros(ne,1);
        % Calculate depth-averaged values for specified range
        for ie = 1:length(pe)
            element = pe(ie);
            id1 = pxyzc(:,1) == element;
            id2 = zp > Depth(1) & zp <= Depth(2);
            id = id1 & id2;
            if sum(id) > 0
                pconc(element,1) = sum(pxyzc(id,end)) / zrange;
            end
        end
    end
        
    % Interpolate results from element centres. This version currently
    % interpolates onto a regular grid, for ease of plotting. Remove the
    % meshgrid command to interpolate onto original nodes [xg yg].
    F.Method = 'natural';
    F = scatteredInterpolant(xc,yc,pconc);
    if nargin > 6
        axlim = varargin{2};
        dx = (axlim(2) - axlim(1)) / 1000;
        dy = (axlim(4) - axlim(3)) / 1000;
        xint = [axlim(1):dx:axlim(2)];
        yint = [axlim(3):dy:axlim(4)];
    else
        dx = (max(xg) - min(xg)) / 1000;
        dy = (max(yg) - min(yg)) / 1000;
        xint = [min(xg):25:max(xg)];
        yint = [min(yg):25:max(yg)];
    end
    [x,y] = meshgrid(xint,yint);
    p2plot = F(x,y);
    xlim = [min(min(x)) max(max(x))];
    ylim = [min(min(y)) max(max(y))];

    if ic == 1
        Fh = scatteredInterpolant(xg,yg,h);
        hplot = F(x,y);
    end    

    % Modify units for plots
    if exist('units') == 1
        if strcmp(units,'ng/L')
            p2plot = p2plot * 1e9;
        elseif strcmp(units,'ug/L')
            p2plot = p2plot * 1e6;
        elseif strcmp(units,'mg/L')
            p2plot = p2plot * 1e3;
        end
    else
        units = 'kg m^{-3}';
    end

    % plot log value
    pconc_log = log10(p2plot);
    id = isfinite(pconc_log);
    clmin = floor(min(pconc_log(id)));
    clmax = ceil(max(pconc_log(id)));
    cl = [clmin:clmax]';
    cl = [-4:1]';
    if length(cl) < 3
        dcl = (clmax - clmin)/5;
        cl = [clmin:dcl:clmax]';
    end
    pconc_log(pconc_log < min(cl))=min(cl);
    pconc_log(hplot < min(h)) = min(cl);
    cllab = 10.^cl;
    clim = [min(cl) max(cl)];

%   plot concentrations

    %[cs, phan] = tricontf(xg,yg,ele,pconc_log);
    [cs, phan] = contourf(x,y,pconc_log,cl);
    hold on
    if ic == 1        
        set(gca,'fontsize',14,'linewidth',1, ...
            'position', [0.1300 0.1542 0.5945 0.7314]);
        axis equal
        set(gca,'xlim',xlim,'ylim',ylim);
        cmap = colormap;
        cmap(1,:) = [1 1 1];
        colormap(cmap);
        hold on
        %grid on
        xlabel('Distance (m)','fontsize',14);
        ylabel('Distance (m)','fontsize',14);
        hbar = colorbar('vertical');
        set(hbar,'fontsize',14,'location','EastOutside');
        hbarpos = get(hbar,'position');
        set(hbar,'fontsize',14,'position',[0.6573    0.2481    0.0166    0.4755]);
        %set(hbar,'xaxislocation','bottom');
        %set(hbar,'xtick',[-2 -1 0 1],'xticklabel',hbarlab);
        title(hbar,units,'fontsize',14);
        set(hbar,'ytick',cl(2:end));
        set(hbar,'yticklabel',num2str(cllab(2:end)));
        %set(gca,'xlim',[1000 5000],'ylim',[5000 9000]);
    end

        % plot bathymetry
        id1 = xg >= xlim(1) & xg <= xlim(2);
        id2 = yg >= ylim(1) & yg <= ylim(2);
        id = id1 & id2;
        hplot = h(id);
        hcont = [round(min(hplot)/10)+1:round(max(hplot)/10)-1]*10;
        if length(hcont) > 6
            hcont = hcont(1:2:end);
        end
        [cs2, han2] = tricontf(xg,yg,ele,h,hcont,'--');
        set(han2,'facecolor','none');
        caxis(clim);
        [cs3, han3] = tricontf(xg,yg,ele,h,[max(h) max(h)],'-');
        set(han3,'linewidth',2,'facecolor','none','edgecolor','k');
   
        % Set axes
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
        xtext = xlim(1) + 0.75 * (xlim(2) - xlim(1));
        ytext = ylim(1) + 0.02 * (ylim(2) - ylim(1));
        text(xtext,ytext,['Depth = ',num2str(abs(Depth(2))),' - ',num2str(abs(Depth(1))), ...
            ' m'], 'fontsize',14);
        hold on
    
    %caxis([min(cl) max(cl)]);
    caxis(clim);

    % add title
    plotdate = pTime(ic);
    plottime = pTime(ic)/60; punits = ' minutes';
    if plottime > 120
        plottime = plottime / 60;
        punits = ' hours';
    end
    title(['Predicted Concentrations after ',num2str(plottime),punits],'fontsize',16);

    % save plot to file
    %set(gcf, 'InvertHardCopy', 'off');
    %print('-dpng','-r600',['pconc_',tstring]);
    if nFrame > 1
        frame = getframe(gcf);
        writeVideo(OBJ,frame);
    end

end
if nFrame > 1
    close(OBJ);
end

end


