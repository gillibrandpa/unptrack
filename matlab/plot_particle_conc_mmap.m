function [p2plot, x, y, hg] = plot_particle_conc_mmap(map,gridFileName,pTime,Depth,pStage,units,varargin)

%
% Plot horizontal distributions of particle concentrations.
%
% Usage:    [xg,yg,pconc] = plot_particle_conc_mmap(map,gridFileName,pTime,Depth,pStage,plotYN,units);
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

if nargin > 6
    plotYN = varargin{1};
else
    plotYN = 'n';
end
if strcmpi(plotYN,'n'); return; end

% Plot particle concentrations

% Draw map in m_map;
eval(map);

% Number of frames/images to plot
nFrame = length(pTime);

% If more than one frame, make movie
if nFrame > 1;
    OBJ = VideoWriter('pConcPlayBack.avi');
    dt = pTime(2) - pTime(1);
    OBJ.FrameRate = 7200 / dt;
    open(OBJ);
end

% Make a regular grid for plotting
if nargin > 7
    lllim = varargin{2};
    [axlim(1), axlim(3)] = WGS84toNGR(lllim(3),lllim(1));
    [axlim(2), axlim(4)] = WGS84toNGR(lllim(4),lllim(2));
    dx = (axlim(2) - axlim(1)) / 1000;
    dy = (axlim(4) - axlim(3)) / 1000;
    xint = [axlim(1):dx:axlim(2)];
    yint = [axlim(3):dy:axlim(4)];
else
    dx = (max(xg) - min(xg)) / 1000;
    dy = (max(yg) - min(yg)) / 1000;
    xint = [min(xg):dx:max(xg)];
    yint = [min(yg):dy:max(yg)];
    [lllim(1), lllim(3)] = NGRtoWGS84(min(xg),min(yg));    
    [lllim(2), lllim(4)] = NGRtoWGS84(max(xg),max(yg));    
end
[x,y] = meshgrid(xint,yint);

% Convert x and y to lat and long
dlon = (lllim(2) - lllim(1)) / 1000;
dlat = (lllim(4) - lllim(3)) / 1000;
lon = lllim(1):dlon:lllim(2);
lat = lllim(3):dlat:lllim(4);
[lon,lat] = meshgrid(lon,lat);

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
        
    % Modify units for plots
    if exist('units','var')
        if strcmp(units,'ng/L')
            pconc = pconc * 1e9;
        elseif strcmp(units,'ug/L')
            pconc = pconc * 1e6;
        elseif strcmp(units,'mg/L')
            pconc = pconc * 1e3;
        end
    else
        units = 'kg m^{-3}';
    end

    % Interpolate results from element centres. This version currently
    % interpolates onto a regular grid, for ease of plotting. Remove the
    % meshgrid command to interpolate onto original nodes [xg yg].
    F.Method = 'natural';
    F = scatteredInterpolant(xc,yc,pconc);
    p2plot = F(x,y);
    
    if ic == 1
        Fh = scatteredInterpolant(xg,yg,h);
        hg = Fh(x,y);
        hg(hg <= min(h)) = NaN;
    end

    % plot log value
    pconc_log = log10(p2plot);
    id = isfinite(pconc_log);
    clmin = floor(min(pconc_log(id)));
    clmax = ceil(max(pconc_log(id)));
    cl = [clmin:clmax]';
    cl = [-4:0]';
    if length(cl) < 3
        dcl = (clmax - clmin)/5;
        cl = [clmin:dcl:clmax]';
    end
    pconc_log(pconc_log < min(cl)) = NaN;
    pconc_log(isnan(hg)) = NaN;
    cllab = 10.^cl;
    clim = [min(cl) max(cl)];

    % Plot concentrations
    %[cs, phan] = m_contourf(lon,lat,pconc_log,cl);
    [cs, phan] = m_contourf(lon,lat,hg);
    colorbar
    return
    caxis(clim);
    cmap = colormap;
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    hbar = colorbar('vertical');
    set(hbar,'fontsize',14,'location','EastOutside');
    hbarpos = get(hbar,'position');
    set(hbar,'fontsize',14,'position',[0.8156 0.2825 0.0166 0.4755]);
    %set(hbar,'xaxislocation','bottom');
    set(hbar,'ytick',cl,'yticklabel',cllab);
    title(hbar,units,'fontsize',14);
    %set(hbar,'ytick',cl(2:end));
    %set(hbar,'yticklabel',num2str(cllab(2:end)));
    %set(gca,'xlim',[1000 5000],'ylim',[5000 9000]);

    % Add text
    xtext = lllim(1) + 0.75 * (lllim(2) - lllim(1));
    ytext = lllim(3) + 0.03 * (lllim(4) - lllim(3));
    m_text(xtext,ytext,['Depth = ',num2str(abs(Depth(2))),' - ',num2str(abs(Depth(1))), ...
           ' m'], 'fontsize',14);
    
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


