function [pxysz, pass, varargout] = plot_particle_loc(pTime,pStage,plotYN,varargin)

% function [pxysz, pass] = plot_particle_loc(pTime,pStage,plotYN,varargin)
% plot particle locations on bathymetry map of model domain.
% Usage:    pTime = time in hours of particle output file
%                   e.g. pTime = 1; file = ptrack_000001.dat
%           pStage: 0 = plot all Stages, 1 = plot nauplii,
%                  2 = plot copepodids only.
%           plotYN: switch 'y' or 'n' to plot results
%           varargin can include:
%               1. gridFileName, containing grid coordinate [default =
%               'grid_coordinates.dat']
%               2. axis limits [xmin xmax ymin ymax];

% Plot time
tstring = num2str(pTime(1));

% load particle data
fileno = num2str(pTime(1) + 100000000);
fileno(1) = '0';
filename = ['ptrack-',fileno,'.dat'];
disp(['Loading ',filename])

% Load data
data = load(filename);
pxysz = data(:,1:4);
pass = data(:,5:7);

if plotYN == 'n'; return; end

% Plot particle positions

% load grid
if nargin <= 3
    gridFileName = 'grid_coordinates.dat';
else
    gridFileName = varargin{1};
end
xyz = load(gridFileName);
np = xyz(1,1);
ne = xyz(1,2);
nz = xyz(1,3);
x = xyz(2:np+1,1);
y = xyz(2:np+1,2);
h = xyz(2:np+1,3);
h(find(h == 0)) = NaN;
ele = xyz(np+2:end,:);

% set figure properties
figure(1);
clf(1);
set(gcf,'position',get(0,'ScreenSize'),'renderer','zbuffer');
set(gcf,'papertype','a4');
orient portrait;
colormap('default');

% plot bathymetry if not constant
hmin = min(h); hmax = max(h);
if hmax > hmin
    [CS,fhan]=tricontf(x,y,ele,h);
    set(fhan,'edgecolor','none');
    %cmap = colormap;
    %cmap = flipud(cmap(1:20,:));
    %colormap(cmap);
    hold on
end

% Number of frames/images to plot
nFrame = length(pTime);

% If more than one frame, make movie
if nFrame > 1
    OBJ = VideoWriter('ParticlePlayBack.avi');
    OBJ.FrameRate = 2;
    open(OBJ);
end

% Loop through files, plotting particles and capturing image
for ic = 1:nFrame
    
    % Delete particles from previous step
    if ic > 1
        delete(phan1);
        if pStage == 0
            delete(phan0);
            delete(phan2);
            delete(phan3);
        end        
        % load next data file
        fileno = num2str(pTime(ic) + 100000000);
        fileno(1) = '0';
        filename = ['ptrack-',fileno,'.dat'];
        disp(['Loading ',filename])
        data = load(filename);
        pxysz = data(:,1:4);
        pass = data(:,5:7);
    end
    
    % if large number of particles, do not plot all
    if pStage > 0
        id = pass(:,3) == pStage;
        pxysz = pxysz(id,:);
        pass = pass(id,:);
        npart = sum(id);
    else
        npart = size(pxysz,1);
    end
    if npart > 1e5
        pxysz = pxysz(10:10:end,:);
        pass = pass(10:10:end,:);
    end
    
    % Plot stages
    if pStage > 0
        phan1 = plot(pxysz(:,1),pxysz(:,2),'r.','markersize',3);
    else
        id = pass(:,3) == 0;
        phan0 = plot(pxysz(id,1),pxysz(id,2),'.','markersize',2,'markerface','k','markeredge','k');
        id = pass(:,3) == 1;
        phan1 = plot(pxysz(id,1),pxysz(id,2),'.','markersize',2,'markerface','r','markeredge','r');
        hold on
        id = pass(:,3) == 2;
        phan2 = plot(pxysz(id,1),pxysz(id,2),'.','markersize',2,'markerface','y','markeredge','y');
        id = pass(:,3) == 3;
        phan3 = plot(pxysz(id,1),pxysz(id,2),'.','markersize',2,'markerface','g','markeredge','g');
    end
    
    % Set axis properties
    if ic == 1
        if nargin > 4
            xylim = varargin{2};
            xlim = xylim(1:2);
            ylim = xylim(3:4);
        else
            xlim = [min(x) max(x)];
            ylim = [min(y) max(y)];
        end
        axis equal
        set(gca,'xlim',xlim,'ylim',ylim,'linewidth',1,'fontsize',14);
    end

    % add date string
    title(['Predicted particle locations (Time = ',num2str(pTime(ic)/3600),' hrs)'],'fontsize',16);
    pause(0.2);

    % Save plot to file
    %set(gcf, 'InvertHardCopy', 'off');
    %saves(gcf,['images\p_',fileno,'.png']); 
    if ic > 1
        frame = getframe(gcf);
        writeVideo(OBJ,frame);
    end
    
    % Calculate some stats
    if pStage > 0
        id = pass(:,3) == pStage;
    else
        id = pass(:,3) > 0;
    end    
    pstats.xmean(ic,1) = mean(pxysz(id,1));
    pstats.xstd(ic,1) = std(pxysz(id,1));
    pstats.xvar(ic,1) = var(pxysz(id,1));
    pstats.ymean(ic,1) = mean(pxysz(id,2));
    pstats.ystd(ic,1) = std(pxysz(id,2));
    pstats.yvar(ic,1) = var(pxysz(id,2));
    pstats.zmean(ic,1) = mean(pxysz(id,4));
    pstats.zstd(ic,1) = std(pxysz(id,4));
    pstats.zvar(ic,1) = var(pxysz(id,4));
    Nid = sum(id);
    pstats.pvar(ic,1) = sum((pxysz(id,1)-pstats.xmean(ic)).^2+(pxysz(id,2)-pstats.ymean(ic)).^2)/Nid;
    pstats.pstd = sqrt(pstats.pvar);
    
    
end
if nFrame > 1
    close(OBJ);
end
varargout{1} = pstats;
end