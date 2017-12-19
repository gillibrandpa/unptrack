function [pxysz, pass, varargout] = plot_particle_loc_mmap(map,pTime,pStage,plotYN,flowfile)

% function [pxysz, pass, varargout] = plot_particle_loc_mmap(map,pTime,pStage,plotYN,flowfile)
% plot particle locations on map of model domain.
% Usage:    pTime = time in hours of particle output file
%                   e.g. pTime = 1; file = ptrack_000001.dat
%           pStage: 0 = plot all Stages, 1 = plot nauplii,
%                  2 = plot copepodids only.
%           plotYN: switch 'y' or 'n' to plot results

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

% Draw map in m_map;
eval(map);

% Number of frames/images to plot
nFrame = length(pTime);

% If more than one frame, make movie
if nFrame > 1
    OBJ = VideoWriter('ParticlePlayBackmmap.avi');
    dt = pTime(2) - pTime(1);
    framerate = 7200 / dt;
    OBJ.FrameRate = framerate;
    open(OBJ);
end

% Loop through files, plotting particles and capturing image
for ic = 1:nFrame
    
    % Delete particles from previous step
    if ic > 1
        if exist('hvec','var')
            delete(hvec);
        end
        delete(phan1);
        if pStage == 0
            delete(phan0);
        %    delete(phan2);
        %    delete(phan3);
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
        npart = sum(id);
    else
        id = pass(:,3) > 0;
        npart = sum(id);
    end
    if npart > 1e5
        pxysz = pxysz(id,:);
        pass = pass(id,:);
        pxysz = pxysz(10:10:end,:);
        pass = pass(10:10:end,:);
    end

    % Convert to latitude and longitude
    %OSref = particle2os(pxysz(:,1),pxysz(:,2));
    %longlat = os2longlat(OSref);
    %lon = longlat(:,1) * 180 / pi;
    %lat = longlat(:,2) * 180 / pi;
    [lon, lat] = NGRtoWGS84(pxysz(:,1),pxysz(:,2));
    
    % Plot vectors
    if exist('flowfile','var')
        hvec = m_vec(1,flowfile.lon,flowfile.lat,flowfile.u(:,ic),flowfile.v(:,ic),'headangle',30);
    end
    
    % Plot stages
    if pStage > 0
       id = pass(:,3) == pStage;
       phan1 = m_plot(lon,lat,'r.','markersize',3);
    else
        id = pass(:,3) == 0;
        phan0 = m_plot(lon(id),lat(id),'.','markersize',3,'markerface',[0.5 0.5 0.5],'markeredge',[0.5 0.5 0.5]);
        id = pass(:,3) == 1;
        phan1 = m_plot(lon(id),lat(id),'.','markersize',3,'markerface','r','markeredge','r');
        hold on
        %id = pass(:,3) == 2;
        %phan2 = m_plot(lon(id),lat(id),'.','markersize',3,'markerface','y','markeredge','y');
        %id = pass(:,3) == 3;
        %phan3 = m_plot(lon(id),lat(id),'.','markersize',3,'markerface','g','markeredge','g');
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