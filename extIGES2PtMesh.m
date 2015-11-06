function MeshVertAll = extIGES2PtMesh(igsfile,MatFname)
%% EXTIGES2PTMESH extract 3D mesh points from any given IGES-file
%
% Simple usage:
%
% MeshVertAll = extIGES2PtMesh('GB.igs','TestMesh1.mat')
% 
% Input:
%
% igsfile - Name of the input IGES file
% (optional) MatFname - Name of the .mat file to save
% 
% Output:
%
% MeshVertAll - Vertices of the point mesh indexed by Z cross-sections
%
%
% written by Chen Chen 7-17-2015
%

if nargin < 2
    MatFname = 'SurfMesh-Fine.mat';
end

%% Global Control Parameters
% useHighResol = 0; % Decide whether or not to use high resolution mesh
isFinMesh = 0;    % = 1 if this is for fin mesh
switchAxis = 0;

deltaX = 0.25; % in mm
deltaY = 10;    % in mm
deltaZ = 25;   % in mm

%% Add MinBoundSuite Toolbox
% addpath('./MinBoundSuite');

%% Global Parameters
% Spatial Resolution, unit: percentage of range
% if isFinMesh
%     % Fin Mesh
%     if ~useHighResol
%         % Low Resolution
%         dX = 1 / 10;
%         dY = 1 / 10;
%         dZ = 1 / 20; % dZ is the volume axis
%     else
%         % High Resolution
%         dX = 1 / 25;
%         dY = 1 / 50;
%         dZ = 1 / 100; % dZ is the volume axis
%     end
% else
%     % Body Mesh
%     if ~useHighResol
%         % Low Resolution
%         dX = 1 / 10;
%         dY = 1 / 10;
%         dZ = 1 / 20; % dZ is the volume axis
%     else
%         % High Resolution
%         dX = 1 / 100;
%         dY = 1 / 200;
%         dZ = 1 / 500; % dZ is the volume axis
%     end
% end

%% Read and Extract Parameters from IGES file
if isFinMesh
    % Fin Mesh
    [igs.ParameterData,igs.EntityType,igs.numEntityType,...
        igs.unknownEntityType,igs.numunknownEntityType] = iges2matlab(igsfile);
    
    % Abstract Vertices
    [verts,vID] = extIGES(igs.ParameterData);
else
    % Body Mesh
    disp('Loading data...');
    load('GB-IGES-SurfaceVerts-Best');
end
%% Preprocess Data
% Get Data range
if isFinMesh || switchAxis
    tmp = verts(:,3);
    verts(:,3) = verts(:,1);
    verts(:,1) = tmp;
    clear tmp
end

rangeX = [min(verts(:,1)),max(verts(:,1))];
rangeY = [min(verts(:,2)),max(verts(:,2))];
rangeZ = [min(verts(:,3)),max(verts(:,3))];

% Get deltaX, deltaY, deltaZ => the actual minimal step and grid size 
% deltaX = dX * range(rangeX);
% deltaY = dY * range(rangeY);
% deltaZ = dZ * range(rangeZ);

dX = deltaX;
dY = deltaY;
dZ = deltaZ;


% Sort vertices data by Z axis ascending
% [~,tmp.ID] = sort(verts(:,3));
% verts = verts(tmp.ID,:);

%% IGES Mech Plot
figure(1);clf;
subplot(1,2,1);
plotIGES(igs.ParameterData,1,1);
view([90,90,45]);
hold(gca,'on');grid on;
hdlPlot.CroSectionAll = plot3(0,0,0,'c-*');...
    hdlPlot.CroSectionAll.MarkerSize = 2;...
    hdlPlot.CroSectionAll.LineWidth = 2;
hdlPlot.CroSection = plot3(0,0,0,'m-*');...
    hdlPlot.CroSection.MarkerSize = 2;...
    hdlPlot.CroSection.LineWidth = 2;

% Release memory
clear igs

%% Sample 2D Cross-section along the Z Axis
dispstat('','init')      %Initialization. Does not print anything.
iZ = rangeZ(1); % iZ = min(Z)
zVal = zeros(length(vID)-1,1,'double');
zId = zeros(length(vID)-1,1);
if isFinMesh
    zThresh = 1e-1;
else
    zThresh = 1e-1;
end
% minThresh = 1e-5;

% Plots
figure(1);subplot(1,2,2);
hdlPlot.Mesh = plot(0,0,'bo','markersize',4);hold(gca,'on');
hdlPlot.Verts = plot(0,0,'mo','markersize',4);
hdlPlot.ConvHull = plot(0,0,'c','LineWidth',8);axis image;
hold(gca,'off');

% Preallocate Memory for Result
MeshVertAll = cell(length(rangeZ(1):deltaZ:rangeZ(2))+1,1);
i = 1;
NumOfPts = 0;
stopFlag = 0;
while ~stopFlag
    if iZ > rangeZ(2)
        iZ = rangeZ(2);
        stopFlag = 1;
    end

    % Find nearby vertices
    sliceVerts = verts(abs(verts(:,3)-iZ)<=zThresh,:);

    % Solve for Convex Hull
    vertsHit = minConvHull(sliceVerts,deltaX,deltaY,iZ,isFinMesh,hdlPlot);
    
    % Push Result
    MeshVertAll{i} = vertsHit;
    i = i + 1;
    NumOfPts = NumOfPts + size(vertsHit,1);
    
    % Update Plot
    if (isFinMesh || switchAxis) && ~isempty(vertsHit)
        hdlPlot.CroSection.XData = vertsHit(:,3);
        hdlPlot.CroSection.YData = vertsHit(:,2);
        hdlPlot.CroSection.ZData = vertsHit(:,1);
        if iZ == rangeZ(1)
            hdlPlot.CroSectionAll.XData = vertsHit(:,3);
            hdlPlot.CroSectionAll.YData = vertsHit(:,2);
            hdlPlot.CroSectionAll.ZData = vertsHit(:,1);
        else
            hdlPlot.CroSectionAll.XData = ...
                [hdlPlot.CroSectionAll.XData,vertsHit(:,3)'];
            hdlPlot.CroSectionAll.YData = ...
                [hdlPlot.CroSectionAll.YData,vertsHit(:,2)'];
            hdlPlot.CroSectionAll.ZData = ...
                [hdlPlot.CroSectionAll.ZData,vertsHit(:,1)'];
        end
    else
        hdlPlot.CroSection.XData = vertsHit(:,1);
        hdlPlot.CroSection.YData = vertsHit(:,2);
        hdlPlot.CroSection.ZData = vertsHit(:,3);
        if iZ == rangeZ(1)
            hdlPlot.CroSectionAll.XData = vertsHit(:,1);
            hdlPlot.CroSectionAll.YData = vertsHit(:,2);
            hdlPlot.CroSectionAll.ZData = vertsHit(:,3);
        else
            hdlPlot.CroSectionAll.XData = ...
                [hdlPlot.CroSectionAll.XData,vertsHit(:,1)'];
            hdlPlot.CroSectionAll.YData = ...
                [hdlPlot.CroSectionAll.YData,vertsHit(:,2)'];
            hdlPlot.CroSectionAll.ZData = ...
                [hdlPlot.CroSectionAll.ZData,vertsHit(:,3)'];
        end
    end
    
%     pause;
    % Creat 2D Meshgrid and interpolate vertices in 2D
    iZ = iZ + deltaZ;
    dispstat(sprintf('Current Z: %3.4f, %.2f%% done!',iZ,100*(iZ-rangeZ(1))/range(rangeZ)));    
end

disp('Saving Data...');
Res.dX = dX;
Res.dY = dY;
Res.dZ = dZ;
Res.rangeX = rangeX;
Res.rangeY = rangeY;
Res.rangeZ = rangeZ;
Res.deltaX = deltaX;
Res.deltaY = deltaY;
Res.deltaZ = deltaZ;
Res.NumOfPts = NumOfPts;
fprintf('Total %d points generated, saving data...\n',NumOfPts);
save(MatFname,'MeshVertAll','Res');


%% Dump
%     for ivID = 1:length(vID)-1
% %         [zVal(ivID),zID(ivID)] = min(abs( verts(vID(ivID):vID(ivID+1),3) - iZ ));
%         dZverts = verts(vID(ivID):vID(ivID+1),3) - iZ;
%         [zVal,zId] = min(abs(dZverts));
%         % If found direct match
%         if ~zVal
%             sliceVerts = [sliceVerts;verts(abs(verts(:,3)-iZ)<=minThresh,:)];
%         else
%             % Couteur segment is far away => ignore
%             if zVal > zThresh
%                 continue;
%             end
%             zValP = min(dZverts(dZverts >= 0));
%             zValN = min(dZverts(dZverts <  0));
%             % LERP
%             if ~isempty(zValN) && ~isempty(zValP) && abs(zValP) <= zThresh .. .
%                     && abs(zValN) <= zThresh
%                 baseVertN = verts(abs(verts(:,3)-(zValN+iZ))<=minThresh,:);
%                 baseVertP = verts(abs(verts(:,3)-(zValP+iZ))<=minThresh,:);
%                 zPercent = -zValN / (zValP - zValN);              
%                 sliceVerts = [sliceVerts;...
%                     baseVertN + zPercent * (baseVertP - baseVertN)];
%             elseif isempty(zValP)
%                 sliceVerts = [sliceVerts;...
%                     verts(abs(verts(:,3)-(zValN+iZ))<=minThresh,:)];
%             else
%                 sliceVerts = [sliceVerts;...
%                     verts(abs(verts(:,3)-(zValP+iZ))<=minThresh,:)];
%             end
%         end
%     end


function vertsHit = minConvHull(verts,deltaX,deltaY,iZ,isFinMesh,varargin)
%% Local Parameters
minDiffZ = 1e-5;
if isempty(varargin)
    PlotSW = 0; % DEFAULT
else
    PlotSW = 1;
    hdlPlot = varargin{1};
end

% Get 2D verts range
xMin = min(verts(:,1));
xMax = max(verts(:,1));
yMin = min(verts(:,2));
yMax = max(verts(:,2));

% Generate Meshgrid
if ~isFinMesh
    % Body Mesh
    [X,Y] = meshgrid(xMin:deltaX:xMax, yMin:deltaY:yMax);
else
    % Fin Mesh
    [X,Y] = meshgrid(0, yMin:deltaY:yMax);
end

% Convert grid vectors to vertices
% meshVerts = [X(:),Y(:)];

% Get Convex Hull Vertices
[K, A] = convhull(verts(:,1),verts(:,2));

% Area of the polygon < 1 => Null convex hull
if A < 1
    % Sort by Z axis
    verts(:,3) = abs(verts(:,3)-iZ);
    [~,tmp.ID] = sort(verts(:,3)-iZ,'ascend');
    verts = verts(tmp.ID,:);
    
    % Get the closest set of vertices
    diffZ = abs(diff(verts(:,3)));
    ID = find(diffZ > minDiffZ,1);
    if isempty(ID)
        ID = size(verts,1);
    end
    verts = verts(1:ID,:);
    verts(:,3) = iZ;
    
    % Get range
    xMin = min(verts(:,1));
    xMax = max(verts(:,1));
    yMin = min(verts(:,2));
    yMax = max(verts(:,2));
    
    % Sample vertices
    if ID < 2
        vertsHit = [];
        return;
    elseif ID == 2
        % Line 
        if range(verts(:,1)) < deltaX || range(verts(:,2)) < deltaY
            vertsHit = verts;
            return;
        end
        % Generate Meshgrid
        if ~isFinMesh
            % Body Mesh
            [X,Y] = meshgrid(xMin:deltaX:xMax, yMin:deltaY:yMax);
        else
            % Fin Mesh
            [X,Y] = meshgrid(0, yMin:deltaY:yMax);
        end
        vertsHit = [X(:),Y(:)];
        vertsHit(:,3) = iZ;
        return;
    elseif ID > 2
        K = convhull(verts(:,1),verts(:,2));
        % Generate Meshgrid
        if ~isFinMesh
            % Body Mesh
            [X,Y] = meshgrid(xMin:deltaX:xMax, yMin:deltaY:yMax);
        else
            % Fin Mesh
            [X,Y] = meshgrid(0, yMin:deltaY:yMax);
        end
    end
end

% Test Polygon inhit
[in,on] = inpolygon(X(:),Y(:),verts(K,1),verts(K,2));

% Save Valid Vertices
if abs(yMin - yMax) < deltaY
    vertsHit = [X(in|on)',Y(in|on)',repmat(iZ,length(X(in|on)),1)];
else
    vertsHit = [X(in|on),Y(in|on)];
    vertsHit(:,3) = iZ;
end

if PlotSW   
    hdlPlot.Mesh.XData = X(:);
    hdlPlot.Mesh.YData = Y(:);
    
    hdlPlot.Verts.XData = X(in|on);
    hdlPlot.Verts.YData = Y(in|on); 
    
    hdlPlot.ConvHull.XData = verts(K,1);
    hdlPlot.ConvHull.YData = verts(K,2);
    
    drawnow;
%     pause;
end



function [verts,vID] = extIGES(ParameterData,fine_flag)
% EXTIGES extracts surfaces, curves and points from IGES-file
%
% Simple usage:
%
% extIGES(ParameterData)
%
% Ordinary usage:
%
% plotIGES(ParameterData,srf,fignr,subd,holdoff_flag)
%
% Input:
%
% ParameterData - Parameter data from IGES file. ParameterData
%                 is the output from IGES2MATLAB
% srf - Flag for surface plotting. (0,1,2)
%        0 (default), no surfaces are plotted,
%        1 the surface is plotted as triangular patches,
%        2 the surface is plotted as triangular mesh
% fignr - Figure number of the plot. 1 default
% subd - Nuber of subdivisions when plotting curves
%        subd is nubmer of subdivisions for each parameter when
%        plotting surfaces. 100 default
% holdoff_flag - Bolean value (1/0). If 1 then hold off the plot
%                when the plot is done. 1 default
% fine_flag - Bolean value (1/0). If 0 the surface will be rough
%             and if 1 the surface will be finer. 0 default
% plotCrvPnts - Flag for curve and point plotting (1/0)
%        0, no curves and points are plotted,
%        1  curves and points are plotted (default)
% srfClr - Color of surface
%
% Output:
%
% handlePlot - plothandle
%
% m-file can be downloaded at
% http://www.mathworks.com/matlabcentral/fileexchange/13253-iges-toolbox
%
% written by Per Bergström 2013-11-14
%

subd=100;

% Plot with fine mesh
if nargin < 2
    fine_flag = 1;
end

if isempty(ParameterData)
    error('Empty ParameterData');
elseif not(iscell(ParameterData))
    error('Invalid ParameterData. Must be a cell array!');
end

subd=subd+1;  % subd now number of points, not number of subintervals

siz = length(ParameterData);


verts = [];
vID = 1;
for i=1:siz
    
    [P,isSCP,isSup]=retSrfCrvPnt(2,ParameterData,1,i,subd,3);
    
    if ~and(isSCP,not(isSup)) && not(isSCP)
        
        [P,isSCP]=retSrfCrvPnt(3,ParameterData,1,i);
        
        if ~isSCP
            
            if fine_flag
                [P,isSCP,isSup,~]=retSrfCrvPnt(1,ParameterData,1,i,subd,1);
            else
                [P,isSCP,isSup,~] = retSrfCrvPnt(1,ParameterData,1,i,subd);
            end
            
            if and(isSCP,not(isSup))
                verts = [verts;P'];
                vID = [vID;size(verts,1)];
            end
        end
    end
end
