function plotGBMesh(finMatFname,bodyMatFname,varargin)
%% plotGBMesh Plots GB Body Mesh in 3D
%
% Simple usage:
%
% plotGBMesh('SurfMesh-Fin-CoarseThin.mat','SurfMesh-Body-Coarse.mat');
% 
% Input:
%
% finMatFname  - Name of the input .mat file
% bodyMatFname - Name of the input .mat file
% (optional) 'SaveFileName' - Name of the target file to save
% 
% Output:
%
%
%
% written by Chen Chen 11-3-2015
%

%% Load .mat file
Fin = load(finMatFname);
Body = load(bodyMatFname);

%% Handle input
saveMatFname = [];
if ~isempty(varargin)
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'SaveFileName'
                saveMatFname = varargin{i+1};
        end
    end
end

% Check saveMatFname
if isempty(strfind(saveMatFname,'.mat'))
    saveMatFname = 'GB-Mesh.mat'; % default
end

%% Correct Fin Offset
finOffect = 81.1183;
tmpFinStartPoint = Body.Res.rangeZ(2) - finOffect;
tmpOffset = tmpFinStartPoint - Fin.Res.rangeZ(2);
Fin.Res.rangeZ = Fin.Res.rangeZ + tmpOffset;

for i = 1:numel(Fin.MeshVertAll)
    Fin.MeshVertAll{i}(:,3) = Fin.MeshVertAll{i}(:,3) + tmpOffset;
end

save('GB-Mesh-Fine.mat')
clear tmpOffset tmpFinStartPoint

%% Plot
totleNumOfPts = Body.Res.NumOfPts + Fin.Res.NumOfPts;
% Assemble Plot Vertices
plotX = zeros(totleNumOfPts,1,'double');
plotY = zeros(totleNumOfPts,1,'double');
plotZ = zeros(totleNumOfPts,1,'double');
tmpInd = 1;
for i = 1:numel(Body.MeshVertAll)
    % Body Verts
    tmpLen = length(Body.MeshVertAll{i}(:,1));
    
    plotX(tmpInd:tmpInd+tmpLen-1) = Body.MeshVertAll{i}(:,1);
    plotY(tmpInd:tmpInd+tmpLen-1) = Body.MeshVertAll{i}(:,2);
    plotZ(tmpInd:tmpInd+tmpLen-1) = Body.MeshVertAll{i}(:,3);
    
    tmpInd = tmpInd + tmpLen;
end
for i = 1:numel(Fin.MeshVertAll)
    % Fin Verts
    tmpLen = length(Fin.MeshVertAll{i}(:,1));
    
    plotX(tmpInd:tmpInd+tmpLen-1) = Fin.MeshVertAll{i}(:,1);
    plotY(tmpInd:tmpInd+tmpLen-1) = Fin.MeshVertAll{i}(:,2);
    plotZ(tmpInd:tmpInd+tmpLen-1) = Fin.MeshVertAll{i}(:,3);
    
    tmpInd = tmpInd + tmpLen;
end

figure(1);clf;
plot3(plotX,plotY,plotZ,'.','MarkerSize',16);grid on;axis image;