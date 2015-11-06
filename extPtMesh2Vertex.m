function extPtMesh2Vertex(matFname,varargin)
%% extPtMesh2Vertex export mesh points to .Vertex file for CFD
%
% Simple usage:
%
% extPtMesh2Vertex('GB-Mesh-Fine-Corrected.mat','GB-v1.vertex');
% 
% Input:
%
% matFname - Name of the input .mat file (which has the Body and Fin mesh)
% (optional) 'SaveFileName' - Name of the target file to save
% 
% Output:
%
% MeshVertAll - Vertices of the point mesh indexed by Z cross-sections
%
%
% written by Chen Chen 10-5-2015
%

%% Load .mat file
load(matFname);

%% Handle input
saveMatFname = [];
for i = 1:2:nargin-1
    switch varargin{i}
        case 'SaveFileName'
            saveMatFname = varargin{i+1};
    end
end

% Check saveMatFname
if isempty(strfind(saveMatFname,'.vertex'))
    saveMatFname = 'GB-v1.vertex'; % default
end

%% Correct Fin Offset
tmpOffset = Body.Res.rangeZ(2) - 81.1183;
tmpOffset = tmpOffset - Fin.Res.rangeZ(2);
Fin.Res.rangeZ = Fin.Res.rangeZ + tmpOffset;

for i = 1:numel(Fin.MeshVertAll)
    Fin.MeshVertAll{i}(:,3) = Fin.MeshVertAll{i}(:,3) + tmpOffset;
end

save('GB-Mesh-Fine.mat')
clear tmpOffset

%% Create target file handle
% Create file handle
fid = fopen(saveMatFname, 'w');

%% Start Write file
% STEP 1 - Count total number of vertices
% Count total points
% NumOfPts = 0;
% NumOfBodyPts = 0;
% NumOfFinPts = 0;
% for i = 1:numel(Body.MeshVertAll)
%     NumOfBodyPts = NumOfBodyPts + size(Body.MeshVertAll{i},1);
% end
% for i = 1:numel(Fin.MeshVertAll)
%     NumOfFinPts = NumOfFinPts + size(Fin.MeshVertAll{i},1);
% end

NumOfBodyPts = Body.Res.NumOfPts;
NumOfFinPts = Fin.Res.NumOfPts;
NumOfPts = NumOfBodyPts + NumOfFinPts;
fprintf(fid,'%d\n',NumOfPts);
fprintf('Total points %d, writing file %s...\n',NumOfPts,saveMatFname);

% STEP 2 - Assemble Mesh Points
for i = 1:numel(Body.MeshVertAll)
    for j = 1:size(Body.MeshVertAll{i},1)
        fprintf(fid,'%.4f %.4f %.4f\n',...   
            Body.MeshVertAll{i}(j,1),... X
            Body.MeshVertAll{i}(j,2),... Y
            Body.MeshVertAll{i}(j,3)); % Z
    end
end
for i = 1:numel(Fin.MeshVertAll)
    for j = 1:size(Fin.MeshVertAll{i},1)
        fprintf(fid,'%.4f %.4f %.4f\n',...   
            Fin.MeshVertAll{i}(j,1),... X
            Fin.MeshVertAll{i}(j,2),... Y
            Fin.MeshVertAll{i}(j,3)); % Z
    end
end