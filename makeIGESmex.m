function makeIGESmex()
% Makefile
% run
% >> makeIGESmex
% in MATLAB to compile the source code in the IGES-toolbox


try
    mex -v nrbevalIGES.c
end

try
    mex -v nrbSrfRegularEvalIGES.c
end

try
    mex -v closestNrbLinePointIGES.c
end

try
    mex -v nrbCrvPlaneIntrsctIGES.c
end

try
    mex -v LScrvApp.c
end

try
    mex -v offsetNURBSsurface.c
end

try
    mex -v createDVGtree.c
end

try
    mex -v icpDVGpnt2pntLS.c
end

try
    mex -v icpDVGpnt2pntTu.c
end
