% exampleProjection.m plots an IGES CAD-object and project points on its surface

% Compile the c-files (if necessary)
% makeIGESmex;

% Load parameter data from IGES-file.
[ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab('example.igs');

% Project a set of points in a regular grid on the surface

% Projection data

R=[2000;-1300;0]; % Grid origin
normal=[2;8;-2.5]; % Direction of projection
normal=normal/norm(normal);
pdir=randn(3,1); % First (primary) direction of grid
pdir=pdir-dot(pdir,normal)*normal;
pdir=pdir/norm(pdir);
sdir=cross(pdir,normal); % Second (secondary) direction of grid
sdir=sdir/norm(sdir);

dp=10; % Grid spacing
nppos=20; % Number of positive grids in primary direction
npneg=40; % Number of negative grids in primary direction
nspos=25; % Number of positive grids in secondary direction
nsneg=35; % Number of negative grids in secondary direction

% Do the projection
[model,UV,srfind,srfDerivind,srfDer,numpoints]=projpartIGES(ParameterData,R,normal,pdir,sdir,dp,nppos,npneg,nspos,nsneg);

% Plots

% Plot the IGES object
figno=1;
plotIGES(ParameterData,1,figno,[],0);
for i=-npneg:nppos
    for j=-nsneg:nspos
        Pnt=R+i*dp*pdir+j*dp*sdir;
        plot3(Pnt(1),Pnt(2),Pnt(3),'.','Color',[0.1 0.5 0.4]);
    end
end

% Plot normal, starting at grid origin
lines=[R R+15*dp*normal];
plot3(lines(1,:),lines(2,:),lines(3,:),'r-');

% Plot pdir, starting at grid origin
lines=[R R+15*dp*pdir];
plot3(lines(1,:),lines(2,:),lines(3,:),'g-');

% Plot sdir, starting at grid origin
lines=[R R+15*dp*sdir];
plot3(lines(1,:),lines(2,:),lines(3,:),'b-');

% Plot projected points
plot3(model(1,srfind>0),model(2,srfind>0),model(3,srfind>0),'.','Color',[0.0 0.8 0.8]);

hold off;
