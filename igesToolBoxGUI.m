function varargout = igesToolBoxGUI(varargin)
% IGESTOOLBOXGUI M-file for igesToolBoxGUI.fig
%      IGESTOOLBOXGUI, by itself, creates a new IGESTOOLBOXGUI or raises the existing
%      singleton*.
%
%      H = IGESTOOLBOXGUI returns the handle to a new IGESTOOLBOXGUI or the handle to
%      the existing singleton*.
%
%      IGESTOOLBOXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IGESTOOLBOXGUI.M with the given input arguments.
%
%      IGESTOOLBOXGUI('Property','Value',...) creates a new IGESTOOLBOXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before igesToolBoxGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to igesToolBoxGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help igesToolBoxGUI

% Last Modified by GUIDE v2.5 26-Jan-2012 12:23:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @igesToolBoxGUI_OpeningFcn, ...
    'gui_OutputFcn',  @igesToolBoxGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles)
% If the plotdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.


axes(handles.axesCAD);
plot3(0,0,0,'k.','MarkerSize',1e-3);

handles.plotdata.filename = '';
set(handles.textFile,'String','');

handles.plotdata.igesfile = false;
handles.plotdata.stlfile = false;

handles.plotdata.az = 0;
handles.plotdata.el = 0;
handles.plotdata.rot = 0;

handles.plotdata.xl = [0,1];
handles.plotdata.yl = [0,1];
handles.plotdata.zl = [0,1];

handles.plotdata.dx = 1;
handles.plotdata.dy = 1;
handles.plotdata.dz = 1;

handles.plotdata.hxyz=zeros(3,1);
handles.plotdata.hTextX=0;
handles.plotdata.hTextY=0;
handles.plotdata.hTextZ=0;

handles.plotdata.saveas=false;

view(handles.plotdata.az,handles.plotdata.el);

xlim(handles.plotdata.xl);
ylim(handles.plotdata.yl);
zlim(handles.plotdata.zl);

set(handles.checkboxXYZdir,'Value',0);
set(handles.checkboxAxis,'Value',1);
set(handles.checkboxTransp,'Value',0);

set(handles.editAz,'String', num2str(handles.plotdata.az));
set(handles.sliderAz,'Value',handles.plotdata.az);

set(handles.editEl,'String', num2str(handles.plotdata.el));
set(handles.sliderEl,'Value',handles.plotdata.el);

set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);


% Update handles structure

guidata(handles.figure1, handles);



% --- Executes just before igesToolBoxGUI is made visible.
function igesToolBoxGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to igesToolBoxGUI (see VARARGIN)

% Choose default command line output for igesToolBoxGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initialize_gui(hObject, handles);

% UIWAIT makes igesToolBoxGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = igesToolBoxGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliderEl_Callback(hObject, eventdata, handles)
% hObject    handle to sliderEl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.plotdata.el=get(hObject,'Value');

handles.plotdata.rot = 0;
set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);

set(handles.editEl,'String', num2str(handles.plotdata.el));

axes(handles.axesCAD);
view(handles.plotdata.az,handles.plotdata.el);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderEl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderEl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderAz_Callback(hObject, eventdata, handles)
% hObject    handle to sliderAz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


handles.plotdata.az=get(hObject,'Value');

handles.plotdata.rot = 0;
set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);

set(handles.editAz,'String', num2str(handles.plotdata.az));

axes(handles.axesCAD);
view(handles.plotdata.az,handles.plotdata.el);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function sliderAz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderAz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function editAz_Callback(hObject, eventdata, handles)
% hObject    handle to editAz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editAz as text
%        str2double(get(hObject,'String')) returns contents of editAz as a double

az=str2double(get(hObject, 'String'));

if isnan(az)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    az=0;
end

if or(az<-180,az>=180)
    az=az-floor((az+180)/360)*360;
    set(hObject, 'String', num2str(az));
end

set(handles.sliderAz,'Value',az);

handles.plotdata.az = az;

handles.plotdata.rot = 0;
set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);

axes(handles.axesCAD);
view(handles.plotdata.az,handles.plotdata.el);

% Update handles structure

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editAz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editAz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editEl_Callback(hObject, eventdata, handles)
% hObject    handle to editEl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEl as text
%        str2double(get(hObject,'String')) returns contents of editEl as a double

el=str2double(get(hObject, 'String'));

if isnan(el)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    el=0;
end

if or(el<-180,el>=180)
    el=el-floor((el+180)/360)*360;
    set(hObject, 'String', num2str(el));
end

set(handles.sliderEl,'Value',el);

handles.plotdata.rot = 0;
set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);

handles.plotdata.el = el;

axes(handles.axesCAD);
view(handles.plotdata.az,handles.plotdata.el);

% Update handles structure

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editEl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxXYZdir.
function checkboxXYZdir_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxXYZdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxXYZdir

hoval=get(hObject,'Value');

le=0.12*min([handles.plotdata.dx,handles.plotdata.dy,handles.plotdata.dz]);
fosi=12;

tepo=1.3;

if hoval
    axes(handles.axesCAD);
    
    hold on
    
    handles.plotdata.hxyz(:)=plot3([handles.plotdata.xl(1) handles.plotdata.xl(1)+le],[handles.plotdata.yl(1) handles.plotdata.yl(1)],[handles.plotdata.zl(1) handles.plotdata.zl(1)],'r-',[handles.plotdata.xl(1) handles.plotdata.xl(1)],[handles.plotdata.yl(1) handles.plotdata.yl(1)+le],[handles.plotdata.zl(1) handles.plotdata.zl(1)],'g-',[handles.plotdata.xl(1) handles.plotdata.xl(1)],[handles.plotdata.yl(1) handles.plotdata.yl(1)],[handles.plotdata.zl(1) handles.plotdata.zl(1)+le],'b-');

    handles.plotdata.hTextX=text('String','x','Position',[handles.plotdata.xl(1)+tepo*le,handles.plotdata.yl(1),handles.plotdata.zl(1)],'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fosi,'Color','r');
    handles.plotdata.hTextY=text('String','y','Position',[handles.plotdata.xl(1),handles.plotdata.yl(1)+tepo*le,handles.plotdata.zl(1)],'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fosi,'Color','g');
    handles.plotdata.hTextZ=text('String','z','Position',[handles.plotdata.xl(1),handles.plotdata.yl(1),handles.plotdata.zl(1)+tepo*le],'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',fosi,'Color','b');
    
    xlim(handles.plotdata.xl);
    ylim(handles.plotdata.yl);
    zlim(handles.plotdata.zl);
    
    % Update handles structure
    guidata(handles.figure1, handles);
    
elseif handles.plotdata.hTextX>0
    
    delete(handles.plotdata.hxyz);
    delete(handles.plotdata.hTextX);
    delete(handles.plotdata.hTextY);
    delete(handles.plotdata.hTextZ);
    
    handles.plotdata.hxyz(:)=0;
    
    handles.plotdata.hTextX=0;
    handles.plotdata.hTextY=0;
    handles.plotdata.hTextZ=0;
    
end


% --- Executes on button press in checkboxAxis.
function checkboxAxis_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxAxis

hoval=get(hObject,'Value');

axes(handles.axesCAD);
if hoval
    axis on
else
    axis off
end

% --------------------------------------------------------------------
function menuFile_Callback(hObject, eventdata, handles)
% hObject    handle to menuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuHelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menuProductHelp_Callback(hObject, eventdata, handles)
% hObject    handle to menuProductHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('No help is available for the moment','Product Help','help');

% --------------------------------------------------------------------
function menuAbout_Callback(hObject, eventdata, handles)
% hObject    handle to menuAbout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msgbox('Written by Per Bergström 2012-03-01','About','help');

% --------------------------------------------------------------------
function menuOpenIGES_Callback(hObject, eventdata, handles)
% hObject    handle to menuOpenIGES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile({'*.igs';'*.iges'},'Select the IGES-file');

handles.plotdata.filename = file;

set(handles.textFile,'String',file);

handles.plotdata.saveas=false;

if ~isequal(file, 0)
    
    % Load parameter data from IGES-file
    [ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab([path,file]);
    
    % Plot the IGES object
    axes(handles.axesCAD);
    cla;
    hold on
    
    handles.plotdata.handlePlot=plotIGES(ParameterData,1,0);
    
    handles.plotdata.az = 0;
    handles.plotdata.el = 0;
    handles.plotdata.rot = 0;
    
    view(handles.plotdata.az,handles.plotdata.el);
    
    set(handles.editAz,'String', num2str(handles.plotdata.az));
    set(handles.sliderAz,'Value',handles.plotdata.az);
    
    set(handles.editEl,'String', num2str(handles.plotdata.el));
    set(handles.sliderEl,'Value',handles.plotdata.el);
    
    set(handles.editRot,'String', num2str(handles.plotdata.rot));
    set(handles.sliderRotation,'Value',handles.plotdata.rot);
    
    axis on
    set(handles.checkboxXYZdir,'Value',0);
    set(handles.checkboxAxis,'Value',1);
    set(handles.checkboxTransp,'Value',0);
    
    axis tight
    
    sc=0.2;
    
    xl=xlim;
    dx=xl(2)-xl(1);
    xl(1)=xl(1)-sc*dx;
    xl(2)=xl(2)+sc*dx;
    xlim(xl);
    
    yl=ylim;
    dy=yl(2)-yl(1);
    yl(1)=yl(1)-sc*dy;
    yl(2)=yl(2)+sc*dy;
    ylim(yl);
    
    zl=zlim;
    dz=zl(2)-zl(1);
    zl(1)=zl(1)-sc*dz;
    zl(2)=zl(2)+sc*dz;
    zlim(zl);
    
    handles.plotdata.igesfile = true;
    
    handles.plotdata.xl = xl;
    handles.plotdata.yl = yl;
    handles.plotdata.zl = zl;
    
    handles.plotdata.dx = dx;
    handles.plotdata.dy = dy;
    handles.plotdata.dz = dz;
    
    handles.plotdata.ParameterData = ParameterData;
    handles.plotdata.EntityType = EntityType;
    handles.plotdata.numEntityType = numEntityType;
    handles.plotdata.unknownEntityType = unknownEntityType;
    handles.plotdata.numunknownEntityType = numunknownEntityType;
    
    % Update handles structure
    guidata(handles.figure1, handles);
    
end


% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to menuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuTools_Callback(hObject, eventdata, handles)
% hObject    handle to menuTools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuRotate3D_Callback(hObject, eventdata, handles)
% hObject    handle to menuRotate3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rotate3d(handles.axesCAD)


% --------------------------------------------------------------------
function menuPan_Callback(hObject, eventdata, handles)
% hObject    handle to menuPan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pan(handles.axesCAD)


% --------------------------------------------------------------------
function menuOpenData_Callback(hObject, eventdata, handles)
% hObject    handle to ToolbarOpenData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.mat'},'Select file');

load([pathname,filename],'plotdata');

handles.plotdata=plotdata;

set(handles.textFile,'String',handles.plotdata.filename);

if handles.plotdata.igesfile
    
    % Plot the IGES object
    axes(handles.axesCAD);
    cla;
    hold on
    
    handles.plotdata.handlePlot=plotIGES(handles.plotdata.ParameterData,1,0);
    
elseif handles.plotdata.stlfile
    
    % Plot the STL object
    axes(handles.axesCAD);
    cla;
    hold on
    
    handles.plotdata.handlePlot=plotSTL(handles.plotdata.VertexData,[],0);
    
end

axis on
set(handles.checkboxXYZdir,'Value',0);
set(handles.checkboxAxis,'Value',1);
set(handles.checkboxTransp,'Value',0);

view(handles.plotdata.az,handles.plotdata.el);

xlim(handles.plotdata.xl);
ylim(handles.plotdata.yl);
zlim(handles.plotdata.zl);

handles.plotdata.rot = 0;
set(handles.editRot,'String', num2str(handles.plotdata.rot));
set(handles.sliderRotation,'Value',handles.plotdata.rot);

set(handles.editAz,'String', num2str(handles.plotdata.az));
set(handles.sliderAz,'Value',handles.plotdata.az);

set(handles.editEl,'String', num2str(handles.plotdata.el));
set(handles.sliderEl,'Value',handles.plotdata.el);


% Update handles structure

guidata(hObject,handles);

% --------------------------------------------------------------------
function menuSave_Callback(hObject, eventdata, handles)
% hObject    handle to menuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.plotdata.saveas
    plotdata=handles.plotdata;
    save(handles.plotdata.file,'plotdata');
else
    menuSaveDataAs_Callback(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function menuSaveDataAs_Callback(hObject, eventdata, handles)
% hObject    handle to menuSaveDataAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.plotdata.saveas=true;

[filename, pathname] = uiputfile({'*.mat';'*.*'},'Save as');

handles.plotdata.file=[pathname,filename];

plotdata=handles.plotdata;

save([pathname,filename],'plotdata');

% --- Executes on slider movement.
function sliderRotation_Callback(hObject, eventdata, handles)
% hObject    handle to sliderRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


rot=get(hObject,'Value');

set(handles.editRot,'String', num2str(rot));

axes(handles.axesCAD);

camroll(handles.axesCAD,rot-handles.plotdata.rot);

% Update handles structure

handles.plotdata.rot = rot;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sliderRotation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderRotation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editRot_Callback(hObject, eventdata, handles)
% hObject    handle to editRot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRot as text
%        str2double(get(hObject,'String')) returns contents of editRot as a double

rot=str2double(get(hObject, 'String'));

if isnan(rot)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
    rot=0;
end

if or(rot<-180,rot>=180)
    rot=rot-floor((rot+180)/360)*360;
    set(hObject, 'String', num2str(rot));
end

set(handles.sliderRotation,'Value',rot);

axes(handles.axesCAD);

camroll(handles.axesCAD,rot-handles.plotdata.rot);

% Update handles structure

handles.plotdata.rot = rot;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editRot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxTransp.
function checkboxTransp_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxTransp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxTransp

hoval=get(hObject,'Value');

if handles.plotdata.igesfile
    if hoval
        set(handles.plotdata.handlePlot{3},'FaceAlpha',0.4);
    else
        set(handles.plotdata.handlePlot{3},'FaceAlpha',1);
    end
elseif handles.plotdata.stlfile
    if hoval
        set(handles.plotdata.handlePlot,'FaceAlpha',0.4);
    else
        set(handles.plotdata.handlePlot,'FaceAlpha',1);
    end
end


% --------------------------------------------------------------------
function menuZoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to menuZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom(handles.axesCAD)

% --------------------------------------------------------------------
function menuZoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to menuZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

zoom(handles.axesCAD)


% --------------------------------------------------------------------
function menuExportFig_Callback(hObject, eventdata, handles)
% hObject    handle to menuExportFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uiputfile({'*.eps';'*.*'},'Export as');

saveas(handles.axesCAD,[pathname,filename]);
%saveas(handles.figure1,[pathname,filename]);


% --------------------------------------------------------------------
function menuNew_Callback(hObject, eventdata, handles)
% hObject    handle to menuNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

initialize_gui(hObject, handles);

% --------------------------------------------------------------------
function ToolbarOpenData_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolbarOpenData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menuOpenSTL_Callback(hObject, eventdata, handles)
% hObject    handle to menuOpenSTL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile({'*.stl'},'Select the STL-file');

handles.plotdata.filename = file;

set(handles.textFile,'String',file);

handles.plotdata.saveas=false;

if ~isequal(file, 0)
    
    % Load data from STL-file
    [VertexData,FVCD,isBinary]=stl2matlab([path,file]);
    
    % Plot the STL object
    axes(handles.axesCAD);
    cla;
    hold on
    
    handles.plotdata.handlePlot=plotSTL(VertexData,[],0);
    
    handles.plotdata.az = 0;
    handles.plotdata.el = 0;
    handles.plotdata.rot = 0;
    
    view(handles.plotdata.az,handles.plotdata.el);
    
    set(handles.editAz,'String', num2str(handles.plotdata.az));
    set(handles.sliderAz,'Value',handles.plotdata.az);
    
    set(handles.editEl,'String', num2str(handles.plotdata.el));
    set(handles.sliderEl,'Value',handles.plotdata.el);
    
    set(handles.editRot,'String', num2str(handles.plotdata.rot));
    set(handles.sliderRotation,'Value',handles.plotdata.rot);
    
    axis on
    set(handles.checkboxXYZdir,'Value',0);
    set(handles.checkboxAxis,'Value',1);
    set(handles.checkboxTransp,'Value',0);
    
    axis tight
    
    sc=0.2;
    
    xl=xlim;
    dx=xl(2)-xl(1);
    xl(1)=xl(1)-sc*dx;
    xl(2)=xl(2)+sc*dx;
    xlim(xl);
    
    yl=ylim;
    dy=yl(2)-yl(1);
    yl(1)=yl(1)-sc*dy;
    yl(2)=yl(2)+sc*dy;
    ylim(yl);
    
    zl=zlim;
    dz=zl(2)-zl(1);
    zl(1)=zl(1)-sc*dz;
    zl(2)=zl(2)+sc*dz;
    zlim(zl);
    
    handles.plotdata.stlfile = true;
    
    handles.plotdata.xl = xl;
    handles.plotdata.yl = yl;
    handles.plotdata.zl = zl;
    
    handles.plotdata.dx = dx;
    handles.plotdata.dy = dy;
    handles.plotdata.dz = dz;
    
    handles.plotdata.VertexData = VertexData;
    handles.plotdata.FVCD = FVCD;
    handles.plotdata.isBinary = isBinary;
    
    % Update handles structure
    guidata(handles.figure1, handles);
    
end
