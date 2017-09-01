function varargout = properties_polarization_map(varargin)
% PROPERTIES_POLARIZATION_MAP MATLAB code for properties_polarization_map.fig
%      PROPERTIES_POLARIZATION_MAP can not be called by itself.
%      POLARIZATION_MAP calls this GUI.
%
%      In PROPERTIES_POLARIZATION_MAP you can change the parameters that will be
%      used by POLARIZATION_MAP. When satisfied with the changes
%      press "Save". To cancel press the "Cancel"-button or just close the
%      window.
%
%      Abbreviations:
%      QWP:         Quarter-wave plate
%      HWP:         Half-wave plate
%      WP:          Waveplates
%      LP:          Linear polarizer
%      CP:          Circular polarization
%      LP:          Linear polarization (Poor choise of abbreviations,
%                   however it will be clear from surrounding text)
%      EM:          Power meter
%
%      This code has to be altered if a different brand of power meter or 
%      rotation motors are used and they have other properties that should 
%      be adjustable by the user.
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @properties_polarization_map_OpeningFcn, ...
                   'gui_OutputFcn',  @properties_polarization_map_OutputFcn, ...
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

% --- Executes just before properties_polarization_map is made visible.
function properties_polarization_map_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to properties_polarization_map (see VARARGIN)

% Choose default command line output for properties_polarization_map
handles.output = hObject;

% Takes in the current values of the properties:
prop = varargin{1};
% Enters the serial numbers:
for i = 1:3
    handles.table_SN.Data{i,2} = prop.HWSerialNum(i);
end
% Enters the analyzer step size:
handles.table_LP.Data{2} = prop.LP_stepsize;
% Enters the effect meter properties:
handles.table_EM.Data{3,2} = prop.wavelength;
handles.table_EM.Data{1,2} = sum(prop.meas_time.*[3600 60 1]);
handles.table_EM.Data{2,2} = prop.meas_freq;
% Enters the range of the WPs: (HWP: half-wave plate, QWP: quarter-wave
% plate)
handles.table_WP_range.Data(1,:) = prop.HWP_range;
handles.table_WP_range.Data(2,:) = prop.QWP_range;
% Enters the resolution of the WPs:
handles.table_WP_resolution.Data{1,2} = prop.HWP_resolution;
handles.table_WP_resolution.Data{2,2} = prop.QWP_resolution;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes properties_polarization_map wait for user response
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = properties_polarization_map_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.figure1);

% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Stores the serial numbers:
for i = 1:3
    out.HWSerialNum(i,1) = handles.table_SN.Data{i,2};
end
% Stores the analyzer step size:
out.LP_stepsize = handles.table_LP.Data{2};
% Stores the effect meter paramters:
out.wavelength = handles.table_EM.Data{3,2};
out.meas_time = floor(mod(handles.table_EM.Data{1,2},[0 3600 60])...
    ./[3600 60 1]);
out.meas_freq = handles.table_EM.Data{2,2};
% Stores the range of the WPs:
out.HWP_range = handles.table_WP_range.Data(1,:);
out.QWP_range = handles.table_WP_range.Data(2,:);
% Stores the resolution of the WPs:
out.HWP_resolution = handles.table_WP_resolution.Data{1,2};
out.QWP_resolution = handles.table_WP_resolution.Data{2,2};

% Makes the stored data the output:
handles.output = out;

% Update handles structure
guidata(hObject, handles);

% The user responce
uiresume(handles.figure1);

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% When using the exit button there should be send no data back
handles.output = [];

% Update handles structure
guidata(hObject, handles);

% The user responce
uiresume(handles.figure1);

% --- Executes when user attempts to close properties_polarization_map.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to properties_polarization_map window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% When using the exit button there should be send no data back:
handles.output = [];

% Update handles structure
guidata(hObject, handles);

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT
    uiresume(hObject);
else
    % The GUI is no longer waiting
    delete(hObject);
end
