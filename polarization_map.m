function varargout = polarization_map(varargin)
% POLARIZATION_MAP MATLAB code for polarization_map.fig
%      POLARIZATION_MAP, by itself, creates a new 
%      POLARIZATION_MAP or raises the existing singleton.
%
%      POLARIZATION_MAP takes no input parameters.
%
%      POLARIZATION_MAP opens an additional window containing the
%      activeX controls for the motors and power meter. Do not close this
%      window, it will automatically close when closing 
%      POLARIZATION_MAP.
%
%      "File, properties" in POLARIZATION_MAP opens 
%      PROPERTIES_POLARIZATION_MAP. Here various properties of the program 
%      can be changed.
%
%      Push the "connect"-button to connect to the power meter and motors.
%      Push the "disconnect"-button if any manual operation of the power 
%      meter is desired
%
%      After connecting start the polarization measurements by pushing the 
%      "start"-button in POLARIZATION_MAP.
%
%      The "stop"-button will not make the program stop instantly. The code
%      checks regularly if the button has been pushed. The code will not
%      stop untill it reaches such a check-point.
%
%      When the measurements are finished the "Calibration"-button is
%      enabled. Pushing this button fits a theoretical model to the data
%      and determines the QWP and HWP angles that provide linearly and
%      circularly polarized light. The results are displayed in the figure
%      and tables in the lower half of the GUI.
%
%      After the polarization measurements are finished "file, save" is 
%      enabled, which allows saving of the data and the calibration if the 
%      calibration was performed.
%
%      Previous measurements can be loaded by selecting "file, open".
%
%      To do a single polarization measurement enter the desired QWP and
%      HWP angle and push the "Single measurement"-button. Thereafter the
%      measurement can be saved by selecting "file, save single".
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
%      Symbols:
%      %!           Commands specific for our system (motors or power 
%                   meter) and HAS to be addapted to your system
%      %p           If a different brand of power meter is used these
%                   commands have to be addapted
%      %m           If a different brand of rotation motors are used these
%                   commands have to be addapted
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @polarization_map_OpeningFcn, ...
                   'gui_OutputFcn',  @polarization_map_OutputFcn, ...
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

% --- Executes just before polarization_map is made visible.
function polarization_map_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to polarization_map (see VARARGIN)

% Choose default command line output for polarization_map
handles.output = hObject;

% Opens a window containing all the activex controls necessary:
handles.actx_figure = figure;
% Top right is the power meter:
handles.EM = actxcontrol('LabMaxLowLevelControl.LabMaxLowLevCtl.1', [300 200 300 200]); %p
% Top left is the analyzer:
handles.LP = actxcontrol('MGMOTOR.MGMotorCtrl.1', [0 200 300 200]); %m
% Bottom left is the HWP:
handles.HWP = actxcontrol('MGMOTOR.MGMotorCtrl.1', [300 0 300 200]); %m
% Bottom right is the QWP:
handles.QWP = actxcontrol('MGMOTOR.MGMotorCtrl.1', [0 0 300 200]); %m

% Labels the axes of the figure that will contain the results:
axes(handles.axes_result)
set(handles.axes_result,'XAxisLocation','top','fontsize',8,'box','on');
xlabel('QWP [deg]')
ylabel('HWP [deg]')

% Sets the series number for the motors and energy meter:
handles.LP.HWSerialNum = 83845971; %!
handles.HWP.HWSerialNum = 83854989; %!
handles.QWP.HWSerialNum = 83846117; %!

% Sets the stepsize of the analyzer [deg]:
handles.LP_stepsize = 30;

% Sets the wavelength of the laser [nm]:
handles.wavelength = 890;
% Sets the acquisition time of each measurement with the power meter 
% [h m s]:
handles.meas_time = [0 0 1];
% Sets the frequency of the measurements with the power meter [1/s]:
handles.meas_freq = 10;

% The range and resolution of the HWP and QWP [deg] (this will determine 
% for which WP angles the polarization will be measured):
handles.HWP_range = [0, 60];
handles.QWP_range = [0, 150];
handles.HWP_resolution = 30;
handles.QWP_resolution = 30;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = polarization_map_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
% ------------------- FUNCTIONS IN THE GUI: -------------------------------
% --- Executes on button press in push_disconnect.

% --- Executes on button press in push_connect.
function push_connect_Callback(hObject, eventdata, handles)
% hObject    handle to push_connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Connects to the motors and power meter:
[handles, c] = connect_actx(handles);
% If something went wrong:
if ~c
    % Update handles structure
    guidata(hObject, handles);
    return;
end

% Enables the start, stop and single measurement button:
set(handles.push_start,'Enable','on');
set(handles.push_stop,'Enable','on');
set(handles.push_single,'Enable','on');

% Disables to the connect button and the properties menu:
set(hObject,'Enable','off');
set(handles.menu_prop,'Enable','off');
% Enables to the disconnect button:
set(handles.push_disconnect,'Enable','on')

% Update handles structure
guidata(hObject, handles);

function push_disconnect_Callback(hObject, eventdata, handles)
% hObject    handle to push_disconnect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disconnects from the motors and power meter:
handles = disconnect_actx(handles);

% Disables the start, stop and single measurement button:
set(handles.push_start,'Enable','off');
set(handles.push_stop,'Enable','off');
set(handles.push_single,'Enable','off');

% Enables to the connect button and properties menu:
set(handles.push_connect,'Enable','on');
set(handles.menu_prop,'Enable','on');
% Disables the disconnect button:
set(hObject,'Enable','off');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_start.
function push_start_Callback(hObject, eventdata, handles)
% hObject    handle to push_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Removes all previous results:
handles = clear_results(handles);

% Sets the stop variable to false:
setappdata(0 , 'stop', 0);

% Creates a HWP and QWP matrix containing all combinations of QWP and HWP 
% angles:
[handles.QWP_values, handles.HWP_values] = ndgrid(...
    handles.QWP_range(1):handles.QWP_resolution:handles.QWP_range(2),...
    handles.HWP_range(1):handles.HWP_resolution:handles.HWP_range(2));
[k,l] = size(handles.QWP_values);
% Creates a matrix that will contain all the analyzer and intensity
% measurements:
handles.LP_values = 0:handles.LP_stepsize:180;
n = length(handles.LP_values);
handles.intensity_measurements = zeros(k,l,n);
% The fitted parameters from the ellipticity analysis will be stored in
% these parameters:
handles.fitted_parameters.E_max = zeros(k,1);
handles.fitted_parameters.E_min = zeros(k,1);
handles.fitted_parameters.pol_angle = zeros(k,1);

% Creates a figure into which all the fittings will be plotted, so the user
% can see the progression of the measurements:
handles.overview_figure = figure('position',[10,50,(k*75)+75,(l*75)+75]);
for x = 1:k
    for y = 1:l
        % Rotates the QWP and HWP to specific angles:
        handles.QWP.MoveAbsoluteRot(0,handles.QWP_values(x,y),0,3,1); %m
        handles.HWP.MoveAbsoluteRot(0,handles.HWP_values(x,y),0,3,1); %m
        
        % Measures the polarization:
        [handles, I] = ellipticity_analyzer(handles);
        handles.intensity_measurements(x,y,:) = I;
        
        % Plots the results:
        [handles,~] = plot_polarization(handles, x, y);
        
        % Gives the GUI time to update:
        pause(0.1)
        % In order for the stop button to function we have to check 
        % throughout the code if it has been pushed:
        if getappdata(0,'stop')
            % Update handles structure
            guidata(hObject, handles);
            % Ends the calibration:
            return;
        end
    end
end

% Enables the save and calibration option:
set(handles.menu_save,'Enable','on');
set(handles.push_cali,'Enable','on');
% Dissables the "save single"-button:
set(handles.menu_save_single,'Enable','off');

% Enables the posibility to view different parameters of the measured
% polarization map:
set(handles.popup_result,'Enable','on');
handles = image_result(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_stop.
function push_stop_Callback(hObject, eventdata, handles)
% hObject    handle to push_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Sets the stop variable to true, this variable is accessible from all
% functions:
setappdata(0 , 'stop', 1);

% --- Executes on button press in push_single.
function push_single_Callback(hObject, eventdata, handles)
% hObject    handle to push_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Retrieves the HWP and QWP angle for the polarization measurement:
d = handles.table_single_settings.Data;
QWP = d(1); HWP = d(2);
% Rotates the QWP and HWP to specific angles:
handles.QWP.MoveAbsoluteRot(0,QWP,0,3,1); %m
handles.HWP.MoveAbsoluteRot(0,HWP,0,3,1); %m
% The analyzer angles:
handles.LP_values = 0:handles.LP_stepsize:180;
% Does the intensity measurement:
[handles, I] = ellipticity_analyzer(handles);

% Stores the single measurement:
handles.intensity_single = I;
handles.QWP_single = QWP;
handles.HWP_single = HWP;
handles.LP_single = handles.LP_values;

% Plots the result:
handles = plot_single(handles, I);

% Enables the "save single":
set(handles.menu_save_single,'Enable','on');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_cali.
function push_cali_Callback(hObject, eventdata, handles)
% hObject    handle to push_cali (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Does the phase and reflection fitting:
handles = find_phase_and_reflection(handles);
% Determines the QWP and HWP angles that give linear and circular
% polarization at the imaging plane:
handles = solve_for_LP_and_CP(handles);

% Plots the data:
handles = image_result(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popup_result.
function popup_result_Callback(hObject, eventdata, handles)
% hObject    handle to popup_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Images the polarization map as selected by the popup menu:
handles = image_result(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popup_result_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
% ------------------- CALIBRATION FUNCTIONS: ------------------------------

% --- Measures the ellipticity of the EM wave at the detector:
function [handles,I] = ellipticity_analyzer(handles)
% handles    structure with handles and user data (see GUIDATA)
% I          the measured intensity [W] for all analyzer angles stored as
%            an array

% The angles of the analyzer for which the intensity will be measured:
LP = handles.LP_values;
% Checks the start position:
start_pos = handles.LP.GetPosition_Position(0); %m
% For a linear polarizer 180 degrees = 0 degrees:
if start_pos >= 180
    start_pos = 360-start_pos;
end
% Changes the order if the start position is closer to 180deg than 0deg:
if start_pos > 90
    LP = LP(end:-1:1);
    change_order = 1;
else
    change_order = 0;
end

% Empty array to enter the measured intensity:
I = zeros(1,length(LP));

% Time required to take one measurement:
t = sum(handles.meas_time.*[3600 60 1]);

% Measures the intensity as a function of the angles of the analyzer:
for i = 1:length(LP)
    % Moves the analyzer to the correct angle:
    handles.LP.MoveAbsoluteRot(0,LP(i),0,3,1); %m
    % Asks the power meter to take a measurement:
    handles.EM.SendCommandOrQuery(0,'CONFigure:STATistics:STARt'); %p
    % Wait for the power meter to finish:
    pause(t);
    % Acquires the measurement (will be an array containing: AV, DEV, MIN, 
    % MAX, DOSE, MISSED, SEQ):
    handles.EM.SendCommandOrQuery(0,'STATistics:FETCh:NEXT?'); %p
    r = str2num(handles.EM.GetNextString(0)); %p
    % Check if the measurement has been received:
    while isempty(r)
        handles.EM.SendCommandOrQuery(0,'STATistics:FETCh:NEXT?'); %p
        r = str2num(handles.EM.GetNextString(0)); %p
    end
    % Sets the intensity equal to the average:
    I(i) = r(1);
end
% If the LP order was inverted:
if change_order == 1
    % Inverts the order of the measured intensities:
    I = I(end:-1:1);
end

% --- Plots the polarization of the laser for a given HWP and QWP angle
% --- combination as part of the polarization map:
function [handles,p] = plot_polarization(handles, i, j)
% handles    structure with handles and user data (see GUIDATA)
% i, j       the index of the QWP and HWP of the current measurement
% p          the parameters from the fitting

% Retrieves the correct measurement based on the current QWP and HWP angle:
I = squeeze(handles.intensity_measurements(i,j,:))';
LP = handles.LP_values;
% Adds the inverted datapoints such that the measurement is
% centro-symmetric:
I = [I I];
LP = [LP LP+180];
% Converts the intensity and angles from polar coordinates to cartesian
% coordinates:
[x,y] = pol2cart(deg2rad(LP),I);
% Does a curve fitting to determine the maximum and minimum electric
% amplitude and the rotation of the distribution:
p = lsqcurvefit(@intensity_fitting,[1,0.5,0],LP,I);
% Creates the fitted curve:
LP_fit = 0:360;
I_fit = (p(1)*cosd(LP_fit-p(3))).^2 + (p(2)*sind(LP_fit-p(3))).^2;
[x_fit,y_fit] = pol2cart(deg2rad(LP_fit),I_fit);

% Plots the fitted curve and the data points:
axes(handles.axes_single);
plot(x,y,'b.',x_fit,y_fit,'r')
axis('equal')

% Stores the fitted parameters:
E = sort(p(1:2));
handles.fitted_parameters.E_max(i,j) = E(2);
handles.fitted_parameters.E_min(i,j) = E(1);
% The fitted angle is always the angle between the x-axis and p(1), while 
% we want the angle between the x-axis and the maximum component of the 
% electric field. Therefor if p(1) is the minimum component we should add 
% 90 degrees to the angle:
if E(1) == p(1)
    p(3) = p(3)+90;
end
% Converts the angle to the range 0-180degrees:
p(3) = p(3) - floor(p(3)/180)*180;
handles.fitted_parameters.pol_angle(i,j) = p(3);

% Enters the measurement in the single measurement table:
handles.table_single_settings.Data = [handles.QWP_values(i,j) ...
    handles.HWP_values(i,j)];
handles.table_single_results.Data = [E(1)^2 E(2)^2 round(p(3))];

% Plots the polarization in the overview image:
figure(handles.overview_figure)
[k,l] = size(handles.QWP_values);
h = subplot('Position',[i/(k+1),(l-j)/(l+1),1/(k+1),1/(l+1)]);
plot(x,y,'b.',x_fit,y_fit,'r');
axis([-0.5 0.5 -0.5 0.5])
axis('equal')
set(h,'xtick',[],'ytick',[],'box','off','xcolor','w','ycolor','w')
if j == 1 && i == ceil(k/2)
    title({'QWP [deg]';num2str(handles.QWP_values(i,j))},'fontweight','normal')
elseif j == 1
    title(num2str(handles.QWP_values(i,j)),'fontweight','normal')
end
if i == 1 && j == ceil(l/2)
    ylabel({'HWP [deg]';num2str(handles.HWP_values(i,j))})
    set(h,'ycolor','k')
elseif i == 1
    ylabel(num2str(handles.HWP_values(i,j)))
    set(h,'ycolor','k')
end

% --- Plots the polarization of the laser for a single measurement:
function [handles] = plot_single(handles, I)
% handles    structure with handles and user data (see GUIDATA)
% I          the intensity measurements

% The analyzer angles:
LP = handles.LP_values;
% Adds the inverted datapoints such that the measurement is
% centro-symmetric:
I = [I I];
LP = [LP LP+180];
% Converts the intensity and angles from polar coordinates to cartesian
% coordinates:
[x,y] = pol2cart(deg2rad(LP),I);
% Does a curve fitting to determine the maximum and minimum electric
% amplitude and the tilt of the distribution:
p = lsqcurvefit(@intensity_fitting,[1,0.5,0],LP,I);
% Creates the fitted curve:
LP_fit = 0:360;
I_fit = (p(1)*cosd(LP_fit-p(3))).^2 + (p(2)*sind(LP_fit-p(3))).^2;
[x_fit,y_fit] = pol2cart(deg2rad(LP_fit),I_fit);

% Plots the fitted curve and the data points:
axes(handles.axes_single);
plot(x,y,'b.',x_fit,y_fit,'r')
axis('equal')

% The maximum and minimum electric field component:
E = sort(p(1:2));
% The fitted angle is always the angle between the x-axis and p(1), while 
% we want the angle between the x-axis and the maximum component of the 
% electric field. Therefor if p(1) is the minimum component we should add 
% 90 degrees to the angle:
if E(1) == p(1)
    p(3) = p(3)+90;
end
% Converts the angle to the range 0-180degrees:
p(3) = p(3) - floor(p(3)/180)*180;

% Enters the results in the table:
handles.table_single_results.Data = [E(1)^2 E(2)^2 round(p(3))];

% --- Determines the phase and reflection introduced by optical elements in
% --- the microscope:
function [handles] = find_phase_and_reflection(handles)
% handles    structure with handles and user data (see GUIDATA)

% Gathers the input variables for the fitting:
[k,l] = size(handles.QWP_values);
LP = repmat(reshape(handles.LP_values,1,1,[]),k,l);
QWP = repmat(handles.QWP_values,1,1,length(handles.LP_values));
HWP = repmat(handles.HWP_values,1,1,length(handles.LP_values));
X = [HWP(:) QWP(:) LP(:)];
I = handles.intensity_measurements;
% Fits the polarization measurements:
P = lsqcurvefit(@polarization_fitting,[0.1,1,0,0,0,0],...
    X,I(:));

% The found parameters:
% Max intensity:
handles.fitted_parameters.I0 = P(1);
% The reflection coefficient:
handles.fitted_parameters.gamma = P(2);
% The angular difference between the reflection axis and the HWP (at 0deg):
handles.fitted_parameters.HWP_diff = P(3);
% The angular difference between the reflection axis and the QWP (at 0deg):
handles.fitted_parameters.QWP_diff = P(4);
% The angular difference between the reflection axis and the LP (at 0deg):
handles.fitted_parameters.LP_diff = P(5);
% The induced phase along the reflection axis:
handles.fitted_parameters.delta = P(6);

% --- Finds the HWP and QWP combinations that create LP and CP light:
function [handles] = solve_for_LP_and_CP(handles)
% handles    structure with handles and user data (see GUIDATA)

% The parameters from the phase and reflection fitting:
delta = handles.fitted_parameters.delta;
gamma = handles.fitted_parameters.gamma;
QWP_diff = handles.fitted_parameters.QWP_diff;
HWP_diff = handles.fitted_parameters.HWP_diff;

% Solves theoretical model for circularly polarized light:
% Equation 1:
syms phi;
x = atan(-gamma*sind(delta)/(gamma*cosd(delta)*tan(phi)-1));
y = atan((gamma*cosd(delta)-tan(phi))/(gamma*sind(delta)*tan(phi)));
phi = solve(x == y,phi,'IgnoreAnalyticConstraints',true); k = 1; 
phi = eval(phi);
WP_diff(1:2) = rad2deg(eval(x));
QWP_rel(1:2) = rad2deg(phi);
% Equation 2:
syms phi;
x =  atan(-gamma*sind(delta)/(gamma*cosd(delta)*tan(phi)+1));
y = atan((gamma*cosd(delta)+tan(phi))/(gamma*sind(delta)*tan(phi)));
phi = solve(x == y,phi,'IgnoreAnalyticConstraints',true); k = 1; 
phi = eval(phi);
WP_diff(3:4) = rad2deg(eval(x));
QWP_rel(3:4) = rad2deg(phi);
HWP_rel = (WP_diff + QWP_rel)/2;
% The results:
HWP_circ = HWP_rel + HWP_diff;
QWP_circ = QWP_rel + QWP_diff;
HWP_circ = (HWP_circ-90*floor(HWP_circ/90));
QWP_circ = (QWP_circ-180*floor(QWP_circ/180));
% Enters the results in appropriate table:
handles.table_CP.Data = [QWP_circ;HWP_circ];

% Solves theoretical model for LP light:
WP_diff = (0:0.5:180);
QWP_rel = zeros(2,length(WP_diff));
for i = 1:length(WP_diff)
    QWP_rel(1,i) = 0.5*atand(-tand(delta)*sind(2*WP_diff(i)));
    QWP_rel(2,i) = 0.5*(180+atand(-tand(delta)*sind(2*WP_diff(i))));
end
HWP_rel = 0.5*([WP_diff;WP_diff]+QWP_rel);
% The results:
HWP_LP = HWP_rel(:) + HWP_diff;
QWP_LP = QWP_rel(:) + QWP_diff;
HWP_LP = (HWP_LP-90*floor(HWP_LP/90));
QWP_LP = (QWP_LP-180*floor(QWP_LP/180));

% The polarization angle and intensity is also of interest for the settings 
% that create linearly polarized light:
I_LP = zeros(length(HWP_LP),1);
angle_LP = zeros(length(HWP_LP),1);
% For all HWP and QWP combinations that provide linearly polarized light:
for i = 1:length(HWP_LP)
    % The variables:
    LP = 0:10:360; n = length(LP);
    QWP = repmat(QWP_LP(i),1,n);
    HWP = repmat(HWP_LP(i),1,n);
    X = [HWP(:) QWP(:) LP(:)];
    % The parameters:
    P(1) = handles.fitted_parameters.I0;
    P(2) = handles.fitted_parameters.gamma;
    P(3) = handles.fitted_parameters.HWP_diff;
    P(4) = handles.fitted_parameters.QWP_diff;
    P(5) = handles.fitted_parameters.LP_diff;
    P(6) = handles.fitted_parameters.delta;
    % Finds the intensity distribution as a function of the analyzer (LP):
    I = polarization_fitting(P,X);
    % Does a fitting of the intensity to get more precise results:
    p = lsqcurvefit(@intensity_fitting,[1,0.5,0],LP,I');
    % Finds the maximum intensity:
    E = sort(p(1:2));
    I_LP(i) = E(2)^2;
    % The fitted angle is always the angle between the x-axis and p(1),  
    % while we want the angle between the x-axis and the maximum component 
    % of the electric field. Therefor if p(1) is the minimum component we 
	% should add 90 degrees to the angle:
    if E(1) == p(1)
        p(3) = p(3)+90;
    end
    p(3) = p(3) - floor(p(3)/180)*180;
    angle_LP(i) = round(p(3));
end
% Sorts the values based on the polarization angle:
[angle_LP,i] = sort(angle_LP);
% Stores the linearly polarized data in the appropriate table:
handles.table_LP.Data = [QWP_LP(i) HWP_LP(i) I_LP(i) angle_LP];

% --- Deletes all previous results:
function [handles] = clear_results(handles)
% handles    structure with handles and user data (see GUIDATA)

% Sets all table values to zero:
handles.table_single_results.Data = [];
handles.table_CP.Data = [];
handles.table_LP.Data = [];
handles.fitted_parameters = [];

% Removes the previous fitted ellipses and other results:
cla(handles.axes_single)
cla(handles.axes_result)

% Disables the appropriate buttons and options:
set(handles.menu_save,'Enable','off');
set(handles.push_cali,'Enable','off');
set(handles.popup_result,'Enable','off');

% --- Images a calibrated parameter in the result window:
function [handles] = image_result(handles)
% handles    structure with handles and user data (see GUIDATA)

% Collects which result should be imaged:
v = get(handles.popup_result,'value');

switch v
    case 1 % Intensity:
        E_min = handles.fitted_parameters.E_min;
        E_max = handles.fitted_parameters.E_max;
        r = E_min.^2 + E_max.^2;
        c = [0 max(r(:))];
        cm = 'jet';
    case 2 % Minimum intensity:
        E_min = handles.fitted_parameters.E_min;
        r = E_min.^2;
        c = [0 max(r(:))];
        cm = 'jet';
    case 3 % Maximum intensity:
        E_max = handles.fitted_parameters.E_max;
        r = E_max.^2;
        c = [0 max(r(:))];
        cm = 'jet';
    case 4 % Ratio between min and max I:
        E_min = handles.fitted_parameters.E_min;
        E_max = handles.fitted_parameters.E_max;
        r = E_min./E_max;
        c = [0 1];
        cm = 'jet';
    case 5 % Maximum intensity angle:
        r = handles.fitted_parameters.pol_angle;
        c = [0 180];
        cm = 'hsv';
end
[k,l] = size(r);
% Plots the result:
axes(handles.axes_result)
x = handles.QWP_values(:,1)';
y = handles.HWP_values(1,:);
imagesc(x,y,r',c)
axis equal; axis tight;
xlabel('QWP [deg]'); ylabel('HWP [deg]');
set(gca,'XAxisLocation','top','XTick',x(1:2:end),'YTick',y(1:2:end))
colormap(cm); colorbar;

% Checks if the calibration has been done:
if isempty(handles.table_CP.Data) || isempty(handles.table_LP.Data)
    % If no, the function is terminated here
    return;
end

% Plots the data for circularly and linearly polarized light:
D_CP = handles.table_CP.Data;
D_LP = handles.table_LP.Data;
hold on;
plot(D_CP(1,:),D_CP(2,:),'ko', D_LP(:,1), D_LP(:,2), 'kx','linewidth',2);
axis([-5 185 -5 95])
set(gca,'XTick',0:20:180,'YTick',0:20:90)
legend('Circular polarization','Linear polarization');
hold off;

% -------------------------------------------------------------------------
% ------------------- Functions that are used for fitting -----------------

% --- The intensity fitting to determine the polarization:
function I = intensity_fitting(p,a)
% p         The parameters that will be fitted. Here p(1) and p(2) are two
%           perpendicular electric amplitudes. Where the one is the maximum 
%           and the other the minimum amplitude. p(3) is the angle between 
%           the x-axis and the direction of p(1) (independent of if this is 
%           the minimum or maximum).
% a         An array containing an overview of the analyzer angles for 
%           which the intensity was measured 
% I         An array containing the intensities that correspond to the 
%           given input angles and fitted parameters.

I = (p(1)*cosd(a-p(3))).^2 + (p(2)*sind(a-p(3))).^2;

% --- The relation between intensity, analyzer, HWP and QWP angles and
% --- microscope specific parameters:
function I = polarization_fitting(P,X)
% P         Microscope parameters
% X         The variables (i.e. HWP, QWP and analyzer angles)
% I         An array containing the intensities that correspond to the 
%           given input angles and parameters.

% The variables:
HWP = X(:,1); QWP = X(:,2); LP = X(:,3);

% The paramters:
I0 = P(1); gamma = P(2); HWP_diff = P(3); 
QWP_diff = P(4); LP_diff = P(5); delta = P(6);

% The relative angles:
theta = HWP-HWP_diff; phi = QWP-QWP_diff; alpha = LP-LP_diff;

% The function describing how the intensity depends on the analyzer, QWP
% and HWP angle:
d1 = -gamma*(cosd(delta)*sind(phi).*sind(2*theta-phi) + ...
    sind(delta)*cosd(phi).*cosd(2*theta-phi));
d2 = -gamma*(sind(delta)*sind(phi).*sind(2*theta-phi) - ...
    cosd(delta)*cosd(phi).*cosd(2*theta-phi));
d3 = sind(phi).*cosd(2*theta-phi);
d4 = cosd(phi).*sind(2*theta-phi);

I = I0*((d1.^2 + d2.^2).*(cosd(alpha).^2) + ...
    (d3.^2 + d4.^2).*(sind(alpha).^2) + ...
    2*(d1.*d3 + d2.*d4).*sind(alpha).*cosd(alpha));

% -------------------------------------------------------------------------
% ----------- Connecting and disconnecting the activeX controls: ----------

% --- Connects to the motors and power meter:
function [handles, connected] = connect_actx(handles)
% handles    structure with handles and user data (see GUIDATA)
% connected  boolean stating if the connection was successfull

% Connects to the power meter:
handles.EM.Initialize(); %p

% Sets the communication mode to USB:
set(handles.EM,'CommunicationMode',1); %p
% Checks the serial number:
if isempty(handles.EM.SerialNumber(0)) %p
    handles = disconnect_actx(handles);
    warning('Could not connect to the power meter.');
    connected = 0;
    return
else
    connected = 1;
end

% Sets the units to watt:
handles.EM.SendCommandOrQuery(0,'CONFigure:MEASure W'); %p
% Sets the laser wavelength:
handles.EM.SendCommandOrQuery(0,'CONFigure:WAVElength:CORRection ON'); %p
handles.EM.SendCommandOrQuery(0,['CONFigure:WAVElength:WAVElength ' ...
    num2str(handles.wavelength)]); %p
% Sets the acquisition time of the power meter:
handles.EM.SendCommandOrQuery(0,['CONFigure:STATistics:BSIZe:TIME ' ...
    num2str(handles.meas_time(1),'%02.0f') ...
    ':' num2str(handles.meas_time(2),'%02.0f') ...
    ':' num2str(handles.meas_time(3),'%02.0f')]); %p
% Sets the frequency of the power meter:
handles.EM.SendCommandOrQuery(0,['CONFigure:STATistics:RATE:FREQuency ' ...
    num2str(handles.meas_freq)]); %p
% Determines if the data sent from the machine should have any headers:
handles.EM.SendCommandOrQuery(0,'CONFigure:READings:HEADers OFF'); %p
% Sets that we want the last available record and not a stream:
handles.EM.SendCommandOrQuery(0,'CONFigure:READings:CONTinuous LAST'); %p

% Changes the display to the statistics mode:
handles.EM.SendCommandOrQuery(0,'CONFigure:DISPlay:STATistics'); %p
% Enables the collection of statistical data for sending to the host:
handles.EM.SendCommandOrQuery(0,'STATistics:INITiate'); %p
% Collects all data in case there is any registered:
handles.EM.SendCommandOrQuery(0,'STATistics:FETCh:ALL?'); %p
handles.EM.GetNextString(0); %p

% Initializes the motors:
handles.LP.StartCtrl; %m
handles.HWP.StartCtrl; %m
handles.QWP.StartCtrl; %m
% Gives the GUI time to update:
pause(0.5)

% Homes the motors:
c(1) = handles.LP.MoveHome(0,1); %m
c(2) = handles.HWP.MoveHome(0,1); %m
c(3) = handles.QWP.MoveHome(0,1); %m

if sum(c) ~= 0
    handles = disconnect_actx(handles);
    warning('Could not connect to the motors.');
    connected = 0;
    return
end

% --- Disconnects the motors and power meters:
function [handles] = disconnect_actx(handles)
% handles    structure with handles and user data (see GUIDATA)

% Stops the motors:
handles.LP.StopCtrl; %m
handles.HWP.StopCtrl; %m
handles.QWP.StopCtrl; %m

% Disables the collection of data:
handles.EM.SendCommandOrQuery(0,'STATistics:ABORt'); %p
% Discontect the meter:
handles.EM.DeInitialize(); %p

% -------------------------------------------------------------------------
% -------------------------- MENU FUNCTIONS: ------------------------------

% --------------------------------------------------------------------
function menu_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clears all previous results:
handles = clear_results(handles);

% Lets the user select the matlab file that will be opened:
[file, path] = uigetfile('*.mat','Select the MATLAB file containing the results');
% Checks if the action is canceled:
if file == 0
    return;
end

% Loads the variables to the workspace:
load([path file]);
try
    handles.QWP_values = measured_parameters.QWP;
    handles.HWP_values = measured_parameters.HWP;
    handles.intensity_measurements = measured_parameters.I;
    handles.LP_values = measured_parameters.LP;
    handles.fitted_parameters = fitted_parameters.polarization;
    if isfield(fitted_parameters,'map')
        if isfield(fitted_parameters.map,'I0')
            handles.fitted_parameters.I0 = fitted_parameters.map.I0;
            handles.fitted_parameters.gamma = fitted_parameters.map.gamma;
            handles.fitted_parameters.HWP_diff = fitted_parameters.map.HWP_diff;
            handles.fitted_parameters.QWP_diff = fitted_parameters.map.QWP_diff;
            handles.fitted_parameters.LP_diff = fitted_parameters.map.LP_diff;
            handles.fitted_parameters.delta = fitted_parameters.map.delta;
            handles = solve_for_LP_and_CP(handles);
        end
    end
catch
    disp('Could not load the data');
    return;
end

% Enables the save option:
set(handles.menu_save,'Enable','on');
set(handles.push_cali,'Enable','on');
% Enables the posibility to image results:
set(handles.popup_result,'Enable','on');
handles = image_result(handles);

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lets the user choose file name and folder:
[file, path] = uiputfile('*.mat','Save Results As');
% Checks if the action is canceled:
if file == 0
    return;
end

% Collects the variables that will be saved:
measured_parameters.QWP = handles.QWP_values;
measured_parameters.HWP = handles.HWP_values;
measured_parameters.I = handles.intensity_measurements;
measured_parameters.LP = handles.LP_values;
fitted_parameters.polarization.E_max = handles.fitted_parameters.E_max;
fitted_parameters.polarization.E_min = handles.fitted_parameters.E_min;
fitted_parameters.polarization.pol_angle = handles.fitted_parameters.pol_angle;
if isfield(handles.fitted_parameters,'I0')
    if ~isempty(handles.fitted_parameters.I0)
        fitted_parameters.map.I0 = handles.fitted_parameters.I0;
        fitted_parameters.map.gamma = handles.fitted_parameters.gamma;
        fitted_parameters.map.HWP_diff = handles.fitted_parameters.HWP_diff;
        fitted_parameters.map.QWP_diff = handles.fitted_parameters.QWP_diff;
        fitted_parameters.map.LP_diff = handles.fitted_parameters.LP_diff;
        fitted_parameters.map.delta = handles.fitted_parameters.delta;
    end
end
% Saves the results:
save([path file],'measured_parameters','fitted_parameters');

if ~isempty(handles.table_CP.Data) && ~isempty(handles.table_LP.Data)
    % Opens a file to write the HWP and QWP values to:
    FID = fopen([path file(1:end-4) '.txt'],'w');
    % Stores the circular polarization settings:
    t = 'Circular polarization:\nQWP [deg] HWP [degrees]\n';
    fprintf(FID, t);
    d = handles.table_CP.Data;
    f = '%d %d\n';
    fprintf(FID,f,d);
    % Stores the linear polarization settings:
    t = 'Linear polarization:\nQWP [deg] HWP [degrees] Intensity [W] Polarization angle [deg]\n';
    fprintf(FID, t);
    d = handles.table_LP.Data;
    f = '%d %d %f %d\n';
    fprintf(FID,f,d');
    fclose(FID);
end

% --------------------------------------------------------------------
function menu_save_single_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lets the user choose file name and folder:
[file, path] = uiputfile('*.mat','Save Single Measurement As');

% Collects the variables that will be saved:
I = handles.intensity_single;
LP = handles.LP_single;
QWP = handles.QWP_single;
HWP = handles.HWP_single;
% Saves the variables:
save([path file],'I','LP','QWP','HWP');
  
% --------------------------------------------------------------------
function menu_prop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Collects all the parameters that can be changed:
in.HWSerialNum(1,1) = handles.LP.HWSerialNum;
in.HWSerialNum(2,1) = handles.HWP.HWSerialNum;
in.HWSerialNum(3,1) = handles.QWP.HWSerialNum;
in.LP_stepsize = handles.LP_stepsize;
in.wavelength = handles.wavelength;
in.meas_time = handles.meas_time;
in.meas_freq = handles.meas_freq;
in.HWP_resolution = handles.HWP_resolution;
in.QWP_resolution = handles.QWP_resolution;
in.HWP_range = handles.HWP_range;
in.QWP_range = handles.QWP_range;

% Calls the properties editor:
out = properties_polarization_map(in);
% Checks if the action was canceled:
if isempty(out)
    return
end

% Stores all the data that could be changed:
handles.LP.HWSerialNum = out.HWSerialNum(1);
handles.HWP.HWSerialNum = out.HWSerialNum(2);
handles.QWP.HWSerialNum = out.HWSerialNum(3);
handles.LP_stepsize = out.LP_stepsize;
handles.wavelength = out.wavelength;
handles.meas_time = out.meas_time;
handles.meas_freq = out.meas_freq;
handles.HWP_resolution = out.HWP_resolution;
handles.QWP_resolution = out.QWP_resolution;
handles.HWP_range = out.HWP_range;
handles.QWP_range = out.QWP_range;

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% -------------------------------------------------------------------------
% ------------------- WHEN THE GUI IS CLOSED: -----------------------------

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disconnects all activex controls if they are not already deactivated:
handles = disconnect_actx(handles);
% Closes the figure containing the activex controls:
close(handles.actx_figure);

% Hint: delete(hObject) closes the figure
delete(hObject);
