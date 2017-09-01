function varargout = polarization_imaging(varargin)
% POLARIZATION_IMAGING MATLAB code for polarization_imaging.fig
%      POLARIZATION_IMAGING, by itself, creates a new POLARIZATION_IMAGING 
%      or raises the existing singleton.
%
%      POLARIZATION_IMAGING takes no input arguments.
%
%      POLARIZATION_IMAGING opens an additional window containing the
%      activeX controls for the motors. Do not close this window, it will 
%      automatically close when closing POLARIZATION_IMAGING
%
%      Start by opening a previously measured polarization map 
%      (calibration), "file,open". 
%
%      The "Calculate"-button will calculate the a given number of 
%      quarter-wave plate (QWP) and half-wave plate (HWP) combinations that 
%      provide linearly polarized (LP) light. The number of
%      combinations can be changed by altering the number in the edit box.
%      The combinations that provide circularly polarized (CP) light are 
%      calculated as well.
%
%      It is possible to change directory for where the polarization 
%      measurements will be saved.
%
%      Select the desired polarization measurements from the list of QWP 
%      and HWP combinations.
%
%      When ready to take the polarization measurements press the 
%      "start"-button. 
%
%      The "stop"-button will not make the program stop instantly. The code
%      checks regularly if the button has been pushed. The code will not
%      stop untill it reaches such a check-point.
%
%      Abbreviations:
%      QWP:         Quarter-wave plate
%      HWP:         Half-wave plate
%      WP:          Waveplates
%      LP:          Linear polarizer
%      CP:          Circular polarization
%      LP:          Linear polarization (Poor choise of abbreviations,
%                   however it will be clear from surrounding text)
%
%      Symbols:
%      %!           Commands specific for our system (motors or microscope) 
%                   and HAS to be addapted to your system
%      %M           Commands specific for a Leica TCS SP8 confocal 
%                   microscope
%      %m           If a different brand of rotation motors are used these
%                   commands have to be addapted
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @polarization_imaging_OpeningFcn, ...
                   'gui_OutputFcn',  @polarization_imaging_OutputFcn, ...
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

% --- Executes just before polarization_imaging is made visible.
function polarization_imaging_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to polarization_imaging (see VARARGIN)

% Choose default command line output for polarization_imaging
handles.output = hObject;

% Opens a window containing all the activex controls necessary:
handles.fig = figure;
p = get(handles.fig,'Position');
p(3) = p(3)/2;
set(handles.fig,'Position',p);
% Top is the HWP:
handles.HWP = actxcontrol('MGMOTOR.MGMotorCtrl.1', [0 200 300 200]);%m
% Bottom is the QWP:
handles.QWP = actxcontrol('MGMOTOR.MGMotorCtrl.1', [0 0 300 200]);%m

% Sets the series number for the motors:
handles.HWP.HWSerialNum = 83854989;%!
handles.QWP.HWSerialNum = 83846117;%!

% Creates a tcpip object for the microscope:
handles.microscope = tcpip('localhost',8895,'NetworkRole','client');%M

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = polarization_imaging_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% -------------------------------------------------------------------------
% ------------------- FUNCTIONS IN THE GUI: -------------------------------

% --- Executes on button press in push_calc.
function push_calc_Callback(hObject, eventdata, handles)
% hObject    handle to push_calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Calculates the appropriate QWP and HWP values for CP and LP light:
handles = solve_for_LP_and_CP(handles);

% Images the results in the axis:
handles = image_ratio(handles);

% Stores the data in the table:
n = length(handles.CP.QWP);
for i = 1:n
    D{i,1} = false;
    D{i,2} = 'CP';
    D{i,3} = handles.CP.HWP(i);
    D{i,4} = handles.CP.QWP(i);
    D{i,5} = '-';
end
m = length(handles.LP.QWP);
for i = 1:m
    D{i+n,1} = false;
    D{i+n,2} = handles.LP.pol_angle(i);
    D{i+n,3} = handles.LP.HWP(i);
    D{i+n,4} = handles.LP.QWP(i);
    D{i+n,5} = handles.LP.I(i);
end
handles.table_data.Data = D;

% Enables the start button:
set(handles.push_start,'enable','on');
set(handles.push_stop,'enable','on');

% Update handles structure
guidata(hObject, handles);

function edit_num_angels_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num_angels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_num_angels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num_angels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in push_dir.
function push_dir_Callback(hObject, eventdata, handles)
% hObject    handle to push_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Asks for the file and path name:
[file_name, path_name] = uiputfile({'*.*','Base name'}, ...
    'Select the location to save the images and information file');
% If nothing was chosen:
if ~file_name
    return
end
% Enters the location and base name in the text box:
set(handles.text_dir,'string',[path_name file_name]);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_start.
function push_start_Callback(hObject, eventdata, handles)
% hObject    handle to push_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Sets the stop variable to false:
setappdata(0 , 'stop', 0);

% Connects to the motors and microscope:
[handles, c] = connect_all(handles);
if ~c
    % Update handles structure
    guidata(hObject, handles);
    return;
end
% Update handles structure
guidata(hObject, handles);

% Gives the GUI time to update:
pause(0.1)
% In order for the stop button to function we have to check throughout the
% code if it has been pushed:
if getappdata(0,'stop')
    % Disconnects the motors and microscope:
    handles = disconnect_all(handles);
    % Update handles structure
    guidata(hObject, handles);
    % Ends the calibration:
    return;
end

% Opens a file to enter the information about the images:
FID = fopen([get(handles.text_dir,'string') '.txt'],'w');
% Makes the titles:
t = 'Image number [#] Polarization angle [deg] HWP [deg] QWP [degrees] Intensity [W] \n';
fprintf(FID, t);

% The total number of polarization measurements:
[k,~] = size(handles.table_data.Data);
im = 1;
for i = 1:k
    if handles.table_data.Data{i,1}
        % The WP angles:
        HWP = handles.table_data.Data{i,3};
        QWP = handles.table_data.Data{i,4};
        % Moves the WP motors to the correct angles:
        handles.HWP.MoveAbsoluteRot(0,HWP,0,3,1); %m
        handles.QWP.MoveAbsoluteRot(0,QWP,0,3,1); %m


        % Makes the code to do a scan:
        code = [handles.code_pre 'startscan']; %M
        code = double(code'); %M
        % Sends the code:
        fwrite(handles.microscope,code); %M
        % Gets the answer:
        m = ' ';
        t = 1;
        while isempty(strfind(m,'/inf:scanfinished')) && t<100 %M
            while handles.microscope.BytesAvailable == 0 %M
                pause(0.05)
            end
            pause(3)
            m = fread(handles.microscope,handles.microscope.BytesAvailable); %M
            m = char(m');
            t = t+1;
        end
        % If the image could not be made:
        if t==100
            warning('Could not make an image.');
            % Disconnects from motors and microscope:
            handles = disconnect_all(handles);
            % Update handles structure
            guidata(hObject, handles);
            return
        end
        % The place the images are stored:
        path_pre = 'D:\MatrixScreenerImages\'; %M
        if isempty(strfind(m,'.tif'))
            warning('Could not find the image.');
            % Disconnects the motors and microscope:
            handles = disconnect_all(handles);
            % Update handles structure
            guidata(hObject, handles);
            return
        end
        si = strfind(m,'relpath:') + 8; %M
        sj = strfind(m,'.tif') + 4; %M
        path_file = m(si:sj);
        source = [path_pre path_file]; %M
        % Copies the image to the desired folder:
        destination = [get(handles.text_dir,'string') '_' num2str(im) '.tif'];
        copyfile(source,destination)

        % Stores the information about the image to a text file:
        if ischar(handles.table_data.Data{i,2})
            d = [num2str(im) ' ' ...
                handles.table_data.Data{i,2} ' ' ...
                num2str(handles.table_data.Data{i,3}) ' '...
                num2str(handles.table_data.Data{i,4}) ' '...
                handles.table_data.Data{i,5} '\n'];
            c = c+1;
        else
            d = [num2str(im) ' ' ...
                num2str(handles.table_data.Data{i,2}) ' '...
                num2str(handles.table_data.Data{i,3}) ' '...
                num2str(handles.table_data.Data{i,4}) ' '...
                num2str(handles.table_data.Data{i,5}) '\n'];
        end
        fprintf(FID,d);

        % Gives the GUI time to update:
        pause(0.1)
        % Checks if the stop button has been pushed:
        if getappdata(0,'stop')
            % Stops the for-loop:
            break
        end
        im = im+1;
        pause(0.05)
    end
end
% Closes the text file:
fclose(FID);

% Disconnects the motors and microscope:
handles = disconnect_all(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in push_stop.
function push_stop_Callback(hObject, eventdata, handles)
% hObject    handle to push_stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Sets the stop variable to true, this variable is accessable from all
% functions:
setappdata(0 , 'stop', 1);

% -------------------------------------------------------------------------
% -------------------------- MENU FUNCTIONS: ------------------------------

% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Asks for the file and path name:
[file, path] = uigetfile({'*.mat','MAT files (*.mat)'}, 'Select the data file');
% If nothing was chosen:
if ~file
    return
end

% Clears all previous results:
[handles] = clear_results(handles);

% Opens the data:
load([path file]);
try
    handles.QWP_values = measured_parameters.QWP;
    handles.HWP_values = measured_parameters.HWP;
    handles.intensity_measurements = measured_parameters.I;
    handles.LP_values = measured_parameters.LP;
    handles.fitted_parameters.E_max = fitted_parameters.polarization.E_max;
    handles.fitted_parameters.E_min = fitted_parameters.polarization.E_min;
    handles.fitted_parameters.pol_angle = fitted_parameters.polarization.pol_angle;
catch
    disp('Could not load the data');
    return;
end

% Does the fitting:
handles = find_phase_and_reflection(handles);

% Enables the calculate button:
set(handles.push_calc,'Enable','on')

% Update handles structure
guidata(hObject, handles);

% -------------------------------------------------------------------------
% -------------------------- Sub-fuctions: --------------------------------

% --- Finds the HWP and QWP combinations that create LP and CP light:
function [handles] = solve_for_LP_and_CP(handles)
% handles    structure with handles and user data (see GUIDATA)

% The parameters from the phase and reflection fitting:
delta = handles.fitted_parameters.delta;
gamma = handles.fitted_parameters.gamma;
QWP_diff = handles.fitted_parameters.QWP_diff;
HWP_diff = handles.fitted_parameters.HWP_diff;

% Solves for circular polarized light:
syms phi;
% Equation 1:
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
% Stores the QWP and HWP values for CP light:
handles.CP.QWP = QWP_circ;
handles.CP.HWP = HWP_circ;

% Solves for LP light:
n = str2num(get(handles.edit_num_angels,'string'));
WP_diff = (0:180/n:180);
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

% For all HWP and QWP combinations that provide linear polarized light:
I_LP = zeros(length(HWP_LP),1);
angle_LP = zeros(length(HWP_LP),1);
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
    % Does the fitting of the intensity to get more precise results:
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
% Sorts the values based on the polarization angle and removes duplicates:
D = [angle_LP,QWP_LP,HWP_LP,I_LP];
D = unique(D, 'rows');
% Stores the linear polarized data:
handles.LP.pol_angle = D(:,1);
handles.LP.QWP = D(:,2);
handles.LP.HWP = D(:,3);
handles.LP.I = D(:,4);

% --- Deletes all previous results:
function [handles] = clear_results(handles)
% handles    structure with handles and user data (see GUIDATA)

% Sets all table values to zero:
handles.table_data.Data = [];
handles.CP = [];
handles.LP = [];
handles.fitted_parameters = [];
handles.measured_parameters = [];

% Removes the previous fitted ellipses and other results:
cla(handles.axes_ratio)

% Disables the appropriate puttons and options:
set(handles.push_calc,'Enable','off');
set(handles.push_start,'Enable','off');
set(handles.push_stop,'Enable','off');

% --- Maps the ratio of min and max component of the electric field with 
% --- the HWP and QWP combinations that provide linear and circular
% --- polarized light
function [handles] = image_ratio(handles)
% handles    structure with handles and user data (see GUIDATA)

% Calculates the ratio:
E_min = handles.fitted_parameters.E_min;
E_max = handles.fitted_parameters.E_max;
r = E_min./E_max;

% Plots the result:
axes(handles.axes_ratio)
x = handles.QWP_values(:,1)';
y = handles.HWP_values(1,:);
imagesc(x,y,r',[0 1])
axis equal; axis tight;
xlabel('QWP [deg]'); ylabel('HWP [deg]');
set(gca,'XAxisLocation','top','XTick',x(1:2:end),'YTick',y(1:2:end))
colormap('jet'); colorbar;

% Plots the data for circular and linear polarized light:
hold on;
plot(handles.CP.QWP,handles.CP.HWP,'ko', ...
    handles.LP.QWP, handles.LP.HWP, 'kx','linewidth',2);
axis([-5 185 -5 95])
set(gca,'XTick',0:20:180,'YTick',0:20:90)
legend('Circular polarization','Linear polarization');
hold off;

% -------------------------------------------------------------------------
% ------------------- Functions that are used for fitting -----------------

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
% The angular difference between the reflection axis and the HWP (at 0 deg)
handles.fitted_parameters.HWP_diff = P(3);
% The angular difference between the reflection axis and the QWP (at 0 deg)
handles.fitted_parameters.QWP_diff = P(4);
% The angular difference between the reflection axis and the LP (at 0 deg):
handles.fitted_parameters.LP_diff = P(5);
% The induced phase along the reflection axis:
handles.fitted_parameters.delta = P(6);

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
% ---------- Connecting and disconnecting the WPs and microscope: ---------

% --- Connects to the motors and microscope:
function [handles, connected] = connect_all(handles)
% handles    structure with handles and user data (see GUIDATA)
% connected  boolean stating if the connection was successfull

% Initializes the motors:
handles.HWP.StartCtrl; %m
handles.QWP.StartCtrl; %m
% Gives the GUI time to update:
pause(0.5)

% Homes the motors:
c(1) = handles.HWP.MoveHome(0,1); %m
c(2) = handles.QWP.MoveHome(0,1); %m

% Gives a message and stops if the motors were not connected:
if sum(c) ~= 0
    handles = disconnect_all(handles);
    warning('Could not connect to the motors.');
    connected = 0;
    return
end

% Tries to connect to the microscope:
try
    fopen(handles.microscope); %M
catch
    handles = disconnect_all(handles);
    warning('Could not connect to the microscope.');
    connected = 0;
    return
end

% Receives an initial message where the microscope welcomes you and
% mentions your client name:
t = 0;
while handles.microscope.BytesAvailable == 0 && t<100 %M
    pause(0.05)
    handles.microscope.BytesAvailable %M
    t = t+1;
end
% If the message never came:
if t == 100
    handles = disconnect_all(handles);
    warning('Could not connect to the microscope.');
    connected = 0;
    return
end
m = fread(handles.microscope,handles.microscope.BytesAvailable); %M
m = char(m');
disp(m)
% The client name and the text that every message to the microscope should
% start with:
handles.client = m(18:end-1); %M
handles.code_pre = ['/cli:' handles.client ' /app:matrix /cmd:']; %M
connected = 1;


% --- Disconnects the motors and microscope:
function [handles] = disconnect_all(handles)
% handles    structure with handles and user data (see GUIDATA)

% Stops the motors:
handles.HWP.StopCtrl; %m
handles.QWP.StopCtrl; %m

% Disables the communication with the microscope:
fclose(handles.microscope); %M

% -------------------------------------------------------------------------
% ------------------- WHEN THE GUI IS CLOSED: -----------------------------

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Disconnects WPs and microscope if they are not already deactivated:
handles = disconnect_all(handles);
% Closes the figure containing the WP controls:
close(handles.fig);

% Hint: delete(hObject) closes the figure
delete(hObject);
