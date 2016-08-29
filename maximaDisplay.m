function varargout = maximaDisplay(varargin)
% MAXIMADISPLAY MATLAB code for maximaDisplay.fig
%      MAXIMADISPLAY, by itself, creates a new MAXIMADISPLAY or raises the existing
%      singleton*.
%
%      H = MAXIMADISPLAY returns the handle to a new MAXIMADISPLAY or the handle to
%      the existing singleton*.
%
%      MAXIMADISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAXIMADISPLAY.M with the given input arguments.
%
%      MAXIMADISPLAY('Property','Value',...) creates a new MAXIMADISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before maximaDisplay_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to maximaDisplay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help maximaDisplay

% Last Modified by GUIDE v2.5 23-Jun-2016 22:10:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @maximaDisplay_OpeningFcn, ...
                   'gui_OutputFcn',  @maximaDisplay_OutputFcn, ...
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


% --- Executes just before maximaDisplay is made visible.
function maximaDisplay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to maximaDisplay (see VARARGIN)

% Choose default command line output for maximaDisplay
handles.output = hObject;

 % get the handle of Gui1
 h = findobj('Tag','main_GUI');
 % if exists (not empty)
 if ~isempty(h)
    % get handles struct from main GUI
    mainHandles = guidata(h);
    if mainHandles.maximaType == 1
        k = any(strcmp('PSI_num',fieldnames(mainHandles))); % Check if numerical data exists in memory
        if ~k % If it doesn't, display error
            h = errordlg('No numerical data exists in memory, cannot compute maxima.');
        end
        % maybe you want to set the text in Gui2 with that from Gui1
        PSI = mainHandles.PSI_num;
        x = mainHandles.x_num;
        t = mainHandles.t_num;
    elseif mainHandles.maximaType == 2
        PSI = mainHandles.PSI_anal;
        x = mainHandles.x_anal;
        t = mainHandles.t_anal;
    end
    
    maxima = maxima(PSI, x, t);
    
    handles.uitable1.Data = maxima;
 end

 handles.maxima = maxima;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes maximaDisplay wait for user response (see UIRESUME)
% uiwait(handles.maxima_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = maximaDisplay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

lim = str2double(hObject.String);
k = size(handles.maxima); k = k(1);

maxima = handles.maxima;
I = maxima(:, 3) < lim;
maxima(I, :) = [];

handles.uitable1.Data = maxima;


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
