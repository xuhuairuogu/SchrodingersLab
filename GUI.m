function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 16-Oct-2015 17:43:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
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


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function dtEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function NEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function tfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function LEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function aEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LRadioButton.
function LRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to LRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LRadioButton
    set(handles.LEdit, 'Enable', 'on')
    set(handles.aEdit, 'Enable', 'off')


% --- Executes on button press in aRadioButton.
function aRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to aRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aRadioButton
    set(handles.LEdit, 'Enable', 'off')
    set(handles.aEdit, 'Enable', 'on')

% --- Executes during object creation, after setting all properties.
function cubicEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cubicEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function parabolicEdit_Callback(hObject, eventdata, handles)
% hObject    handle to parabolicEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of parabolicEdit as text
%        str2double(get(hObject,'String')) returns contents of parabolicEdit as a double


% --- Executes during object creation, after setting all properties.
function parabolicEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parabolicEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cubicRadioButton.
function cubicRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to cubicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cubicRadioButton
    set(handles.cubicEdit, 'Enable', 'on')
    set(handles.parabolicEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'off')


% --- Executes on button press in parabolicRadioButton.
function parabolicRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to parabolicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parabolicRadioButton
    set(handles.cubicEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'off')
    set(handles.parabolicEdit, 'Enable', 'on')

% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
% hObject    han    5dle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    dt = eval(handles.dtEdit.String);
    N = eval(handles.NEdit.String);
    tf = eval(handles.tfEdit.String);
    
    if strcmp(handles.LEdit.Enable, 'on')
        L = eval(handles.LEdit.String);
    else
        a = eval(handles.aEdit.String);
        L = pi/sqrt(1-2*a);
    end
    
    if strcmp(handles.cubicEdit.Enable, 'on')
        V = sprintf('%s*abs(psi).^2', handles.cubicEdit.String);
    else
        V = sprintf('1/2*%s^2*x.^2', handles.parabolicEdit.String);
    end
    
    psi0_selection = handles.psi0Listbox.Value;
    
    if (psi0_selection == 1)
        A1 = eval(handles.iParamEdit1.String);
        A0 = sqrt(1-2*A1^2);
        omega = 2*sqrt(1-2*a);
        psi_0 = sprintf('%.16f+2*%.30f*cos(%.16f*x)', A0, A1, omega) ;
    elseif (psi0_selection == 2)
        s = eval(handles.iParamEdit1.String);
        n = eval(handles.iParamEdit2.String);
        psi_0 = sprintf('exp(-x.^2/2/%.16f^2).*hermiteH(%.16f, x/%.16f)', s, n, s) ;
    elseif (psi0_selection == 3)
        s = eval(handles.iParamEdit1.String);
        n = eval(handles.iParamEdit2.String);
        psi_0 = sprintf('exp(-x.^2/2/%.16f^2).*besselj(%.16f, x)', s, n) ;
    elseif (psi0_selection == 4)
        a = eval(handles.iParamEdit1.String);
        psi_0 = sprintf('exp(%.16f*x).*airy( x)', a) ;
    end
    
    orderSelection = handles.orderBox.Value;
    if orderSelection == 1
        method = 'T2';
    elseif orderSelection == 2
        method = 'T4_NS';
    elseif orderSelection == 3
        method = 'T4';
    elseif orderSelection == 4
        method = 'T6_NS';
    elseif orderSelection == 5
        method = 'T6';
    elseif orderSelection == 6
        method = 'T8_NS';    
    elseif orderSelection == 7
        method = 'T8';
    end
        
    solve(dt, N, tf, L, V, psi_0, method, handles)

function vOtherEdit_Callback(hObject, eventdata, handles)
% hObject    handle to vOtherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vOtherEdit as text
%        str2double(get(hObject,'String')) returns contents of vOtherEdit as a double


% --- Executes during object creation, after setting all properties.
function vOtherEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vOtherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vOtherRadioButton.
function vOtherRadioButton_Callback(hObject, eventdata, handles)
% hObject    handle to vOtherRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vOtherRadioButton
    set(handles.cubicEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'on')
    set(handles.parabolicEdit, 'Enable', 'off')


% --- Executes during object creation, after setting all properties.
function iLabel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handles.iLabel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function iParamEdit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handles.iParamEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tfEdit_Callback(hObject, eventdata, handles)
% hObject    handle to tfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tfEdit as text
%        str2double(get(hObject,'String')) returns contents of tfEdit as a double

% --- Executes during object creation, after setting all properties.
function psi0Listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in psi0Listbox.
function psi0Listbox_Callback(hObject, eventdata, handles)
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns psi0Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from psi0Listbox
    selection=hObject.Value;
    if (selection == 1)
        handles.iLabel1.String = 'psi0=A0+2*A1*cos(Om*x)';
        handles.iLabel2.String = 'A0 = sqrt(1-2*A1^2)';
        handles.iLabel3.Visible = 'on';
        handles.iLabel3.String = 'Om= 2*sqrt(1-2*a)';
        handles.iParam1.String = 'A1';
        handles.iParamEdit1.String = '1-4';
        handles.iParam2.Visible = 'off';
        handles.iParamEdit2.Visible = 'off';
    elseif (selection == 2)
        handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)H_n(x/s)';
        handles.iLabel2.String = 'H_n = nth order Hermite polynomial';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 's';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '5';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'on';
        handles.iParamEdit2.Visible = 'on';
    elseif (selection == 3)
        handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)J_n(x)';
        handles.iLabel2.String = 'J_n = nth order Bessel polynomial';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 's';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '10';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'on';
        handles.iParamEdit2.Visible = 'on';
    elseif (selection == 4)
        handles.iLabel1.String = 'psi0=exp(a*x)*Ai(x)';
        handles.iLabel2.String = 'Ai(x) = Airy function';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 'a';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '0.1';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'off';
        handles.iParamEdit2.Visible = 'off';
end


function iParamEdit1_Callback(hObject, eventdata, handles)
% hObject    handle to handles.iParamEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of handles.iParamEdit1 as text
%        str2double(get(hObject,'String')) returns contents of handles.iParamEdit1 as a double



function iParamEdit2_Callback(hObject, eventdata, handles)
% hObject    handle to handles.iParamEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of handles.iParamEdit2 as text
%        str2double(get(hObject,'String')) returns contents of handles.iParamEdit2 as a double


% --- Executes during object creation, after setting all properties.
function iParamEdit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handles.iParamEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function orderBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
