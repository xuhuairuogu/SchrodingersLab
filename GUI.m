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

% Last Modified by GUIDE v2.5 15-Jun-2016 20:48:01

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
function GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

% Choose default command line output for GUI
handles.output = hObject;

order =  str2double(handles.iParamEdit1.String);        % Read order
a = str2double(handles.aEdit.String);                   % Read parameter a
L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)   
xj = eval(handles.xjEdit.String);                       % Eval x-shift 
tj = eval(handles.tjEdit.String);                       % Eval t-shift
T = str2double(handles.iParamEdit2.String);             % Eval max time
[PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);  % Prepare initial DT picture

handles.PSI_anal = PSI;                                 % Save result to handles struct
handles.x_anal = x;
handles.t_anal = t;
guidata(hObject, handles);

densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);       % Plot density
PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));             % Calculate spectrum
fourierPlot(PSI_k', t, 7, 1, handles.axes5);               % Plot spectrum
axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2)))); % Prepare titles
axes(handles.axes5); title('Analytical');   % Prepare titles

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.main_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function dtEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to dtEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function NEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to NEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function tfEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to tfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function LEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to LEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function aEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to aEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in LRadioButton.
function LRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to LRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LRadioButton
    set(handles.LEdit, 'Enable', 'on')
    set(handles.aEdit, 'Enable', 'off')


% --- Executes on button press in aRadioButton.
function aRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to aRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of aRadioButton
    set(handles.LEdit, 'Enable', 'off')
    set(handles.aEdit, 'Enable', 'on')

% --- Executes during object creation, after setting all properties.
function potEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to potEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function parabolicEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to parabolicEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cubicRadioButton.
function cubicRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to cubicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cubicRadioButton
    set(handles.potEdit, 'Enable', 'on')
    set(handles.parabolicEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'off')


% --- Executes on button press in parabolicRadioButton.
function parabolicRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to parabolicRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parabolicRadioButton
    set(handles.potEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'off')
    set(handles.parabolicEdit, 'Enable', 'on')

% --- Executes during object creation, after setting all properties.
function vOtherEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to vOtherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vOtherRadioButton.
function vOtherRadioButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to vOtherRadioButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of vOtherRadioButton
    set(handles.potEdit, 'Enable', 'off')
    set(handles.vOtherEdit, 'Enable', 'on')
    set(handles.parabolicEdit, 'Enable', 'off')


% --- Executes during object creation, after setting all properties.
function iLabel1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iLabel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function iParamEdit1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function psi0Listbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in psi0Listbox.
function psi0Listbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to psi0Listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns psi0Listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from psi0Listbox
    selection=hObject.Value;
   
    if (selection ~= 1)
        handles.xjLabel.Visible = 'off';
        handles.tjLabel.Visible = 'off';
        handles.xjEdit.Visible = 'off';
        handles.tjEdit.Visible = 'off';
        handles.intensityEdit.Visible = 'off';
        handles.intensityCheck.Visible = 'off';
        handles.iCheck1.Visible = 'off';
        handles.decimalsEdit.Visible = 'off';
    end
    if (selection ~= 2)
        handles.iParam1.Visible = 'on';
        handles.iParamEdit1.Visible = 'on';
        handles.xjEdit.Visible = 'off';
        handles.tjEdit.Visible = 'off';
        handles.intensityEdit.Visible = 'off';
        handles.intensityCheck.Visible = 'off';
    end
    
    
    if (selection == 1)
        handles.iLabel1.String = 'psi0=A0+\Sum_j 2*Aj*cos(j*w*x)';
        handles.iLabel2.Visible = 'off';
        handles.iParam1.String = 'Order';
        handles.iParamEdit1.String = '2';
        handles.iParam2.Visible = 'on';
        handles.iParamEdit2.Visible = 'on';
        handles.iParam2.String = 'T';
        handles.iParamEdit2.String = '-10';
        handles.intensityEdit.Visible = 'on';
        handles.intensityCheck.Visible = 'on';
        handles.iCheck1.Visible = 'on';
        handles.decimalsEdit.Visible = 'on';
        
        handles.xjLabel.Visible = 'on';
        handles.tjLabel.Visible = 'on';
        handles.xjEdit.Visible = 'on';
        handles.tjEdit.Visible = 'on';
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        
        a = str2double(handles.aEdit.String);
        order =  str2double(handles.iParamEdit1.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k.', t, 7, 1, handles.axes5);
        axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
        axes(handles.axes5); title('Analytical');

    elseif (selection == 2)
        handles.iLabel1.String = 'psi0=A0+2*A1*cos(w*x)';
        handles.iParam1.Visible = 'off';
        handles.iParamEdit1.Visible = 'off';
        handles.iParam2.Visible = 'off';
        handles.iParamEdit2.Visible = 'off';
        handles.iCheck1.Visible = 'off';
        handles.decimalsEdit.Visible = 'off';
    
    elseif (selection == 3)
        handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)H_n(x/s)';
        handles.iLabel2.String = 'H_n = nth order Hermite polynomial';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 's';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '5';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'on';
        handles.iParamEdit2.Visible = 'on';
        handles.iCheck1.Visible = 'off';
        handles.decimalsEdit.Visible = 'off';
        handles.iParam1.Visible = 'on';
        handles.iParamEdit1.Visible = 'on';
    elseif (selection == 4)
        handles.iLabel1.String = 'psi0=exp(-x^2/2/s^2)J_n(x)';
        handles.iLabel2.String = 'J_n = nth order Bessel polynomial';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 's';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '10';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'on';
        handles.iParamEdit2.Visible = 'on';
        handles.iCheck1.Visible = 'off';
        handles.decimalsEdit.Visible = 'off';
        handles.iParam1.Visible = 'on';
        handles.iParamEdit1.Visible = 'on';
    elseif (selection == 5)
        handles.iLabel1.String = 'psi0=exp(a*x)*Ai(x)';
        handles.iLabel2.String = 'Ai(x) = Airy function';
        handles.iLabel3.Visible = 'off';
        handles.iParam1.String = 'a';
        handles.iParam2.String = 'n';
        handles.iParamEdit1.String = '0.1';
        handles.iParamEdit2.String = '0';
        handles.iParam2.Visible = 'off';
        handles.iParamEdit2.Visible = 'off';
        handles.iCheck1.Visible = 'off';
        handles.decimalsEdit.Visible = 'off';
        handles.iParam1.Visible = 'on';
        handles.iParamEdit1.Visible = 'on';
    end


function iParamEdit1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of handles.iParamEdit1 as text
%        str2double(get(hObject,'String')) returns contents of handles.iParamEdit1 as a double

selection=handles.psi0Listbox.Value;
if (selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k', t, 7, 1, handles.axes5);
		axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
		axes(handles.axes5); title('Analytical');
end

function iParamEdit2_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of handles.iParamEdit2 as text
%        str2double(get(hObject,'String')) returns contents of handles.iParamEdit2 as a double
selection=handles.psi0Listbox.Value;
if (selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k', t, 7, 1, handles.axes5);
		axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
		axes(handles.axes5); title('Analytical');
end

% --- Executes during object creation, after setting all properties.
function iParamEdit2_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to handles.iParamEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function orderBox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to orderBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function aLMultEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to aLMultEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in iCheck1.
function iCheck1_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to iCheck1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1
    handles.decimalsEdit.Enable = 'on';
else
    handles.decimalsEdit.Enable = 'off';
end
    
% Hint: get(hObject,'Value') returns toggle state of iCheck1


% --- Executes on button press in analSwitch.
function analSwitch_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to analSwitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = allchild(handles.axes4);
% disp(h.FaceLighting)
% disp(h.AmbientStrength)
% disp(h.DiffuseStrength) 
% disp(h.SpecularStrength) 
% disp(h.SpecularExponent)
% disp(h.BackFaceLighting)
if (size(h) == 1); % Density Plot mode
    axes(handles.axes4);
    colorbar off;
    view(-62,42)
    shading interp
    lightangle(-40,50);
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.5;
    h.SpecularExponent = 3;
    h.BackFaceLighting = 'unlit';
    grid off;
else
    axes(handles.axes4);
    colorbar;
    delete(h(1));
    h = h(2);
    view([0,0,90])
    h.FaceLighting = 'flat';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.6;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 10;
    h.BackFaceLighting = 'reverselit';
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = allchild(handles.axes2);
if (size(h) == 1); % Density Plot mode
    axes(handles.axes2);
    colorbar off;
    view(-62,42)
    shading interp
    lightangle(-40,50);
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.5;
    h.SpecularExponent = 3;
    h.BackFaceLighting = 'unlit';
    grid off;
else
    axes(handles.axes2);
    colorbar;
    delete(h(1));
    h = h(2);
    view([0,0,90])
    h.FaceLighting = 'flat';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.6;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 10;
    h.BackFaceLighting = 'reverselit';
end

function aEdit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to aEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aEdit as text
%        str2double(get(hObject,'String')) returns contents of aEdit as a double
selection=handles.psi0Listbox.Value;
if (selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k', t, 7, 1, handles.axes5);
		axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
		axes(handles.axes5); title('Analytical');
        k = str2double(handles.iParamEdit1.String);
        a = str2double(hObject.String);
        handles.intensityEdit.String = num2str(peakPredict(a, k)^2);
end

% --- Executes on button press in undockDT1.
function undockDT1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockDT1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
h = allchild(handles.axes4);
copyobj(handles.axes4,f); colormap('jet'); colorbar;
%title(' ');
if (size(h) == 1); % Density Plot mode
    colorbar
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
else
    colorbar off;
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
end


% --- Executes on button press in undockDT2.
function undockDT2_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockDT2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
leg = legend(handles.axes5);
copyobj([leg, handles.axes5],f);
%title(' ');
set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);


% --- Executes on button press in undockNum1.
function undockNum1_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockNum1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
h = allchild(handles.axes2);
copyobj(handles.axes2,f); colormap('jet'); colorbar;
%title(' ');
if (size(h) == 1); % Density Plot mode
    colorbar
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
else
    colorbar off;
    set(gca, 'Position', [0.1300    0.1100    0.7750    0.8150]);
end

% --- Executes on button press in undockNum2.
function undockNum2_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to undockNum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f = figure;
leg = legend(handles.axes3);
copyobj([leg, handles.axes3],f);
title(' ');
% if handles.SpectrumPubMode
%     set(gcf,'color','w');
%     f.Position = [488.2000  541.8000  560.0000  220.0000];
%     set(gca, 'Position', [0.1    0.22    0.65    0.7]);
% else
%     set(gca, 'Position', [0.1300    0.1100    0.66    0.8150]);
% end
%ax.Position = [0.1300    0.1948    0.7750    0.7302];
set(gca, 'Position', [0.1300    0.1100    0.66    0.8150]);

function xjEdit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to xjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xjEdit as text
%        str2double(get(hObject,'String')) returns contents of xjEdit as a double
selection=handles.psi0Listbox.Value;
if (selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k', t, 7, 1, handles.axes5);
		axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
		axes(handles.axes5); title('Analytical');
end

% --- Executes during object creation, after setting all properties.
function xjEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to xjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tjEdit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to tjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tjEdit as text
%        str2double(get(hObject,'String')) returns contents of tjEdit as a double
selection=handles.psi0Listbox.Value;
if (selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a);                                     %#ok<NASGU> % Read length L (will be used in shifts)    
        lambda = sqrt(8*a*(1-2*a));                             %#ok<NASGU> % Growth factor (used in shifts)
        Omega = 2*sqrt(1-2*a);                                  %#ok<NASGU> % Fundamental wavenumber (used in shifts)  
        xj = eval(handles.xjEdit.String);
        tj = eval(handles.tjEdit.String);
        T = str2double(handles.iParamEdit2.String);
        [PSI, x, t] = calcDarboux(order, a, 256, 300, T, xj, tj);
        
        handles.PSI_anal = PSI;                                 % Save result to handles struct
        handles.x_anal = x;
        handles.t_anal = t;
        guidata(hObject, handles);
        
        densityPlot(abs(PSI).^2, x, t, 1, 1, handles.axes4);
        PSI_k = log(abs(fft(PSI'))/length(PSI(1, :)));
        fourierPlot(PSI_k', t, 7, 1, handles.axes5);
		axes(handles.axes4); title(sprintf('Analytical. Max = %.3f', max(max(abs(PSI).^2))));
		axes(handles.axes5); title('Analytical');
end

% --- Executes during object creation, after setting all properties.
function tjEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to tjEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function decimalsEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to decimalsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Data', cell(6, 1));
set(hObject, 'RowName', {'A0', 'A1', 'A2', 'A3', 'A4', 'A5'});

% --- Executes on button press in pubSizeButton.
function pubSizeButton_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pubSizeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = handles.axes2;
axes(ax);
maxY = ax.YLim(2);
ylim([maxY/2-2, maxY/2+2]);
maxX = ax.XLim(2);
if maxX > 4;
    xlim([-4, 4]);
end
title('');

% --------------------------------------------------------------------
function menuPreferences_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to menuPreferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Preferences;

% --- Executes during object creation, after setting all properties.
function intensityEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to intensityEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function potMenu_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to potMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns potMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from potMenu
str = get(handles.potMenu, 'String');
val = get(handles.potMenu,'Value');
% Set current data to the selected data set.
switch str{val};
case 'g*|psi|^2 (cubic)' % User selects peaks.
    handles.potEdit.String = num2str(-1);
case 'alpha*x^2 (parabolic)' % User selects membrane.
    handles.potEdit.String = num2str(0.3);
end

% --- Executes during object creation, after setting all properties.
function potMenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to potMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = handles.axes3;
axes(ax);
h = allchild(ax);
if handles.SpectrumPubMode
CO =     [         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
    handles.SpectrumPubMode = 0;
    title('Numerical');
    width = 1;
    maxX = ax.XLim(2);
    minX = ax.XLim(1);
    xlim([0, minX + maxX]);
    ylim ([-10, 2]);
    ax.XTick = 0:(minX+maxX)/4:minX+maxX;
    ax.YTick = -10:2:2;
else
CO=[0.0314    0.1882    0.4196; 
    0         0.5       0.1412;  
    0.8941    0.1020    0.1098; 
    1.0000    0.4980    0     ; 
    0.3059    0.7020    0.8275;
    0.6824    0.0039    0.4941;
    0.4980         0         0];
    handles.SpectrumPubMode = 1;
    title('');
    width = 1.5;
    ylim([-8, 2]);
    maxX = ax.XLim(2);
    xlim([maxX/2-2, maxX/2+2]);
    ax.XTick = maxX/2-2:1:maxX/2+2;
    ax.YTick = -8:2:2;
end

nLines = length(h);
for i = 1:nLines
	h(i).Color = CO(nLines-i+1, :);
    h(i).LineWidth = width;
end

guidata(hObject,handles)

% --------------------------------------------------------------------
function vidSettings_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to vidSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
videoSettings


% --- Executes on button press in numericalMaxima.
function numericalMaxima_Callback(hObject, ~, handles)
% hObject    handle to numericalMaxima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maximaType = 1;     % This indicates function should use numerical data to compute maxima
guidata(hObject, handles);
maximaDisplay;

% --- Executes on button press in analyticalMaxima.
function analyticalMaxima_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to analyticalMaxima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maximaType = 2;         % This indicates function should use DT data to compute maxima
guidata(hObject, handles);
maximaDisplay;

% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN BUTTON CALL BACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in runButton.
function runButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    han    5dle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    dt = str2double(handles.dtEdit.String);                    % Read dt from edit
    Nx =  str2double(handles.NEdit.String);                    % read Nx from edit
    Tmax = str2double(handles.tfEdit.String);                  % Read Tmax from edit
    Nt = Tmax/dt;                                              % Calculate Nt
    mult = str2double(handles.aLMultEdit.String);              % Read mult from edit box
    
    if strcmp(handles.LEdit.Enable, 'on')
        L = str2double(handles.LEdit.String);
    else
        a = str2double(handles.aEdit.String);
        L = pi/sqrt(1-2*a)*mult;
    end
    
    xj = eval(handles.xjEdit.String);
    tj = eval(handles.tjEdit.String);
    
    dx = L/Nx;                             % Spatial step size
    x = (-Nx/2:1:Nx/2-1)'*dx;
    
    % Determine the selected data set.
    str = get(handles.potMenu, 'String');
    val = get(handles.potMenu,'Value');
    % Set current data to the selected data set.
    switch str{val};
    case 'g*|psi|^2 (cubic)' % User selects peaks.
        g = str2double(handles.potEdit.String);
        V = @(psi, x) (g*abs(psi).^2);
    case 'alpha*x^2 (parabolic)' % User selects membrane.
        alpha = str2double(handles.potEdit.String);
        V = @(psi, x) 1/2*alpha*x.^2;
    end
    
    psi0_selection = handles.psi0Listbox.Value;
    
    if (psi0_selection == 1)
        order =  str2double(handles.iParamEdit1.String);
        T = str2double(handles.iParamEdit2.String);
        Omega = 2*sqrt(1-2*a);
        [psi_dt, x_dt, ~] = calcDarboux(order, a, Nx, 0, T, xj, tj);
        fitP = genCoeff(x_dt, psi_dt, a, order, 1, handles.uipanel2);
        
        coeff = coeffvalues(fitP{1}) + 1i*coeffvalues(fitP{2});
        A = abs(coeff(2:end));          
                            
        phi_0 = angle(coeff(1));
        phi = zeros(1,order);
        for i=1:order
            phi(i) = angle(coeff(i+1));
            A(i) = A(i)*exp(1i*(phi(i) - phi_0));
        end
        
        if handles.iCheck1.Value == 1
            d = str2double(handles.decimalsEdit.String);
            Dr = d - ceil(log10(abs(real(A))));
            Di = d - ceil(log10(abs(imag(A))));
            D = zeros(1, order);
            for i=1:order
                if Dr > Di
                    D(i) = Di(i);
                else
                    D(i) = Dr(i);
                end
            end
            Ar = zeros(1, order); Ai = zeros(1, order);
            for i=1:order
                Ar(i) = round(real(A(i)), D(i));
                Ai(i) = round(imag(A(i)), D(i));
            end
            A = Ar + 1i*Ai;
        end   
%         format longe
%         disp(A.');
        
        A02 = 1;
        for i = 1:order
            A02 = A02 - 2*abs(A(i))^2;
        end
        A0 = sqrt(A02);

        psi_0 = A0;
        for i=1:order
            psi_0 = psi_0 + 2*A(i)*cos(i*Omega*x);
        end
%         norm = dx*(psi_0')*(psi_0)/L;
%         disp(['Norm is: ', num2str(norm)]);
        handles.uitable1.Data = [A0, A].';
        
    elseif (psi0_selection == 2)
        
        A = handles.uitable1.Data;
        if iscell(A)
            for i = 1:length(A)
                if isempty(A{i})
                    B(i) = 0;
                else
                    B(i) = A{i};
                end
            end
            A = B.';
        end
        A = A(2:end);
        
        A02 = 1;
        for i = 1:length(A)
            A02 = A02 - 2*abs(A(i))^2;
        end
        A0 = sqrt(A02);
        
        Omega = 2*sqrt(1-2*a);
        psi_0 = A0;
        for i=1:length(A)
            psi_0 = psi_0 + 2*A(i)*cos(i*Omega*x);
        end
%        norm = dx*(psi_0')*(psi_0)/L;
%         disp(['Norm is: ', num2str(norm)]);
        handles.uitable1.Data = [A0; A];
    elseif (psi0_selection == 3)
        s = str2double(handles.iParamEdit1.String);
        n = str2double(handles.iParamEdit2.String);
        psi_0 = exp(-x.^2/2/s^2).*hermiteH(n, x/s);
    elseif (psi0_selection == 4)
        s = str2double(handles.iParamEdit1.String);
        n = str2double(handles.iParamEdit2.String);
        psi_0 = exp(-x.^2/2/s^2).*besselj(n, x);
    elseif (psi0_selection == 5)
        a = str2double(handles.iParamEdit1.String);
        psi_0 = exp(a*x).*airy(x);
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
    [PSI, x, t] = solve(dt, Nx, Tmax, L, mult, V, psi_0, method);  

%    initialPlot(PSI(1, :), x, handles.axes1);
    densityPlot(abs(PSI).^2, x, t, ceil(Nt/1000), ceil(Nx/256), handles.axes2);
    PSI_k = log(abs(fft(PSI.'))/Nx);
    cla(handles.axes3);
    fourierPlot(PSI_k', t, 7, 1, handles.axes3);

    % Some special purpose functions
    % [~,~,shift] = ab(PSI, x, t, max(max(abs(psi).^2)), 3/8);      
    % energy(PSI, t, k2, Nx, V, dt);
    
    axes(handles.axes2); title(sprintf('Numerical. Max = %.3f', max(max(abs(PSI).^2))));
    axes(handles.axes3); title('Numerical'); 
    handles.PSI_num = PSI;           % Save result to handles structure.
    handles.x_num = x;
    handles.t_num = t;
    guidata(hObject,handles)
    
    %figure;
    %h = axes();
    %PSI_k = fftshift(log(abs(fft(PSI'))/Nx), 1);
    %densityPlot(PSI_k.', x, t, Nt/1000, 1, h); colormap('jet'); 
    
    handles.SpectrumPubMode = 0;
    handles.IntensityPubMode = 0;
    guidata(hObject,handles)
    


% --------------------------------------------------------------------
function menuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to menuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
