function varargout = analProp(varargin)
% ANALPROP MATLAB code for analProp.fig
%      ANALPROP, by itself, creates a new ANALPROP or raises the existing
%      singleton*.
%
%      H = ANALPROP returns the handle to a new ANALPROP or the handle to
%      the existing singleton*.
%
%      ANALPROP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALPROP.M with the given input arguments.
%
%      ANALPROP('Property','Value',...) creates a new ANALPROP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analProp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analProp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analProp

% Last Modified by GUIDE v2.5 30-Jun-2016 21:41:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analProp_OpeningFcn, ...
                   'gui_OutputFcn',  @analProp_OutputFcn, ...
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


% --- Executes just before analProp is made visible.
function analProp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analProp (see VARARGIN)

% Choose default command line output for analProp
handles.output = hObject;
% get the handle of Gui1
h = findobj('Tag','main_GUI');
% if exists (not empty)
if ~isempty(h)
% get handles struct from main GUI
mainHandles = guidata(h);
    if length(mainHandles.a) == 1
        a(1) = mainHandles.a(1);
        for k = 1:mainHandles.order
            a(k) = k^2*(a(1)-1/2)+1/2;
        end
    else
        a = mainHandles.a;
    end
end
%handles.uitable1.Data = cell(mainHandles.order);

Omega = 2*sqrt(1-2*a);
lambda = sqrt(8*a.*(1-2*a));
L_pi = 1./sqrt(1-2*a);
L = pi./sqrt(1-2*a);
disp(size(a));

data = [a; Omega; lambda; L; L_pi];
handles.uitable1.Data = data';
    

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analProp wait for user response (see UIRESUME)
% uiwait(handles.analFig);


% --- Outputs from this function are returned to the command line.
function varargout = analProp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
