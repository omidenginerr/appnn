function varargout = appnn_gui(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @appnn_gui_OpeningFcn, ...
                       'gui_OutputFcn',  @appnn_gui_OutputFcn, ...
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
    

function appnn_gui_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject, handles);

function varargout = appnn_gui_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;


function input_sequences_Callback(hObject, eventdata, handles)


function input_sequences_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), ...
                       get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
function text_output_Callback(hObject, eventdata, handles)



function text_output_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    
function sequence_list_Callback(hObject, eventdata, handles)
    global sequences residuewise;
    index = get(handles.sequence_list,'value');
    if not(isempty(index))
        x = 1:length(sequences{index});
        plot(x,residuewise{index});
        ylim([-0.5,1.5]);
        xlim([1,length(sequences{index})]);
        line('XData', [1,length(sequences{index})], 'YData', [0.5 0.5], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color','m');
    end

function sequence_list_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

function run_bttn_Callback(hObject, eventdata, handles)
    global sequences residuewise ids;
    data = get(handles.input_sequences,'String');
    data = cellstr(data);
    try
        [ids,sequences] = cell2seq(data);
    catch
        msgbox('Please insert sequences in fasta format')
    end
    [results,hotspots,residuewise] = appnn(sequences);    
    textoutput(ids,sequences,results,hotspots,handles);
    set(handles.sequence_list,'String',ids)
    set(handles.sequence_list,'value',1)
    sequence_list_Callback(hObject, eventdata, handles)

function reset_bttn_Callback(hObject, eventdata, handles)
    set(handles.text_output,'String','');
    set(handles.input_sequences,'String','');
    set(handles.sequence_list,'String','');
    cla(handles.chart)

function open_bttn_Callback(hObject, eventdata, handles)
    [stmfile, stmpath] = uigetfile({'*.fasta','Fasta file'});
    fid = fopen(fullfile(stmpath,stmfile),'r');
    if fid ~= -1
        text = textscan(fid,'%s');
        text = char(text{1});
        set(handles.input_sequences,'String',text);
        fclose(fid);
    end
    

function export_Callback(hObject, eventdata, handles)
    global residuewise ids;
    [stmfile, stmpath] = uiputfile({'*.csv'});
    fid = fopen(fullfile(stmpath,stmfile),'w');
    if fid ~= -1
        text = ids{1};
        for i = 2:length(ids)
           text = sprintf('%s,%s',text,ids{i});
        end
        fprintf(fid,sprintf('%s\n',text));
        for i = 1:max(cellfun(@length,residuewise))
            text = '';
            for j = 1:length(ids)
                if length(residuewise{j})>=i
                    text = sprintf('%s,%0.4f',text,residuewise{j}(i));
                else
                    text = sprintf('%s,',text);
                end
            end
            text = text(2:length(text));
            fprintf(fid,sprintf('%s\n',text));
        end
        fclose(fid);
    end
    
function [identifiers,sequences] = cell2seq(data)
    sequences = {};
    identifiers = {};
    isid = @(x)(x(1)=='>');
    map = cellfun(isid,data);
    sequence = '';
    for i = 1:length(map)
        if map(i)
            if length(sequence)>0
                sequences = vertcat(sequences,{sequence});
            end
            identifiers = vertcat(identifiers,strrep(data{i},'>',''));
            sequence = '';
        else
            sequence = horzcat(sequence,data{i});
        end
    end
    sequences = vertcat(sequences,{sequence});

function textoutput(ids, seqs,results,hotspots,handles)
    text = '';
    for i = 1:length(seqs)
        text = sprintf('%sSequence: %s\nAmyloidogenic: %0.4f\nHotspots:\n',text,ids{i},results(i));
        htspts = hotspots{i};
        if length(htspts)>1
            for j = 1:length(htspts(:,1))
                s = htspts(j,1);
                e = htspts(j,2);
                text = sprintf('%s\t(%i) %s (%i)\n',text,s,seqs{i}(s:e),e);
            end
        else
            text = sprintf('%s\tNone\n',text);
        end
        text = sprintf('%s\n\n',text);
    end
    set(handles.text_output,'String',text);
