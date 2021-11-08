function fn_paramPreviewFileNumCheck(handles,filenames)

allStr = get(handles.fnSelectionMenu,'String');
currStr = allStr{get(handles.fnSelectionMenu,'Value')};
switch currStr
    case 'getTuning'
        if length(filenames) ~=1 
            msgbox('WARNING -- Please select only one file for getTuning');
        end
        %filename = filenames{1};  
end

end