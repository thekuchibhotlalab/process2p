function fn_paramPreviewFileNumCheck(handles,filenames)

allStr = get(handles.fnSelectionMenu,'String');
currStr = allStr{get(handles.fnSelectionMenu,'Value')};
switch currStr
    case 'getTuning'
        if length(filenames) ~=1 
            disp('ERROR: Please select only one file');
        end
        %filename = filenames{1};  
end

end