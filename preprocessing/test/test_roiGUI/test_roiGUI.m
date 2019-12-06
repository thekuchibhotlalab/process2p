screenSize = get(0,'screensize'); 
f_selection = figure; set(gcf, 'Position',  [100, 100, 200, screenSize(4)- 200]);

figureH = screenSize(4)-200;
items = 20;
itemH = figureH/(2*items+3);
itemHPos = linspace(0,figureH, 2*items+4);

buttonHPos = itemHPos(2);
barHPos = itemHPos(4:2:end);


%% get all the txts
global tempIsCell
global roi_redrawn
roi_redrawn = {[],[]};
nCell = 1;
startLoc = 10;
leng = 40;

for k = 1:length(barHPos)-1
c = uicontrol(f_selection,'Style','text'); c.Position = [startLoc barHPos(end-k) leng itemH];c.String = ['day' int2str(k)];
end
%c = uicontrol(f,'Style','pushbutton');c.Position = [startLoc buttonHPos leng itemH];c.String = 'Next';

%% get all the editable 
startLoc = 60;
leng = 60;
edit_box = {};
for k = 1:length(barHPos)-1
c = uicontrol(f_selection,'Style','edit'); c.Position = [startLoc barHPos(end-k) leng itemH];c.String = '1';
edit_box{k} = c;
end
c = uicontrol(f_selection,'Style','pushbutton','CallBack',{@continue_callback,f_selection, edit_box});c.Position = [startLoc buttonHPos leng itemH];c.String = 'next cell';

%% get all the redraw
startLoc = 130;
leng = 60;
redraw_button = {};
for k = 1:length(barHPos)-1
c = uicontrol(f_selection,'Style','togglebutton','CallBack',{@redraw_callback, k,nCell}); c.Position = [startLoc barHPos(end-k) leng itemH];c.String = 'redraw';
redraw_button{k} = c;
end
c = uicontrol(f_selection,'Style','pushbutton','CallBack',{@reject_callback,f_selection,items});c.Position = [startLoc buttonHPos leng itemH];c.String = 'reject all';


function  reject_callback(hObject,eventdata,f_selection,items)
    global tempIsCell 
    tempIsCell = zeros(1,items);
    disp(tempIsCell)
    disp('reject')
    close(f_selection)
end


function  continue_callback(hObject,eventdata,f_selection, edit_box)
    global tempIsCell
    tempIsCell = [];
    for i = 1:length(edit_box)
        tempStr = get(edit_box{i}, 'String');
        tempIsCell(i) = str2double(tempStr);
    end
    disp('continue')
    disp(tempIsCell)
    close(f_selection)
end

function redraw_callback(hObject,eventdata,nDay,nCell)
    global roi_redrawn
    button_state = get(hObject,'Value');
    if button_state % if triggered from off --> on, redraw roi
        h_redraw = figure;
        imagesc(rand(10,10));
        h_roi = imfreehand;
        roiMask = h_roi.createMask;
        [roix, roiy] = find(roiMask);
        roixy = [roix,roiy];
        disp(nDay)
        uiwait(h_redraw);
        
        tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        if sum(existCellDay) >= 1 % if the roi has been redrawn, repalce it
            %idx = find(existCellDay==1);
            disp('repeated roi in:')
            disp(roi_redrawn{existCellDay,1})
            roi_redrawn{existCellDay,2} = roixy;
        else % if roi not redrawn, create it
            roi_redrawn = [roi_redrawn;{tempCellDay,roixy}];
        end
        
    else % if triggered from on --> off, delete redrawn roi
        tempCellDay = [nCell nDay];tempC = cell(size(roi_redrawn,1),1); tempC(:) = {tempCellDay};
        existCellDay = cellfun(@isequal,roi_redrawn(:,1),tempC,'UniformOutput',true);
        if sum(existCellDay) >= 1 % if the roi has been redrawn, delete it
            roi_redrawn(existCellDay,:) = [];
        end
    end
end