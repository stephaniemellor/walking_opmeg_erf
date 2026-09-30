
% --- Preamble ----
addpath('C:\Users\Stephanie Mellor\Documents\GitHub\BrewerMap')
colormap123 = colormap(flipud(brewermap(64,'RdBu')));

results_dir = 'C:\Users\Stephanie Mellor\Documents\Data\Auditory\test_revision3_data\results';
sub = {'sub-001', 'sub-002', 'sub-003'};
final_fileloc = 'C:\Users\Stephanie Mellor\Documents\ReadingWriting\Writing\Self\neuro1_walking\third_submission\figures\';

%% Evoked response figures
clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration ---
numRows = 4;
numCols = 6;
task = {'seatedClosed', 'seatedOpen', 'walkingClosed', 'walkingOpen'};
row_folder_name = {'no_amm', 'hfc', 'amm_spatial', 'amm'};
colorbar_filepath = fullfile(final_fileloc, 'colourbars');

% --- Layout Geometry ---
figWidth = 1534; % Figure width in pixels
figHeight = 700; % Figure height in pixels

% Margins and spacing (normalized 0 to 1)
leftMargin   = 0.15;
rightMargin  = 0.07;
topMargin    = 0.05;
bottomMargin = 0.15;

colSpacing = 0.04; % Space between columns
rowSpacing = 0.04; % Space between rows

% Calculate width and height of individual panels (normalised)
panelWidth  = (1 - leftMargin - rightMargin - (numCols - 1) * colSpacing) / numCols;
panelHeight = (1 - topMargin - bottomMargin - (numRows - 1) * rowSpacing) / numRows;

% --- Pre-calculate Y-Coordinates for Top-Alignment ---
% To keep the grid perfectly aligned, we establish a global top Y-coordinate for each row.
rowYTop = zeros(numRows, 1);
rowYTop(1) = 1 - topMargin;

for r = 2:numRows
    rowYTop(r) = 1 - topMargin - (r-1)*(panelHeight + rowSpacing);
end

% Define the filenames of your 24 sub-panels (Row x Column)
% Adjust these paths/names to match your actual files
imageFiles = cell(numRows, numCols);

% Loop over participants
for ss = 1:length(sub)

    % Loop over recordings
    for tt = 1:length(task)

        % Loop over runs
        if strcmp(sub{ss}, 'sub-002') && startsWith(task{tt}, 'walking')
            run = {'_run-001', '_run-002'};
        else
            run = {''};
        end

        for rr = 1:length(run)
            for r = 1:numRows
                imageFiles{r, 1} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_average_time_series.fig', sub{ss}, task{tt}, run{rr})); 
                imageFiles{r, 2} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_t_stat_time_series.fig', sub{ss}, task{tt}, run{rr})); 
                imageFiles{r, 3} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_anti_averaging_topo_80_to_120_ms_anysig.fig', sub{ss}, task{tt}, run{rr})); % Ultimately add _test to end of filename
                imageFiles{r, 4} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_ROI.fig', sub{ss}, task{tt}, run{rr}));
                imageFiles{r, 5} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_min_norm_pow_left.fig', sub{ss}, task{tt}, run{rr}));
                imageFiles{r, 6} = fullfile(results_dir, sub{ss}, row_folder_name{r}, sprintf('%s_task-%s%s_meg_min_norm_pow_right.fig', sub{ss}, task{tt}, run{rr}));
            end

            f1 = figure('Position', [0, 50, figWidth, figHeight], 'Color', 'w');
            
            % --- Main Plotting Loop ---
            ax = gobjects(numRows, numCols);
            for r = 1:numRows
                for c = 1:numCols
                    
                    % 2. Load and Crop Image
                    if isfile(imageFiles{r, c})
                        f = openfig(imageFiles{r,c}, 'invisible');
                        a = f.Children;

                        % Calculate Position
                        xPos = leftMargin + (c - 1) * (panelWidth + colSpacing);
                        if c == 3
                            xPos = xPos-0.015;
                        end
                        % Align the tops of the plots; the extra bottom height will extend downward
                        yPos = rowYTop(r) - panelHeight; 
                        
                        if c == 3 || c == 4 || c == 5 || c == 6
                            ax(r,c) = copyobj(a(end), f1);
                            if c == 3
                                colormap(ax(r,c), colormap123);
                                ax(r,c).CameraViewAngle = 4.67;
                                if r == numRows
                                    cb = colorbar(ax(r,c), "southoutside");
                                    ylabel(cb, 'B (fT)')
                                end
                            else
                                colormap(ax(r,c), 'hot');
                            end
                            if c == 5 || c == 6
                                if length(a) == 3
                                    ax(r,c).CameraViewAngle = 4.0797;
                                else
                                    if ss == 1
                                        ax(r,c).CameraViewAngle = 6.89;
                                    else
                                        if r == 2
                                            ax(r,c).CameraViewAngle = 4.12;
                                        else
                                            ax(r,c).CameraViewAngle = 6.40;
                                        end
                                    end
                                end
                            end
                        else
                            ax(r, c) = copyobj(a, f1);
                        end
                        ax(r,c).Position = [xPos, yPos, panelWidth, panelHeight];
                        set(ax(r,c), 'FontSize', 16);
                        if r == numRows && c == 3
                            pos = cb.Position;
                            cb.Position(3) = 0.8*panelWidth;
                            cb.Position(4) = 1.2*cb.Position(4);
                            cb.Position(1) = pos(1) + pos(3)/2 - cb.Position(3)/2;
                            cb.Position(2) = ax(r,c).Position(2) - cb.Position(4) - 0.02;
                            cb.FontSize = 16;
                        end

                        axis tight
                        
                        ylabel(ax(r,c), '');
                        if r < numRows
                            set(ax(r,c), 'XTickLabel', []);
                            xlabel(ax(r,c), '');
                        end


                    else
                        % If file is missing
                        continue;
                    end
                    axis off;
                end
            end
            
            % --- Add Overall Labels (using normalized annotations) ---

            figure(f1);
            
            % Y-Axis Labels
            % B (fT) centered vertically on Column 1
            yCenter_col1 = (ax(2,1).Position(2) + ax(2,1).Position(4) + ax(3,1).Position(2))/2;
            annotation('textarrow', [1 1]*leftMargin-colSpacing*1.2, [yCenter_col1 yCenter_col1], ...
                'String', 'B (fT)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
                'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
                'VerticalAlignment','middle');
            
            % t-stat centered vertically on Column 2
            annotation('textarrow', [1 1]*leftMargin+panelWidth+0.009, [yCenter_col1 yCenter_col1], ...
                'String', 't-stat', 'HeadStyle', 'none', 'LineStyle', 'none', ...
                'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
                'VerticalAlignment', 'middle');
            
            % Current (nAm) centered vertically on Column 4
            xPos_col4 = leftMargin + 3 * (panelWidth + colSpacing) - 0.035;
            annotation('textarrow', [xPos_col4 xPos_col4], [yCenter_col1 yCenter_col1], ...
                'String', 'Current (nAm)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
                'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
                'VerticalAlignment', 'middle');

            % Add colourbars
            cb2_ax = colorbar(ax(end,5));
            cb2_ax.Position = [ax(1,5).Position(1) + ax(1,5).Position(3)+0.005, ...
                yCenter_col1-panelHeight*2/2, 0.01, panelHeight*2];

            cb3_ax = colorbar(ax(end,6));
            cb3_ax.Position = [ax(1,6).Position(1) + ax(1,6).Position(3)+0.005, ...
                yCenter_col1-panelHeight*2/2, 0.01, panelHeight*2];

            annotation('textarrow', (cb3_ax.Position(1) + cb3_ax.Position(3) + 0.025)*[1, 1], ...
                (cb3_ax.Position(2) + cb3_ax.Position(4)/2)*[1 1], ...
                'String', 'Minimum Norm Power (a.u.)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
                'TextRotation', -90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 15, ...
                'VerticalAlignment', 'middle');
            
            % Add legend
            lgd = legend(ax(1, 4), {'Left', 'Right'});
            lgd.Position = [ax(1,4).Position(1)+ax(1,4).Position(3)-panelWidth*0.45, ...
                    ax(1,4).Position(2)+ax(1,4).Position(4)-0.04, lgd.Position(3), lgd.Position(4)];

            % Add Column Numbering Headers at the Top
            columnHeaders = {'A) i)', 'ii)', 'iii)', 'B) i)', 'ii)', 'iii)'};
            for c = 1:numCols
                xCenter = ax(1,c).Position(1) - 0.002;
                % Position the text 1.5% of figure height above the top of the first row of panels
                yPos_header = rowYTop(1) + 0.025; 
                
                annotation('textarrow', [xCenter xCenter], [yPos_header yPos_header], ...
                    'String', columnHeaders{c}, 'HeadStyle', 'none', 'LineStyle', 'none', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18);
            end

            % Add Wrapped Row Names on the Left
            % Nested cells control multi-line wrapping perfectly
            rowNames = {
                'No filter', ...
                'HFC', ...
                {'Spatial', 'only AMM'}, ...
                {'AMM with', 'temporal', 'extension'}
            };
            
            for r = 1:numRows
                boxWidth = leftMargin*0.65; 
                boxHeight = 0.12;
                boxX = 0;                  % Placed on the far-left edge
                rowYCenter = ax(r, 1).Position(2) + ax(r,1).Position(4)/2;
                boxY = rowYCenter - (boxHeight / 2);
                
                annotation('textbox', [boxX, boxY, boxWidth, boxHeight], ...
                    'String', rowNames{r}, ...
                    'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'middle', ...
                    'FontName', 'Arial', ...
                    'FontSize', 16, ...
                    'FontWeight', 'bold', ...
                    'FontAngle', 'italic');
            end
            
            % Save output
            exportgraphics(gcf, fullfile(final_fileloc, sprintf('combined_%s_%s%s.png', sub{ss}, task{tt}, run{rr})), 'Resolution', 300);
            close all;
        end
    end
end

%% MMN Figures

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration ---
numRows = 3;
numCols = 4;
open_closed = {'Closed', 'Open'};
spatial_filter_name = {'amm'};

% --- Layout Geometry ---
figWidth = 1250; 
figHeight = 650; 

% Margins and spacing (normalized 0 to 1)
leftMargin   = 0.04;
rightMargin  = 0.06;
topMargin    = 0.075;
bottomMargin = 0.1;

colSpacing = 0.14; % Space between columns
rowSpacing = 0.07; % Space between rows

% Standard width for all panels in a column
panelWidth = (1 - leftMargin - rightMargin - (numCols - 1) * colSpacing) / numCols;

panelHeaders = {
    '(1) i)', 'ii)', '(1) i)', 'ii)';
    '(2) i)', 'ii)', '(2) i)', 'ii)';
    '(3) i)', 'ii)', '(3) i)', 'ii)'
};

% --- Step 1: Pre-calculate Y-Coordinates for Top-Alignment ---
% Establishing a global top-alignment Y-coordinate for each row
rowYTop = zeros(numRows, 1);
rowYTop(1) = 1 - topMargin;

% Calculate reference height based on Column 2 (Line plot) standard crop
panelHeight = (1 - topMargin - bottomMargin - (numRows - 1) * rowSpacing) / numRows;

for r = 2:numRows
    rowYTop(r) = rowYTop(r-1) - panelHeight - rowSpacing;
end

for oc = 1:length(open_closed)
    for sf = 1:length(spatial_filter_name)        
    
        % Define the filenames of your 12 sub-panels (Row x Column)
        % Adjust these paths/names to match your actual files
        imageFiles = cell(numRows, numCols);
        for ss = 1:length(sub)
            if strcmp(sub{ss}, 'sub-002')
                run = '_run-002';
            else
                run = '';
            end
            imageFiles{ss, 1} = fullfile(results_dir, sub{ss}, spatial_filter_name{sf}, sprintf('%s_task-seated%s_meg_MMN_topography_allstandards.fig', sub{ss}, open_closed{oc})); 
            imageFiles{ss, 2} = fullfile(results_dir, sub{ss}, spatial_filter_name{sf}, sprintf('%s_task-seated%s_meg_MMN_sensor_level_allstandards.fig', sub{ss}, open_closed{oc})); 
            imageFiles{ss, 3} = fullfile(results_dir, sub{ss}, spatial_filter_name{sf}, sprintf('%s_task-walking%s%s_meg_MMN_topography_allstandards.fig', sub{ss}, open_closed{oc}, run)); 
            imageFiles{ss, 4} = fullfile(results_dir, sub{ss}, spatial_filter_name{sf}, sprintf('%s_task-walking%s%s_meg_MMN_sensor_level_allstandards.fig', sub{ss}, open_closed{oc}, run)); 
        end

        f1 = figure('Position', [100, 100, figWidth, figHeight], 'Color', 'w');

        % --- Step 2: Main Plotting Loop ---
        ax = gobjects(numRows, numCols);
        for r = 1:length(sub)
            for c = 1:numCols
                if ~isfile(imageFiles{r, c})
                    % Placeholder box if the file is missing
                    xPos = leftMargin + (c - 1) * (panelWidth + colSpacing);
                    yPos = rowYTop(r) - panelHeight;
                    ax(r,c) = axes('Position', [xPos, yPos, panelWidth, panelHeight]);
                    text(0.5, 0.5, sprintf('Row %d, Col %d', r, c), 'HorizontalAlignment', 'center');
                    xlim([0 1]); ylim([0 1]); box on;
                    continue;
                end

                f = openfig(imageFiles{r,c}, 'invisible');
                a = f.Children;
                a = a(end);

                if c == 1 || c == 3
                    comment = findall(f, 'Type', 'TextBox');
                    comment = comment.String;
                end

                if isa(a, 'matlab.graphics.layout.TiledChartLayout')
                    a = a.Children;
                end
                ax(r,c) = copyobj(a(end), f1);
                
                xPos = leftMargin + (c - 1) * (panelWidth + colSpacing);
                % Align top-edges; any extra label height extends downward cleanly
                yPos = rowYTop(r) - panelHeight; 
                
                ax(r,c).Position = [xPos, yPos, panelWidth, panelHeight];
                set(ax(r,c), 'FontSize', 16);
                colormap(ax(r,c), colormap123);

                % Add Panel Header Inside Axis (Top-Left)
                if c == 1 || c == 3
                    ax(r,c).CameraViewAngle = 5.41;
                    annotation(f1, 'textbox', [ax(r,c).Position(1)-0.02, ax(r,c).Position(2)-0.05, 0.6, 0.05], ...
                        'string', comment, 'FontSize', 14, 'EdgeColor', 'None');
                    text(ax(r,c), 0.0001, 1.1, panelHeaders{r, c}, 'Units', 'normalized', ...
                        'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold', ...
                        'FontAngle', 'italic', 'VerticalAlignment', 'top', 'Color', 'k');
                else
                    text(ax(r,c), -0.5, 1.1, panelHeaders{r, c}, 'Units', 'normalized', ...
                        'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold', ...
                        'FontAngle', 'italic', 'VerticalAlignment', 'top', 'Color', 'k');
                    ylabel(ax(r,c),'t-stat');
                    if r < numRows
                        set(ax(r,c), 'XTickLabel', []);
                        xlabel(ax(r,c), '');
                    end
                    xl = xlim(ax(r,c));
                    yl = ylim(ax(r,c));
                    axis tight
                    xlim(ax(r,c), xl);
                    ylim(ax(r,c), yl);
                end
                
            end
        end

        % Add legend
        lgd = legend(ax(1, 4), a(2).String);
        lgd.Position = [ax(1,4).Position(1)+ax(1,4).Position(3)-panelWidth*0.5, ...
                ax(1,4).Position(2)+ax(1,4).Position(4)-panelHeight*0.3, lgd.Position(3), lgd.Position(4)];

        % Place colourbars
        for r = 1:numRows
            for c = [1, 3]
                if r == 2
                    cb = colorbar(ax(r,c));
                    pos = cb.Position;
                    cb.Position = [pos(1)+0.02, pos(2)+pos(4)/2-0.3/2, pos(3), 0.3];
                    ylabel(cb, 't-stat (MMN peak)');
                end
            end
        end

        % --- Step 3: Add Faint Vertical Dividing Line ---
        % Calculate the exact midpoint coordinate between Column 2 and Column 3
        xCol2Right = leftMargin + 1 * (panelWidth + colSpacing) + panelWidth;
        xCol3Left  = leftMargin + 2 * (panelWidth + colSpacing);
        xDivider   = (xCol2Right + xCol3Left) / 2;
        
        yDividerBottom = bottomMargin - 0.03;
        yDividerTop    = 1 - topMargin + 0.01;
        
        annotation(f1, 'line', [xDivider xDivider], [yDividerBottom yDividerTop], ...
            'Color', [0.8 0.8 0.8], 'LineStyle', '-', 'LineWidth', 1.2);
        
        % --- Step 4: Add Group Headings ("A) Seated" and "B) Walking") ---
        yGroupHeader = ax(1,1).Position(2) + ax(1,1).Position(4) + 0.03;
        
        % Center of left half (Columns 1 & 2)
        xCenterLeft = (ax(1,1).Position(1) + ax(1,2).Position(1) + ax(1,2).Position(3)) / 2;
        annotation(f1, 'textbox', [xCenterLeft - 0.1, yGroupHeader, 0.2, 0.04], ...
            'String', 'A) Seated', 'FontName', 'Arial', 'FontSize', 18, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'EdgeColor', 'none');
        
        % Center of right half (Columns 3 & 4)
        xCol4Right = leftMargin + 3 * (panelWidth + colSpacing) + panelWidth;
        xCenterRight = (xCol3Left + xCol4Right) / 2;
        annotation(f1, 'textbox', [xCenterRight - 0.1, yGroupHeader, 0.2, 0.04], ...
            'String', 'B) Walking', 'FontName', 'Arial', 'FontSize', 18, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'EdgeColor', 'none');

        % Save output
        exportgraphics(f1, fullfile(final_fileloc, sprintf('mmn_combined_%s_%s_allstandards.png', open_closed{oc}, spatial_filter_name{sf})), 'Resolution', 300);
    end
end

%% Sensor positions figure

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration ---
numRows = 2;
numCols = 3;

% List of your 6 MATLAB figure (.fig) files
figFiles = cell(1,length(sub));
for ss = 1:length(sub)
    figFiles{ss} = fullfile(results_dir, sub{ss}, 'helmet.fig');
end

% Participant titles for each panel
partTitles = {
    'Participant 1', 'Participant 2', 'Participant 3'
};

% --- Layout Geometry ---
figWidth  = 1100; 
figHeight = 470; 

masterFig = figure('Position', [100, 100, figWidth, figHeight], ...
                   'Color', 'w', 'Renderer', 'opengl');

% Normalized margin and grid spacing
leftMargin   = 0.05;
rightMargin  = 0.25; % Extra room on the right for the floating legend
topMargin    = 0.08; 
bottomMargin = 0.06; 

colSpacing = 0.07; 
rowSpacing = 0.0; 

% Panel dimensions
panelWidth  = (1 - leftMargin - rightMargin - (numCols - 1) * colSpacing) / numCols;
panelHeight = (1 - topMargin - bottomMargin - (numRows - 1) * rowSpacing) / numRows;

% --- Main Assembly Loop ---
for c = 1:numCols

    % Open source figure invisibly
    if ~isfile(figFiles{c})
        warning('File %s not found. Skipping panel.', figFiles{c});
        continue;
    end
    
    srcFig = openfig(figFiles{c}, 'invisible');
    srcAx  = gca(srcFig);

    for r = 1:numRows
        
        % Calculate target panel position
        xPos = leftMargin + (c - 1) * (panelWidth + colSpacing);
        yPos = 1 - topMargin - r * panelHeight - (r - 1) * rowSpacing;        
        
        % Create target axes in the master figure and copy all 3D objects 
        % (heads, sensors, lights) from source to target
        if c == 1 && r == 1
            srcLegend = findobj(srcFig, 'Type', 'Legend');
            newLeg = copyobj([srcAx; srcLegend], masterFig);
            newLeg(1).Position = [xPos, yPos, panelWidth, panelHeight];
        else
            targetAx = axes(masterFig, 'Position', [xPos, yPos, panelWidth, panelHeight]);
            childrenObjs = allchild(srcAx);
            copyobj(childrenObjs, targetAx);
        end
        targetAx = gca(masterFig);
        
        % Transfer 3D camera, lighting, and view properties
        set(targetAx, ...
            'View',             get(srcAx, 'View'), ...
            'CameraPosition',   get(srcAx, 'CameraPosition'), ...
            'CameraTarget',     get(srcAx, 'CameraTarget'), ...
            'CameraUpVector',   get(srcAx, 'CameraUpVector'), ...
            'CameraViewAngle',  get(srcAx, 'CameraViewAngle'), ...
            'DataAspectRatio',  get(srcAx, 'DataAspectRatio'), ...
            'PlotBoxAspectRatio', get(srcAx, 'PlotBoxAspectRatio'), ...
            'XLim',             get(srcAx, 'XLim'), ...
            'YLim',             get(srcAx, 'YLim'), ...
            'ZLim',             get(srcAx, 'ZLim'));

        axis(targetAx, 'tight');
        if r == 1
            view(targetAx, [90,0]);
        elseif r == 2
            view(targetAx, [-90,0]);
        end        
        
        % Formatting & Tight Bounding Box
        axis(targetAx, 'off'); % Hide 3D spatial ticks
        
        % Create a clean outer 2D border around the panel canvas
        borderAx = axes(masterFig, 'Position', [xPos-0.015, yPos-0.015, panelWidth+0.03, panelHeight+0.03], ...
                        'Color', 'none', 'XColor', 'k', 'YColor', 'k', ...
                        'LineWidth', 1.2, 'XTick', [], 'YTick', []);
        box(borderAx, 'off');
        set(borderAx, 'Visible', 'off')
        
        % Add "Participant X" Header above the panel
        if r == 1
            t = title(targetAx, partTitles{c}, 'FontName', 'Arial', ...
                'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
            t.Units = 'normalized';
            t.Position(2) = 1.08;
        end
        
        % Add "A" (Anterior/Front) and "P" (Posterior/Back) Labels
        if r == 1
            text(borderAx, 0, 0.52, 'P', 'Units', 'normalized', ...
                'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
                'FontAngle', 'italic', 'HorizontalAlignment', 'center');
                
            text(borderAx, 0.97, 0.52, 'A', 'Units', 'normalized', ...
                'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
                'FontAngle', 'italic', 'HorizontalAlignment', 'center');
        elseif r == 2
            text(borderAx, 0, 0.52, 'A', 'Units', 'normalized', ...
                'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
                'FontAngle', 'italic', 'HorizontalAlignment', 'center');
                
            text(borderAx, 0.97, 0.52, 'P', 'Units', 'normalized', ...
                'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
                'FontAngle', 'italic', 'HorizontalAlignment', 'center');
        end
        
    end
    % Close source figure to free memory
    close(srcFig);
end

% Position Legend
set(newLeg(end), 'Position', [targetAx.Position(1) + panelWidth + 0.05, 0.5-0.15/2, 0.12, 0.15], ...
            'FontName', 'Arial', 'FontSize', 16);

% Save High-Resolution Output
exportgraphics(masterFig, fullfile(final_fileloc, '1_sensor_positions.png'), 'Resolution', 300);

%% Movement trajectories

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Config ---
open_closed = {'Closed', 'Open'};

% 2x2 Grid of MATLAB .fig files
figFiles = cell(2,2);

% Panel header labels
panelLabels = {
    '1)',   '2.1)';
    '2.2)', '3)'
};

% --- Figure Canvas Setup ---
numRows = 2;
numCols = 2;

figWidth  = 760;
figHeight = 800;

% Normalized margin and grid spacing
leftMargin   = 0.15;
rightMargin  = 0.05; % Extra room on the right for the floating legend
topMargin    = 0.08; 
bottomMargin = 0.1; 

colSpacing = 0.13; 
rowSpacing = 0.13; 

% Panel dimensions
plotW  = (1 - leftMargin - rightMargin - (numCols - 1) * colSpacing) / numCols;
plotH = (1 - topMargin - bottomMargin - (numRows - 1) * rowSpacing) / numRows;


for oc = 1:length(open_closed)

    figure('Position', [50, 50, figWidth, figHeight], 'Color', 'w');
    counter = 1;
    for ss = 1:length(sub)
        if strcmp(sub{ss}, 'sub-002')
            run = {'_run-001', '_run-002'};
        else
            run = {''};
        end
        for rr = 1:length(run)
            figFiles{counter} = fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_trajectory.fig', sub{ss}, open_closed{oc}, run{rr}));
            counter = counter + 1;
        end
    end
    figFiles = figFiles';


    for r = 1:2
        for c = 1:2
            % Calculate subpanel placement
            xPos = leftMargin + (c-1)*(plotW + colSpacing);
            yPos = bottomMargin + (2-r)*(plotH + rowSpacing);
            
            targetAx = axes('Position', [xPos, yPos, plotW, plotH]);
            
            figFile = figFiles{r, c};
            if isfile(figFile)
                srcFig = openfig(figFile, 'invisible');
                srcAx  = gca(srcFig);

                % To get axes in the right places, change to 2D plot from
                % above
                lines = allchild(srcAx);
                for ii = 1:length(lines)
                    lines(ii).YData = lines(ii).ZData;
                    lines(ii).ZData = [];
                end
                
                % Copy trajectory lines to target axes
                copyobj(lines, targetAx);

                % Standardize plot formatting
                % Transfer 3D camera, lighting, and view properties
                set(targetAx, ...
                    'DataAspectRatio',  get(srcAx, 'DataAspectRatio'), ...
                    'PlotBoxAspectRatio', get(srcAx, 'PlotBoxAspectRatio'), ...
                    'XLim',             get(srcAx, 'XLim'), ...
                    'YLim',             get(srcAx, 'ZLim'), ...
                    'XDir',             'reverse', ...
                    'FontSize',         16);

                box(targetAx, 'on'); 
                grid(targetAx, 'on');
                close(srcFig);
            end
            
            
            xlabel(targetAx, 'Left-Right (cm)', 'FontName', 'Arial', 'FontSize', 16);
            ylabel(targetAx, 'Forward-Back (cm)', 'FontName', 'Arial', 'FontSize', 16);
            
            % Add panel labels (1), 2.1), 2.2), 3)) outside top-left of each plot
            titleText = panelLabels{r, c};
            text(targetAx, -0.45, 1.02, titleText, 'Units', 'normalized', ...
                'FontName', 'Arial', 'FontSize', 18, 'FontWeight', 'bold');
        end
    end
    
    % Save High-Resolution Compiled Output
    exportgraphics(gcf, fullfile(final_fileloc, sprintf('movement_trajectories_composite_%s.png', open_closed{oc})), 'Resolution', 300);
end

%% Movement and Magnetic Field Change

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration & File Definitions ---
open_closed = {'Closed', 'Open'};

rowLabels = {'1)', '2.1)', '2.2)', '3)'};
colHeadersC = {'Left-Right', 'Forward-Back', 'Yaw'};

letterFontSize = 16;
labelFontSize = 15;
plotFontSize = 14;

% --- Figure Canvas Setup ---
figWidth  = 1400;
figHeight = 800;

% --- Geometry Layout Margins (Normalized 0 to 1) ---
% Left Column (Sections A & B)
leftColX      = 0.06;
leftColWidth  = 0.36;

% Right Column (Section C)
rightColX     = 0.51;
rightColWidth = 0.49;

% Vertical Bounds
topY    = 0.92;
bottomY = 0.08;

for oc = 1:length(open_closed)

    if strcmp(open_closed{oc}, 'Closed')
        file_SectionA = fullfile(results_dir, 'Magnetic_field_histogram.fig');
        file_SectionB = fullfile(results_dir, 'Magnetic_field_time_series.fig'); % File containing the tiledLayout for (B)
    else
        file_SectionA = fullfile(results_dir, 'Magnetic_field_histogram_open_loop.fig');
        file_SectionB = fullfile(results_dir, 'Magnetic_field_time_series_open_loop.fig'); % File containing the tiledLayout for (B)
    end


    % 4x3 Grid of .fig files for Section C (Rows 1-4, Cols 1-3)
    files_SectionC = cell(4,3);
    counter = 1;
    for ss = 1:length(sub)
        if strcmp(sub{ss}, 'sub-002')
            run = {'_run-001', '_run-002'};
        else
            run = {''};
        end
        for rr = 1:length(run)
            files_SectionC(counter,:) = {
                fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_displacement_Left-Right.fig', sub{ss}, open_closed{oc}, run{rr})), ...
                fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_displacement_Forward-Back.fig', sub{ss}, open_closed{oc}, run{rr})), ...
                fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_rotation_Yaw.fig', sub{ss}, open_closed{oc}, run{rr}))
            };
            counter = counter + 1;
        end
    end

    masterFig = figure('Position', [50, 50, figWidth, figHeight], 'Color', 'w');


    % =========================================================================
    % SECTION A: Top Left Plot
    % =========================================================================
    secA_pos = [leftColX, 0.58, leftColWidth, 0.37];
    
    if isfile(file_SectionA)
        srcFigA = openfig(file_SectionA, 'invisible');
        srcAxA  = gca(srcFigA);
        legA = findobj(srcFigA, 'Type', 'Legend');
        
        % Copy content and legend
        axA = copyobj([srcAxA, legA], masterFig);
        set(axA(1), 'Position', secA_pos, 'XLim', get(srcAxA, 'XLim'), 'YLim', get(srcAxA, 'YLim'), ...
                 'Box', 'on', 'FontName', 'Arial', 'FontSize', plotFontSize);
        close(srcFigA);
    else
        % Placeholder frame
        box(axA, 'on');
        text(axA, 0.5, 0.5, 'Section A (.fig)', 'HorizontalAlignment', 'center');
    end
    
    % Label "A)"
    annotation('textbox', [0, 0.95, 0.05, 0.04], ...
        'String', '(A)', 'FontName', 'Arial', 'FontSize', letterFontSize, 'FontWeight', 'bold', 'EdgeColor', 'none');
    
    
    % =========================================================================
    % SECTION B: Bottom Left Tiled Layout
    % =========================================================================
    % Create a uipanel to host the imported tiledLayout perfectly
    secB_pos = [leftColX, 0.07, leftColWidth, 0.4];
    % panelB = uipanel('Parent', masterFig, 'Position', secB_pos, ...
    %                  'BorderType', 'none', 'BackgroundColor', 'w');
    
    if isfile(file_SectionB)
        srcFigB = openfig(file_SectionB);
        tiledObj = findobj(srcFigB, 'Type', 'TiledLayout');
        
        if ~isempty(tiledObj)
            % Copy the tiledLayout directly into panelB
            panelB = copyobj(tiledObj, masterFig);
            panelB.InnerPosition = secB_pos;
            panelBChildren = panelB.Children;
            for ii = 1:length(panelBChildren)
                set(panelBChildren(ii), 'FontSize', plotFontSize);
            end
        else
            warning('No TiledLayout found in %s', file_SectionB);
        end
        close(srcFigB);
    else
        % Placeholder frame
        axes('Parent', panelB, 'Position', [0 0 1 1]); box on;
        text(0.5, 0.5, 'Section B TiledLayout (.fig)', 'HorizontalAlignment', 'center');
    end
    
    % Label "(B)"
    annotation(masterFig,'textbox', [0, 0.44, 0.05, 0.04], ...
        'String', '(B)', 'FontName', 'Arial', 'FontSize', letterFontSize, 'FontWeight', 'bold', 'EdgeColor', 'none');
    
    
    % =========================================================================
    % DIVIDER LINES (Left Column & Center)
    % =========================================================================
    % Horizontal divider line between Section A and Section B
    annotation('line', [0, leftColX + leftColWidth + 0.02], [0.50, 0.50], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);

    % Vertical divider line between Section A/B and Section C
    annotation('line', [leftColX + leftColWidth + 0.02, leftColX + leftColWidth + 0.02], [0, 1], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    
    
    % =========================================================================
    % SECTION C: Right Side 4x3 Grid
    % =========================================================================
    % Label "(C)"
    annotation(masterFig, 'textbox', [leftColX + leftColWidth + 0.03, 0.95, 0.05, 0.04], ...
        'String', '(C)', 'FontName', 'Arial', 'FontSize', letterFontSize, 'FontWeight', 'bold', 'EdgeColor', 'none');
    
    numRowsC = 4;
    numColsC = 3;
    rowSpacingC = 0.03;
    
    colWidthC  = (rightColWidth - 0.06) / numColsC;
    rowHeightC = (0.82 - rowSpacingC*(numRowsC-1)) / numRowsC;
    
    for r = 1:numRowsC
        % Y-Position for current row
        yPos = topY - r * rowHeightC - (r-1)*rowSpacingC;
        
        % Row Labels: "1)", "2.1)", "2.2)", "3)"
        annotation(masterFig,'textbox', [leftColX + leftColWidth + 0.034, yPos + rowHeightC - 0.01, 0.03, 0.04], ...
            'String', rowLabels{r}, 'FontName', 'Arial', 'FontSize', labelFontSize, ...
            'FontWeight', 'bold', 'EdgeColor', 'none', 'HorizontalAlignment', 'left');
            
        for c = 1:numColsC
            xPos = rightColX + (c - 1) * (colWidthC + 0.02);
            
            % targetAx = axes('Position', [xPos, yPos + 0.01, colWidthC, rowHeightC - 0.02]);
            figFile = files_SectionC{r, c};
            
            if isfile(figFile)
                srcFigC = openfig(figFile, 'invisible');
                srcAxC  = gca(srcFigC);
                
                % Copy plot elements
                targetAx = copyobj(srcAxC, masterFig);
                if c~= numColsC
                    targetAx.Position = [xPos, yPos + 0.01, colWidthC, rowHeightC];
                else
                    targetAx.Position = [xPos, yPos + 0.01, colWidthC, rowHeightC-0.04];
                end
                
                % Mirror properties and limits
                if c < 3
                    set(targetAx, 'XLim', get(srcAxC, 'XLim'), 'YLim', get(srcAxC, 'YLim'), ...
                         'Box', 'on', 'FontName', 'Arial', 'FontSize', plotFontSize);
                else
                    set(targetAx, 'FontSize', plotFontSize);
                end
                
                close(srcFigC);
            else
                % Fallback plot placeholders
                box(targetAx, 'on');
                set(targetAx, 'XTick', [], 'YTick', []);
            end
            
            % Column Headings on top row
            if r == 1
                title(targetAx, colHeadersC{c}, 'FontName', 'Arial', 'FontSize', plotFontSize, 'FontWeight', 'bold');
            else
                title(targetAx, '');
            end
            
            % X-Axis labels on bottom row
            if r < 4 && c < 3
                xlabel(targetAx, '');
                set(targetAx, 'XTickLabel', []);
            end
            if c == 2
                ylabel(targetAx, '');
                set(targetAx, 'YTickLabel', []);
            end
            
        end
        
        
    end
    
    % Save High-Resolution Compiled Output
    exportgraphics(masterFig, fullfile(final_fileloc, sprintf('%s_loop_movement_field_change_composite.png', open_closed{oc})), 'Resolution', 300);
end

%% PSD figures

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration & File Definitions ---
open_closed = {'Closed', 'Open'};

for oc = 1:length(open_closed)
    % 1. Configuration
    % Specify the paths/filenames for your four individual .fig files.
    % Replace these placeholders with your actual filenames.
    filename_top_left = fullfile(results_dir, sprintf('PSD_%s_seated.fig', lower(open_closed{oc})));
    filename_top_right = fullfile(results_dir, sprintf('ShieldingFactor_%s_seated.fig', lower(open_closed{oc})));
    filename_bottom_left = fullfile(results_dir, sprintf('PSD_%s_walking.fig', lower(open_closed{oc})));
    filename_bottom_right = fullfile(results_dir, sprintf('ShieldingFactor_%s_walking.fig', lower(open_closed{oc})));

    % 2. Overall Figure Setup
    % Initialize a new, empty figure with white background and standard landscape size.
    figure_handle = figure('Name', 'Combined PSD and Interference Reduction plots', ...
        'Color', 'w', ...
        'Position', [100, 100, 1100, 800]); % Define size: [left bottom width height]
    
    % Use tiledlayout for clean scientific journal spacing (introduced in R2019b).
    % 'TileSpacing', 'tight' and 'Padding', 'tight' minimize empty space.
    layout = tiledlayout(2, 2, 'TileSpacing', 'loose');
    
    % 3. Copy Plots to New Figure and Reformat
    % Prepare a list of destination indices and their source filenames.
    plot_config = {
        1, filename_top_left,  'i)'; ... % Top-left
        2, filename_top_right, 'ii)'; ...% Top-right
        3, filename_bottom_left,  'i)'; ... % Bottom-left
        4, filename_bottom_right, 'ii)';... % Bottom-right
        };
    
    for k = 1:size(plot_config, 1)
        % Parse configuration for the current tile
        tile_index = plot_config{k, 1};
        source_filename = plot_config{k, 2};
        panel_label = plot_config{k, 3};
        
        % Check if source file exists
        if exist(source_filename, 'file') ~= 2
            warning('Source file not found: %s. Skipping tile %d.', source_filename, tile_index);
            continue;
        end
        
        % a) Select the active tile. This creates a placeholder axes object.
        % nexttile returns the handle to this new placeholder axis.
        next_ax = nexttile(tile_index);
        
        % b) Open source figure invisibly to access its contents.
        source_fig_handle = openfig(source_filename, 'invisible');
        source_ax_handle = gca(source_fig_handle); % Get original axes handle
    
        % c) IMPORTANT: Check for an associated legend object in the source figure.
        % copyobj will handle copying the visual aspects, but re-parenting legends
        % in tiledlayouts requires explicit handling.
        source_legend_handle = findobj(source_fig_handle, 'Type', 'Legend');
    
        % d) Perform the whole-object copy.
        % copyobj(obj, new_parent) makes a complete duplicate of source_ax_handle,
        % including all internal properties, settings, colors, children, and association
        % with existing legends, and makes figure_handle its immediate parent.
        new_ax = copyobj([source_ax_handle, source_legend_handle], figure_handle);

        if k == 1 || k == 3
            ylabel(new_ax(1), '$$PSD (fT\sqrt[-1]{Hz})$$','interpreter','latex');
        end
    
        % e) Parent the new axes duplicate to the tiledlayout.
        % This moves the duplicate axis into the grid structure.
        set(new_ax(1), 'Parent', layout);
        
        % Set its position within the layout.
        new_ax(1).Layout.Tile = tile_index;
    
        % f) Finalize visual transfer.
        % next_ax was only a placeholder needed to advance nexttile and set its bounds.
        % Since new_ax (the complete copy) is now correctly positioned, we can delete the placeholder.
        delete(next_ax);
    
        % g) Add panel label (i) or ii).
        % Use text anchored normalized relative to the new axes object, matching 
        % scientific figure panel label style. Ensure font matches original figure style.
        title_label = text(new_ax(1), -0.2, 1.03, panel_label, ...
            'Units', 'normalized', ...
            'FontName', get(new_ax(1), 'FontName'), 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

        set(new_ax(1), 'FontSize', 14);
        
        % h) Close the invisible source figure instance to release memory.
        close(source_fig_handle);
    end
    
    % Find Tile 1 (Row 1) and Tile 3 (Row 2) axes to anchor annotations
    all_axes = findobj(figure_handle, 'Type', 'axes');
    tile1_ax = [];
    tile3_ax = [];
    
    for ax = all_axes'
        if ~isempty(ax.Layout) && isprop(ax.Layout, 'Tile')
            if ax.Layout.Tile == 1, tile1_ax = ax; end
            if ax.Layout.Tile == 3, tile3_ax = ax; end
        end
    end
    
    if ~isempty(tile1_ax) && ~isempty(tile3_ax)
        pos1 = tile1_ax.Position; % [left, bottom, width, height]
        pos3 = tile3_ax.Position;
        
        % --- Add Row Heading "a) Seated" above Row 1 ---
        annotation('textbox', [0.03, pos1(2) + pos1(4) + 0.035, 0.30, 0.04], ...
            'String', 'a) Seated', ...
            'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'FitBoxToText', 'on');
    
        % --- Add Row Heading "b) Walking" above Row 2 ---
        annotation('textbox', [0.03, pos3(2) + pos3(4) + 0.035, 0.30, 0.04], ...
            'String', 'b) Walking', ...
            'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'FitBoxToText', 'on');
    
        % --- Add Horizontal Divider Line between Row 1 and Row 2 ---
        % Midpoint between the bottom of Row 1 and the top of Row 2
        y_line = (pos1(2) + (pos3(2) + pos3(4))) / 2; 
        
        annotation('line', [0.03, 0.96], [y_line, y_line], ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    else
        % Fallback coordinates if tile structure query fails
        annotation('textbox', [0.06, 0.94, 0.30, 0.04], 'String', 'a) Seated', ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', 'EdgeColor', 'none');
        annotation('textbox', [0.06, 0.46, 0.30, 0.04], 'String', 'b) Walking', ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', 'EdgeColor', 'none');
        annotation('line', [0.04, 0.96], [0.50, 0.50], 'Color', 'k', 'LineWidth', 1.5);
    end

    % 4. Save
    exportgraphics(figure_handle, fullfile(final_fileloc, sprintf('%s_psd_composite.png', open_closed{oc})), 'Resolution', 300);
end

%% PSD figures with HFC with linear gradient terms

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration & File Definitions ---
open_closed = {'Closed', 'Open'};

for oc = 1:length(open_closed)
    % 1. Configuration
    % Specify the paths/filenames for your four individual .fig files.
    % Replace these placeholders with your actual filenames.
    filename_top_left = fullfile(results_dir, sprintf('PSD_%s_seated_hfcgrad.fig', lower(open_closed{oc})));
    filename_top_right = fullfile(results_dir, sprintf('ShieldingFactor_%s_seated_hfcgrad.fig', lower(open_closed{oc})));
    filename_bottom_left = fullfile(results_dir, sprintf('PSD_%s_walking_hfcgrad.fig', lower(open_closed{oc})));
    filename_bottom_right = fullfile(results_dir, sprintf('ShieldingFactor_%s_walking_hfcgrad.fig', lower(open_closed{oc})));

    % 2. Overall Figure Setup
    % Initialize a new, empty figure with white background and standard landscape size.
    figure_handle = figure('Name', 'Combined PSD and Interference Reduction plots', ...
        'Color', 'w', ...
        'Position', [100, 100, 1100, 800]); % Define size: [left bottom width height]
    
    % Use tiledlayout for clean scientific journal spacing (introduced in R2019b).
    % 'TileSpacing', 'tight' and 'Padding', 'tight' minimize empty space.
    layout = tiledlayout(2, 2, 'TileSpacing', 'loose');
    
    % 3. Copy Plots to New Figure and Reformat
    % Prepare a list of destination indices and their source filenames.
    plot_config = {
        1, filename_top_left,  'i)'; ... % Top-left
        2, filename_top_right, 'ii)'; ...% Top-right
        3, filename_bottom_left,  'i)'; ... % Bottom-left
        4, filename_bottom_right, 'ii)';... % Bottom-right
        };
    
    for k = 1:size(plot_config, 1)
        % Parse configuration for the current tile
        tile_index = plot_config{k, 1};
        source_filename = plot_config{k, 2};
        panel_label = plot_config{k, 3};
        
        % Check if source file exists
        if exist(source_filename, 'file') ~= 2
            warning('Source file not found: %s. Skipping tile %d.', source_filename, tile_index);
            continue;
        end
        
        % a) Select the active tile. This creates a placeholder axes object.
        % nexttile returns the handle to this new placeholder axis.
        next_ax = nexttile(tile_index);
        
        % b) Open source figure invisibly to access its contents.
        source_fig_handle = openfig(source_filename, 'invisible');
        source_ax_handle = gca(source_fig_handle); % Get original axes handle
    
        % c) IMPORTANT: Check for an associated legend object in the source figure.
        % copyobj will handle copying the visual aspects, but re-parenting legends
        % in tiledlayouts requires explicit handling.
        source_legend_handle = findobj(source_fig_handle, 'Type', 'Legend');
    
        % d) Perform the whole-object copy.
        % copyobj(obj, new_parent) makes a complete duplicate of source_ax_handle,
        % including all internal properties, settings, colors, children, and association
        % with existing legends, and makes figure_handle its immediate parent.
        new_ax = copyobj([source_ax_handle, source_legend_handle], figure_handle);

        if k == 1 || k == 3
            ylabel(new_ax(1), '$$PSD (fT\sqrt[-1]{Hz})$$','interpreter','latex');
        end
    
        % e) Parent the new axes duplicate to the tiledlayout.
        % This moves the duplicate axis into the grid structure.
        set(new_ax(1), 'Parent', layout);
        
        % Set its position within the layout.
        new_ax(1).Layout.Tile = tile_index;
    
        % f) Finalize visual transfer.
        % next_ax was only a placeholder needed to advance nexttile and set its bounds.
        % Since new_ax (the complete copy) is now correctly positioned, we can delete the placeholder.
        delete(next_ax);
    
        % g) Add panel label (i) or ii).
        % Use text anchored normalized relative to the new axes object, matching 
        % scientific figure panel label style. Ensure font matches original figure style.
        title_label = text(new_ax(1), -0.2, 1.03, panel_label, ...
            'Units', 'normalized', ...
            'FontName', get(new_ax(1), 'FontName'), 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

        set(new_ax(1), 'FontSize', 14);
        
        % h) Close the invisible source figure instance to release memory.
        close(source_fig_handle);
    end
    
    % Find Tile 1 (Row 1) and Tile 3 (Row 2) axes to anchor annotations
    all_axes = findobj(figure_handle, 'Type', 'axes');
    tile1_ax = [];
    tile3_ax = [];
    
    for ax = all_axes'
        if ~isempty(ax.Layout) && isprop(ax.Layout, 'Tile')
            if ax.Layout.Tile == 1, tile1_ax = ax; end
            if ax.Layout.Tile == 3, tile3_ax = ax; end
        end
    end
    
    if ~isempty(tile1_ax) && ~isempty(tile3_ax)
        pos1 = tile1_ax.Position; % [left, bottom, width, height]
        pos3 = tile3_ax.Position;
        
        % --- Add Row Heading "a) Seated" above Row 1 ---
        annotation('textbox', [0.03, pos1(2) + pos1(4) + 0.035, 0.30, 0.04], ...
            'String', 'a) Seated', ...
            'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'FitBoxToText', 'on');
    
        % --- Add Row Heading "b) Walking" above Row 2 ---
        annotation('textbox', [0.03, pos3(2) + pos3(4) + 0.035, 0.30, 0.04], ...
            'String', 'b) Walking', ...
            'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', ...
            'EdgeColor', 'none', 'FitBoxToText', 'on');
    
        % --- Add Horizontal Divider Line between Row 1 and Row 2 ---
        % Midpoint between the bottom of Row 1 and the top of Row 2
        y_line = (pos1(2) + (pos3(2) + pos3(4))) / 2; 
        
        annotation('line', [0.03, 0.96], [y_line, y_line], ...
            'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
    else
        % Fallback coordinates if tile structure query fails
        annotation('textbox', [0.06, 0.94, 0.30, 0.04], 'String', 'a) Seated', ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', 'EdgeColor', 'none');
        annotation('textbox', [0.06, 0.46, 0.30, 0.04], 'String', 'b) Walking', ...
            'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', 'EdgeColor', 'none');
        annotation('line', [0.04, 0.96], [0.50, 0.50], 'Color', 'k', 'LineWidth', 1.5);
    end

    % 4. Save
    exportgraphics(figure_handle, fullfile(final_fileloc, sprintf('%s_psd_composite_hfcgrad.png', open_closed{oc})), 'Resolution', 300);
end

%% HFC with gradients evoked responses

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration ---
numRows = 4;
numCols = 6;
task = {'seatedClosed', 'seatedOpen', 'walkingClosed', 'walkingOpen'};
row_folder_name = {'hfc_with_gradients'};

% --- Layout Geometry ---
figWidth = 1534; % Figure width in pixels
figHeight = 700; % Figure height in pixels

% Margins and spacing (normalized 0 to 1)
leftMargin   = 0.18;
rightMargin  = 0.07;
topMargin    = 0.05;
bottomMargin = 0.15;

colSpacing = 0.04; % Space between columns
rowSpacing = 0.04; % Space between rows

% Calculate width and height of individual panels (normalised)
panelWidth  = (1 - leftMargin - rightMargin - (numCols - 1) * colSpacing) / numCols;
panelHeight = (1 - topMargin - bottomMargin - (numRows - 1) * rowSpacing) / numRows;

% --- Pre-calculate Y-Coordinates for Top-Alignment ---
% To keep the grid perfectly aligned, we establish a global top Y-coordinate for each row.
rowYTop = zeros(numRows, 1);
rowYTop(1) = 1 - topMargin;

for r = 2:numRows
    rowYTop(r) = 1 - topMargin - (r-1)*(panelHeight + rowSpacing);
end

% Define the filenames of your 24 sub-panels (Row x Column)
% Adjust these paths/names to match your actual files
imageFiles = cell(numRows, numCols);


% Loop over recordings
for tt = 1:length(task)

    f1 = figure('Position', [0, 50, figWidth, figHeight], 'Color', 'w');
    ax = gobjects(numRows, numCols);

    % Loop over participants
    counter = 0;
    for ss = 1:length(sub)

        % Loop over runs
        if strcmp(sub{ss}, 'sub-002') && startsWith(task{tt}, 'walking')
            run = {'_run-001', '_run-002'};
        else
            run = {''};
        end

        for rr = 1:length(run)
            counter = counter + 1;
            imageFiles{counter, 1} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_average_time_series.fig', sub{ss}, task{tt}, run{rr})); 
            imageFiles{counter, 2} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_t_stat_time_series.fig', sub{ss}, task{tt}, run{rr})); 
            imageFiles{counter, 3} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_anti_averaging_topo_80_to_120_ms_anysig.fig', sub{ss}, task{tt}, run{rr})); % Ultimately add _test to end of filename
            imageFiles{counter, 4} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_ROI.fig', sub{ss}, task{tt}, run{rr}));
            imageFiles{counter, 5} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_min_norm_pow_left.fig', sub{ss}, task{tt}, run{rr}));
            imageFiles{counter, 6} = fullfile(results_dir, sub{ss}, row_folder_name{1}, sprintf('%s_task-%s%s_meg_min_norm_pow_right.fig', sub{ss}, task{tt}, run{rr}));
        end
    end

            
    % --- Main Plotting Loop ---
    for r = 1:counter
        for c = 1:numCols
                
            % 2. Load and Crop Image
            if isfile(imageFiles{r, c})
                f = openfig(imageFiles{r,c}, 'invisible');
                a = f.Children;

                % Calculate Position
                xPos = leftMargin + (c - 1) * (panelWidth + colSpacing);
                if c == 3
                    xPos = xPos-0.015;
                end
                % Align the tops of the plots; the extra bottom height will extend downward
                yPos = rowYTop(r) - panelHeight; 
                
                if c == 3 || c == 4 || c == 5 || c == 6
                    ax(r,c) = copyobj(a(end), f1);
                    if c == 3
                        colormap(ax(r,c), colormap123);
                        ax(r,c).CameraViewAngle = 4.67;
                        if r == counter
                            cb = colorbar(ax(r,c), "southoutside");
                            ylabel(cb, 'B (fT)')
                        end
                    else
                        colormap(ax(r,c), 'hot');
                    end
                    if c == 5 || c == 6
                        if length(a) == 3
                            ax(r,c).CameraViewAngle = 4.0797;
                        else
                            if ss == 1
                                ax(r,c).CameraViewAngle = 6.89;
                            else
                                if r == 3
                                    ax(r,c).CameraViewAngle = 4.12;
                                else
                                    ax(r,c).CameraViewAngle = 6.40;
                                end
                            end
                        end
                        cb3_ax = colorbar(ax(r,c));
                    end
                else
                    ax(r, c) = copyobj(a, f1);
                end
                ax(r,c).Position = [xPos, yPos, panelWidth, panelHeight];
                set(ax(r,c), 'FontSize', 16);

                % Make colorbar smaller
                if c == 5 || c == 6
                    cb3_ax.Position(2) = cb3_ax.Position(2) + cb3_ax.Position(4)/2 - cb3_ax.Position(4)*0.6/2;
                    cb3_ax.Position(4) = cb3_ax.Position(4)*0.6;
                end

                if r == counter && c == 3
                    pos = cb.Position;
                    cb.Position(3) = 0.8*panelWidth;
                    cb.Position(4) = 1.2*cb.Position(4);
                    cb.Position(1) = pos(1) + pos(3)/2 - cb.Position(3)/2;
                    cb.Position(2) = ax(r,c).Position(2) - cb.Position(4) - 0.02;
                    cb.FontSize = 16;
                end

                axis tight
                
                ylabel(ax(r,c), '');
                if r < counter
                    set(ax(r,c), 'XTickLabel', []);
                    xlabel(ax(r,c), '');
                end


            else
                % If file is missing
                continue;
            end
        end
    end
        
    % --- Add Overall Labels (using normalized annotations) ---

    figure(f1);
    
    % Y-Axis Labels
    % B (fT) centered vertically on Column 1
    yCenter_col1 = (ax(2,1).Position(2) + ax(2,1).Position(4) + ax(counter-1,1).Position(2))/2;
    annotation('textarrow', [1 1]*leftMargin-colSpacing*1.2, [yCenter_col1 yCenter_col1], ...
        'String', 'B (fT)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
        'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
        'VerticalAlignment','middle');
    
    % t-stat centered vertically on Column 2
    annotation('textarrow', [1 1]*leftMargin+panelWidth+0.009, [yCenter_col1 yCenter_col1], ...
        'String', 't-stat', 'HeadStyle', 'none', 'LineStyle', 'none', ...
        'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
        'VerticalAlignment', 'middle');
    
    % Current (nAm) centered vertically on Column 4
    xPos_col4 = leftMargin + 3 * (panelWidth + colSpacing) - 0.035;
    annotation('textarrow', [xPos_col4 xPos_col4], [yCenter_col1 yCenter_col1], ...
        'String', 'Current (nAm)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
        'TextRotation', 90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18, ...
        'VerticalAlignment', 'middle');

    % % Add colourbars
    % cb2_ax = colorbar(ax(1,5));
    % cb2_ax.Position = [ax(1,5).Position(1) + ax(1,5).Position(3)+0.005, ...
    %     yCenter_col1-panelHeight*2/2, 0.01, panelHeight*2];
    % 
    % cb3_ax = colorbar(ax(1,6));
    % cb3_ax.Position = [ax(1,6).Position(1) + ax(1,6).Position(3)+0.005, ...
    %     yCenter_col1-panelHeight*2/2, 0.01, panelHeight*2];

    annotation('textarrow', (cb3_ax.Position(1) + cb3_ax.Position(3) + 0.025)*[1, 1], ...
        yCenter_col1*[1 1], ...
        'String', 'Minimum Norm Power (a.u.)', 'HeadStyle', 'none', 'LineStyle', 'none', ...
        'TextRotation', -90, 'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 15, ...
        'VerticalAlignment', 'middle');
    
    % Add legend
    lgd = legend(ax(1, 4), {'Left', 'Right'});
    lgd.Position = [ax(1,4).Position(1)+ax(1,4).Position(3)-panelWidth*0.45, ...
            ax(1,4).Position(2)+ax(1,4).Position(4)-0.04, lgd.Position(3), lgd.Position(4)];

    % Add Column Numbering Headers at the Top
    columnHeaders = {'A) i)', 'ii)', 'iii)', 'B) i)', 'ii)', 'iii)'};
    for c = 1:numCols
        if c == 1 || c== 2 || c == 4
            xCenter = ax(1,c).Position(1) - 0.01;
        else
            xCenter = ax(1,c).Position(1) - 0.002;
        end
        % Position the text 1.5% of figure height above the top of the first row of panels
        yPos_header = rowYTop(1) + 0.025; 
        
        annotation('textarrow', [xCenter xCenter], [yPos_header yPos_header], ...
            'String', columnHeaders{c}, 'HeadStyle', 'none', 'LineStyle', 'none', ...
            'HorizontalAlignment', 'center', 'FontName', 'Arial', 'FontSize', 18);
    end

    % Add Wrapped Row Names on the Left
    % Nested cells control multi-line wrapping perfectly
    if startsWith(task{tt}, 'walking')
        rowNames = {
            'Participant 1', ...
            'Participant 2.1)', ...
            'Participant 2.2)', ...
            'Participant 3' ...
        };
    else
        rowNames = {
            'Participant 1', ...
            'Participant 2', ...
            'Participant 3' ...
        };
    end
    
    for r = 1:length(rowNames)
        boxWidth = leftMargin*0.65; 
        boxHeight = 0.12;
        boxX = 0;                  % Placed on the far-left edge
        rowYCenter = ax(r, 1).Position(2) + ax(r,1).Position(4)/2;
        boxY = rowYCenter - (boxHeight / 2);
        
        annotation('textbox', [boxX, boxY, boxWidth, boxHeight], ...
            'String', rowNames{r}, ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', ...
            'FontSize', 16, ...
            'FontWeight', 'bold', ...
            'FontAngle', 'italic');
    end

    % Save output
    exportgraphics(gcf, fullfile(final_fileloc, sprintf('hfcgrad_combined_%s.png', task{tt})), 'Resolution', 300);
    close all;
    
end


%% Movement speed figures

clearvars -except colormap123 results_dir sub final_fileloc; clc; close all;

% --- Configuration & File Definitions ---
open_closed = {'Closed', 'Open'};

letterFontSize = 16;
labelFontSize = 15;
plotFontSize = 14;


% 1. Configuration: Figure Filenames and Structure
figWidth  = 1600;
figHeight = 550;

% Section A (Left Plot): ~32% total width
secA_x      = 0.06;
secA_y      = 0.12;
secA_width  = 0.28;
secA_height = 0.76;

% Vertical Divider Line Location
divider_x   = secA_width+secA_x+0.01;

% Section B (Right 2x4 Grid): ~55% total width
secB_x          = secA_width+secA_x+0.07;
secB_totalWidth = 0.55;
colSpacing      = 0.02; % Horizontal space between subplots
rowSpacing      = 0.15;  % Vertical space between rows

colWidth  = (secB_totalWidth - 3 * colSpacing) / 4; % Width per subplot (~0.128)
rowHeight = (secA_height - rowSpacing) / 2;         % Height per subplot (~0.34)

for oc = 1:length(open_closed)
    % Filename for Section A (single left plot)
    if strcmp(open_closed{oc}, 'Open')
        file_SectionA = fullfile(results_dir, 'Magnetic_field_change_rate_histogram_open_loop.fig');
    else
        file_SectionA = fullfile(results_dir, 'Magnetic_field_change_rate_histogram.fig');
    end
    
    % Cell array of filenames for Section B (2x4 grid of right plots)
    % Row 1 is 'Left-Right' plots, Row 2 is 'Forward-Back' plots.
    files_SectionB = cell(2, 4);
    
    counter = 0;
    for ss = 1:length(sub)
        if strcmp(sub{ss}, 'sub-002')
            run = {'_run-001', '_run-002'};
        else
            run = {''};
        end
        for rr = 1:length(run)
            counter = counter + 1;
            files_SectionB{1, counter} = fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_speed_Left-Right.fig', sub{ss}, open_closed{oc}, run{rr}));
            files_SectionB{2, counter} = fullfile(results_dir, sub{ss}, sprintf('%s_task-walking%s%s_meg_speed_Forward-Back.fig', sub{ss}, open_closed{oc}, run{rr}));
        end
    end
    
    % Labels for the columns in Section B (corresponding to participant number)
    columnLabelsB = {'1)', '2.1)', '2.2)', '3)'};
    
    % 3. Create Main Combined Figure
    
    % Define a figure with a wide aspect ratio to accommodate the layout
    masterFig = figure('Position', [50, 100, figWidth, figHeight], 'Color', 'w');
    
    if isfile(file_SectionA)
        srcFigA = openfig(file_SectionA, 'invisible');
        srcAxA  = gca(srcFigA);
        srcLegA = findobj(srcFigA, 'Type', 'Legend');
        
        % Copy whole axes directly into master figure
        axA = copyobj([srcAxA, srcLegA], masterFig);
        set(axA(1), 'Position', [secA_x, secA_y, secA_width, secA_height], ...
                 'FontName', 'Arial', 'FontSize', plotFontSize, 'Box', 'on');
             
        close(srcFigA);
    else
        % Fallback frame if file is missing
        axA = axes('Position', [secA_x, secA_y, secA_width, secA_height]);
        box(axA, 'on');
        text(axA, 0.5, 0.5, 'plot\_A.fig', 'HorizontalAlignment', 'center');
    end
    
    % Label "A)"
    annotation('textbox', [secA_x - 0.05, 0.90, 0.05, 0.05], ...
        'String', '(A)', 'FontName', 'Arial', 'FontSize', letterFontSize, ...
        'FontWeight', 'bold', 'EdgeColor', 'none');
    
    
    % =========================================================================
    % DIVIDER LINE (Between Section A and B)
    % =========================================================================
    annotation('line', [divider_x, divider_x], [0.08, 0.94], ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    
    
    % =========================================================================
    % SECTION B: Right 2x4 Grid
    % =========================================================================
    % Label "(B)"
    annotation('textbox', [divider_x + 0.01, 0.90, 0.05, 0.05], ...
        'String', '(B)', 'FontName', 'Arial', 'FontSize', letterFontSize, ...
        'FontWeight', 'bold', 'EdgeColor', 'none');
    
    for r = 1:2
        for c = 1:4
            % Calculate explicit target axes position
            xPos = secB_x + (c - 1) * (colWidth + colSpacing);
            yPos = secA_y + (2 - r) * (rowHeight + rowSpacing);
            
            figFile = files_SectionB{r, c};
            
            if isfile(figFile)
                srcFigB = openfig(figFile, 'invisible');
                srcAxB  = gca(srcFigB);
                
                % Copy axes directly and reposition
                axB = copyobj(srcAxB, masterFig);
                set(axB, 'Position', [xPos, yPos, colWidth, rowHeight], ...
                         'FontName', 'Arial', 'FontSize', plotFontSize, 'Box', 'on');

                if c > 1
                    ylabel(axB, '');
                    set(axB, 'YTickLabel', []);
                end
                     
                close(srcFigB);
            else
                % Fallback frame
                axB = axes('Position', [xPos, yPos, colWidth, rowHeight]);
                box(axB, 'on');
            end
            
            % Add Column Participant Headings above Row 1: 1), 2.1), 2.2), 3)
            if r == 1
                annotation('textbox', [xPos-0.015, yPos + rowHeight + 0.005, colWidth, 0.04], ...
                    'String', columnLabelsB{c}, 'FontName', 'Arial', 'FontSize', labelFontSize, ...
                    'FontWeight', 'bold', 'EdgeColor', 'none', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            end
        end
    end
    
    % Finalize
    exportgraphics(gcf, fullfile(final_fileloc, sprintf('speed_combined_%s.png', open_closed{oc})), 'Resolution', 300);
end
