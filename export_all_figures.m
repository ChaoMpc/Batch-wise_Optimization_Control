function Report = export_all_figures(Figures,out_dir)
% =========================================================================
% export_all_figures
%
% Export the explicitly supplied figures to PNG, PDF, and EPS.
% Record each export result independently.
% =========================================================================

    if exist(out_dir,'dir') ~= 7
        [ok,message] = mkdir(out_dir);

        if ~ok
            error('export_all_figures: cannot create directory: %s',message);
        end
    end

    Report = struct( ...
        'figure_name',{},'format',{},'path',{}, ...
        'success',{},'message',{});

    formats = {'png','pdf','eps'};
    devices = {'-dpng','-dpdf','-depsc2'};

    index = 0;

    % 获取当前屏幕 DPI (用于将像素转换为英寸)
    try
        screen_dpi = get(0, 'ScreenPixelsPerInch');
    catch
        screen_dpi = 96;
    end
    if isempty(screen_dpi) || ~isscalar(screen_dpi) || screen_dpi <= 0
        screen_dpi = 96;
    end

    for j = 1:numel(Figures)
        fh = Figures(j).handle;
        name = Figures(j).name;

        if ~ishandle(fh)
            error('export_all_figures: invalid handle for %s.',name);
        end

        name = regexprep(name,'[^a-zA-Z0-9_-]','_');

        % --- 修复图片被裁剪的问题 ---
        % 动态计算图窗的物理尺寸，并强制设置 PaperSize 以包含整个图表
        pos = get(fh, 'Position');
        W_in = pos(3) / screen_dpi;
        H_in = pos(4) / screen_dpi;
        margin = 0.05;

        set(fh, 'PaperUnits', 'inches');
        set(fh, 'PaperPositionMode', 'manual');
        set(fh, 'PaperPosition', [margin, margin, W_in, H_in]);
        set(fh, 'PaperSize', [W_in + 2*margin, H_in + 2*margin]);
        try, set(fh, 'InvertHardcopy', 'off'); catch, end
        drawnow();
        % --------------------------------

        for k = 1:numel(formats)
            index = index+1;

            path = fullfile(out_dir,[name,'.',formats{k}]);

            item.figure_name = name;
            item.format = formats{k};
            item.path = path;
            item.success = false;
            item.message = '';

            try
                if strcmp(formats{k},'png')
                    print(fh,path,devices{k},'-r300');
                elseif strcmp(formats{k},'eps')
                    print(fh,path,devices{k},'-painters'); % EPS 需要矢量渲染
                else
                    print(fh,path,devices{k});
                end

                file_info = dir(path);

                if isempty(file_info) || file_info.bytes == 0
                    error('No nonempty export file was produced.');
                end

                item.success = true;

                fprintf('  Exported: %s\n',path);

            catch ME
                item.message = ME.message;

                fprintf('  Export failed: %s\n',path);
                fprintf('    %s\n',ME.message);
            end

            Report(index) = item;
        end
    end

    if isempty(Report)
        fprintf('export_all_figures: no figures supplied.\n');
    else
        fprintf('export_all_figures: %d/%d exports succeeded.\n', ...
            sum([Report.success]),numel(Report));
    end
end
