function [target_dates,ib_dates_set] = IB_date(ib_parent_path,years)

    ib_files = dir(fullfile(ib_parent_path, '*.txt'));
    datePattern = '_\d{8}';
    target_dates = {};
    for i = 1:numel(ib_files)
        name = ib_files(i).name;
        m = regexp(name, datePattern, 'match');
        if ~isempty(m)
            ds = m{1}(2:end);
            yy = str2double(ds(1:4)); 
            mm = str2double(ds(5:6));
            % 由于era辅助数据没有4月数据，因此将四月剔除
            if yy > 2010 && ismember(mm, 3)
                target_dates{end+1} = ds;
            end
%             if yy > 2010 && ismember(mm, [3,4])
%                 target_dates{end+1} = ds;
%             end
        end
    end

    %% Group dates by year
    ib_dates_set = cell(1, numel(years));
    for idx = 1:numel(years)
        y = years(idx);
        list = [];
        for k = 1:numel(target_dates)
            ds = target_dates{k};
            date = datetime(ds, 'InputFormat', 'yyyyMMdd');  % 转换为 datetime
            target_day = day(date,'dayofyear');
            end_date = datetime([num2str(y),'0331'], 'InputFormat', 'yyyyMMdd');  % 转换为 datetime
%             end_date = datetime([num2str(y),'0415'], 'InputFormat', 'yyyyMMdd');  % 转换为 datetime
            rule_day = day(end_date, 'dayofyear');  % 计算年积日

            if str2double(ds(1:4)) == y && target_day <= rule_day
                list(end+1) = str2double(ds);
            end
        end
        ib_dates_set{idx} = list;
    end
end

