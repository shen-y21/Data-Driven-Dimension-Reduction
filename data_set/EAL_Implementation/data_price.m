%% 读取PJM电力市场电价数据
% 与STN和RTN实现使用相同的数据源
% 读取2022年7月PJM系统的实时电价数据
% 支持数据缓存，避免重复读取Excel文件

%% 检查是否已有缓存数据
cache_file = 'price_data_cache.mat';
if exist(cache_file, 'file')
    fprintf('发现电价数据缓存文件，正在加载...\n');
    load(cache_file);
    fprintf('✓ 电价数据从缓存加载完成\n');
    fprintf('数据维度: %d x %d (小时 x 天数)\n', size(Price_days, 1), size(Price_days, 2));
    fprintf('电价范围: %.4f - %.4f 元/kWh\n', min(Price_days(:)), max(Price_days(:)));
    fprintf('数据来源: %s\n', data_source);
    fprintf('缓存时间: %s\n', cache_time);
else
    fprintf('未发现缓存文件，正在从Excel读取电价数据...\n');
    
    %% 市场参数
    % 时间间隔长度，1小时
    delta_t = 1;
    day_price = 1; % 选择的天数
    hour_init = 1; % 从这一天第一个时间间隔开始
    NOFSLOTS = 24 / delta_t;

    % 读取这个月（7月）所有天的电价
    Price_days = [];
    filename = '../STN_Implementation/rt_hrl_lmps.xlsx';
    
    % 检查Excel文件是否存在
    if ~exist(filename, 'file')
        error('找不到电价数据文件: %s\n请确保文件路径正确', filename);
    end
    
    fprintf('正在读取Excel文件: %s\n', filename);
    for day_price = 1 : 31
        start_row = (day_price - 1) * 24 + hour_init + 1; % 起始行
        sheet = 'rt_hrl_lmps'; % 工作表名称
        xlRange = "I" + start_row + ":I" + (start_row + NOFSLOTS - 1); % 范围
        Price = xlsread(filename, sheet, xlRange); % 读取电价
        Price_days = [Price_days, Price];
        
        % 显示进度
        if mod(day_price, 5) == 0 || day_price == 31
            fprintf('  已读取 %d/31 天数据\n', day_price);
        end
    end

    % 清除变量
    clear Price filename sheet xlRange start_row hour_init day_price NOFSLOTS delta_t

    % 将电价转换为$/kWh，然后转换为人民币
    Price_days = Price_days * 1e-3;  % 转换为$/kWh
    
    % 汇率转换：从美元转换为人民币
    if exist('EXCHANGE_RATE', 'var')
        exchange_rate = EXCHANGE_RATE;
    else
        % 默认汇率：1美元 ≈ 7.2人民币
        exchange_rate = 7.2;
        fprintf('使用默认汇率: 1美元 = %.1f人民币\n', exchange_rate);
        fprintf('提示: 可运行 exchange_rate_config.m 设置自定义汇率\n');
    end
    Price_days = Price_days * exchange_rate;  % 转换为元/kWh
    
    % 保存缓存数据
    data_source = '../STN_Implementation/rt_hrl_lmps.xlsx';
    cache_time = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    save(cache_file, 'Price_days', 'data_source', 'cache_time');
    
    fprintf('✓ 电价数据读取完成并已缓存\n');
    fprintf('数据维度: %d x %d (小时 x 天数)\n', size(Price_days, 1), size(Price_days, 2));
    fprintf('电价范围: %.4f - %.4f 元/kWh\n', min(Price_days(:)), max(Price_days(:)));
    fprintf('汇率: 1美元 = %.1f人民币\n', exchange_rate);
    fprintf('缓存文件: %s\n', cache_file);
end
