% 清空内存

clear;
tic;
fileID = fopen('实验数据/data_200m_10dB.txt', 'r');

% 预分配内存（提升性能）
max_lines = 8192;
time = zeros(max_lines, 1);
amplitude = zeros(max_lines, 1);

% 定义数据格式（根据实际分隔符调整）
delimiter = ' '; % 可能的分隔符：',', '\t', ';'等
formatSpec = ['%f' delimiter '%f'];

% 逐行读取
line_count = 0;
while line_count < max_lines
    tline = fgetl(fileID);
    if ~ischar(tline), break; end % 文件结束检查
    
    % 解析当前行
    parsed = sscanf(tline, formatSpec);
    if numel(parsed) == 2
        line_count = line_count + 1;
        time(line_count) = parsed(1);
        amplitude(line_count) = parsed(2);
    else
        warning('第%d行格式错误: %s', line_count+1, tline);
    end
end

% 初始化参数
T_chirp=10e-6;                      % 调频周期
c=3e8;                              % 光速
B=300e6;                            % 调频带宽

% 裁剪实际读取量
time = time(1:line_count);
amplitude = amplitude(1:line_count);
fclose(fileID);

% 基础FFT参数设置
N = length(amplitude);
dt=(time(2)-time(1))*1e-5;                 % 时间间隔
fs=1/dt;                            % 采样频率
signal=amplitude-mean(amplitude);   % 去直流分量

% 加窗
window=hann(N);
signal_windowed=signal.*window;

%% ====================== 方法1: 能量重心法 ======================
% 基本FFT分析
signal_energy=fft(signal_windowed); 
P2_energy = abs(signal_energy / N);             % 归一化双边幅度谱
P1_energy = P2_energy(1:N/2+1);                 % 单边频谱
% 生成频率轴
f_energy = fs * (0:N/2) / N;
% 定位频谱最大值
[~,idx_base]=max(P1_energy);

% 核心算法
% 设置能量窗口范围
half_win = 9;
win_start = max(1, idx_base - half_win);
win_end = min(length(P1_energy), idx_base + half_win);
window_indices = win_start:win_end;
% 计算能量重心
energy_sum = sum(P1_energy(window_indices).^2);
weighted_sum = sum(f_energy(window_indices)*(P1_energy(window_indices).^2));
f_peak_centroid = weighted_sum / energy_sum;
% 计算距离
R_energy=(c * T_chirp * f_peak_centroid) / (2 * B);

%% ====================== 方法2: 最大值估计法 ======================
function corrected_freq = max_spectrum_correction(signal, fs, de, max_iter)
% 最大值估计法频谱校正函数
% 参数:
%   signal   - 输入的中频信号 (时域)
%   fs       - 采样频率 (Hz)
%   de       - 允许的最大频率误差阈值 (Hz)
%   max_iter - 最大迭代次数 (默认100)
% 返回:
%   corrected_freq - 校正后的频率估计值 (Hz)

if nargin < 4
    max_iter = 100; % 默认最大迭代次数
end

n = length(signal);
fft_result = fft(signal);
freq_res = fs/n;
freq_axis = (0:n-1)*freq_res; % 生成正频率轴

% 转换为双边谱
if mod(n,2) == 0
    freq_axis(n/2+1:end) = freq_axis(n/2+1:end) - fs;
else
    freq_axis((n+1)/2+1:end) = freq_axis((n+1)/2+1:end) - fs;
end

magnitude = abs(fft_result);

% 找到主瓣内的最大值和次大值谱线
[~, max_idx] = max(magnitude);
search_range = 3; % 主瓣搜索范围（根据窗函数特性调整）
% 确定主瓣区间
start_idx = max(1, max_idx - search_range);
end_idx = min(n, max_idx + search_range);
mainlobe = start_idx:end_idx;
% 在主瓣内寻找次大值
[~, sorted_idx] = sort(magnitude(mainlobe), 'descend');
second_max_rel_idx = sorted_idx(2);
second_max_idx = mainlobe(1) + second_max_rel_idx - 1;
% 初始化参数
ft = freq_axis(max_idx);
fx = freq_axis(second_max_idx);
k1 = magnitude(max_idx);
k2 = magnitude(second_max_idx);
% 获取初始谱线间隔
delta_f = freq_axis(2) - freq_axis(1);
for iter = 1:max_iter
    % 计算斜率参数
    m = (k2 - k1) / (k1 + eps);
    if abs(m - 2) < eps
        break;
    end
    % 计算频率偏差
    e = ((m-1)*(fx - ft)) / (m - 2);
    
    % 动态调整步长（与谱线间隔关联）
    adaptive_de = min(de, delta_f/2); % 确保不超过半根谱线间隔
    
    if abs(e) <= adaptive_de
        break;
    end
    % 调整谱线位置并更新幅度值
    delta = sign(e) * adaptive_de;
    ft = ft + delta;
    fx = fx + delta;
    % 更新幅度值（通过插值获取新位置的幅度）
    [~, idx] = min(abs(freq_axis - ft));
    k1 = magnitude(idx);
    [~, idx] = min(abs(freq_axis - fx));
    k2 = magnitude(idx);
end
corrected_freq = ft;
end

%f_peak_highres=max_spectrum_correction(signal_windowed, fs,300,50);
%R_highres=(c * T_chirp * f_peak_highres) / (2 * B);
toc;
%% 结果对比输出
disp(['运行时间: ',num2str(toc)]);
fprintf('=============== 测距结果对比 ===============\n');
%fprintf('能量重心法:   峰值频率 = %.5f Hz, 距离 = %.5f m\n', f_peak_centroid, R_energy);
%fprintf('最大值估计法: 峰值频率 = %.5f Hz, 距离 = %.5f m\n', f_peak_highres, R_highres);
%{
% 绘制原始信号时域图
figure;
plot(time,amplitude,'LineWidth',0.05);
grid on;
xlabel("time(s)","FontSize",14);
ylabel("amplitude","FontSize",14);
legend('time-domain','FontSize',14);

% 能量重心法频域绘图
figure(2);
plot(f_energy, P1_energy);
xlim([0 f_peak_centroid*2])
grid on;
xlabel('Frequency (Hz)', 'FontSize', 14);
ylabel('Magnitude', 'FontSize', 14);
%}
