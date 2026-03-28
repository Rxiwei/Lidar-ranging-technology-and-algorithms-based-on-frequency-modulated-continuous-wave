clear;
%% 步骤1：读取数据
tic;
fileID = fopen('实验数据/data_200m_20dB.txt', 'r');

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

%% 2. 初始化参数
T_chirp=10e-6;                      % 调频周期
c=3e8;                              % 光速
B=300e6;                            % 调频带宽

% 裁剪实际读取量
time = time(1:line_count);
amplitude = amplitude(1:line_count);
fclose(fileID);

    %% 3. 小波去噪函数
function denoised = wavelet_denoise(signal, wavelet, level)
    % 小波分解
    [C, L] = wavedec(signal, level, wavelet);
    
    % 估计噪声标准差（使用第1层细节系数）
    detail_coeffs = detcoef(C, L, 1);
    sigma = median(abs(detail_coeffs)) / 0.6745;
    
    % 计算通用阈值
    N = length(signal);
    threshold = sigma * sqrt(2*log(N));
    
    % 阈值处理（软阈值）
    C_thresh = wthresh(C, 's', threshold);
    
    % 重构信号
    denoised = waverec(C_thresh, L, wavelet);
    % 
    % 调整信号长度
    if length(denoised) > length(signal)
        denoised = denoised(1:length(signal));
    end
end
%% 4. Zoom-FFT实现
function [f_zoom, X_zoom] = zoom_fft(x, fs, f_center, D)
    % ZOOM_FFT 实现基于频移法的Zoom-FFT频谱细化
    % 输入参数：
    %   x: 输入信号（向量）
    %   fs: 原始采样率（Hz）
    %   f_center: 目标中心频率（Hz）
    %   D: 下采样倍数（整数）
    % 输出参数：
    %   f_zoom: 细化后的频率轴（Hz）
    %   X_zoom: 细化后的频谱幅值

    % 参数校验
    validateattributes(x, {'double'}, {'vector'}, 'zoom_fft', 'x', 1);
    validateattributes(fs, {'double'}, {'positive'}, 'zoom_fft', 'fs', 2);
    validateattributes(f_center, {'double'}, {'nonnegative', '<', fs/2}, 'zoom_fft', 'f_center', 3);
    validateattributes(D, {'double'}, {'integer', '>=', 1}, 'zoom_fft', 'D', 4);
    
    % 核心参数计算
    N = length(x);               % 原始信号长度
    BW = fs / (2 * D);           % 理论最大带宽（Nyquist约束）
    fs_new = fs / D;             % 下采样后的新采样率
    
    % 复调制（频谱搬移）
    t = (0:N-1) / fs;           
    x_mod = x .* exp(-1j * 2 * pi * f_center * t(:)); 
    
    % 抗混叠滤波器设计
    Wn = BW / (fs/2);            
    [b, a] = butter(2, Wn);      
    
    % 零相位滤波
    x_filt = filtfilt(b, a, x_mod);
    
    % 下采样
    x_down = x_filt(1:D:end);    
    N_zoom = length(x_down);     
    
    % FFT分析与频率轴生成
    % 补零至4的幂次长度（提升计算效率）
    M = 2^nextpow2(N_zoom * 4);   
    X_zoom = fft(x_down, M);     
    X_zoom = abs(X_zoom) / N_zoom * 2; 
    
    % 生成精确频率轴
    f_resolution = fs_new / M;  
    f_base = (-fs_new/2 : f_resolution : fs_new/2 - f_resolution); 
    f_zoom = f_center + f_base;  
    
    % 截取有效频段（可选）
    keep_idx = (f_zoom >= f_center - BW) & (f_zoom <= f_center + BW);
    f_zoom = f_zoom(keep_idx);
    X_zoom = X_zoom(keep_idx);
end


%% 5. FFT变换
% 基础FFT参数设置
N = length(amplitude);
dt=(time(2)-time(1))*1e-5;          % 时间间隔
fs=1/dt;                            % 采样频率
signal=amplitude-mean(amplitude);   % 去直流分量
% 应用小波去噪
%signal_denoised = wavelet_denoise(signal, 'db4', 3);
% 加窗
window=hann(N);
signal_windowed=signal.*window;

% 基本FFT分析
signal_fft=fft(signal_windowed); 
P2 = abs(signal_fft / N);                    % 归一化双边幅度谱
P1 = P2(1:N/2+1);                 % 单边频谱
% 生成频率轴
f_signal = fs * (0:N/2) / N;
% 定位频谱最大值
[~,idx_base]=max(P1);
% 频谱峰值频率
f_direct=f_signal(idx_base);

% 设置Zoom参数
D = 420;           % 下采样倍数
f_center = f_direct; % 以初步估计值为中心
% 执行Zoom-FFT
[f_zoom, X_zoom] = zoom_fft(signal_fft, fs, f_center, D);

% 找到精确峰值
[~, idx_zoom] = max(abs(X_zoom));
max_freq = f_zoom(idx_zoom);

%% 6. 输出结果
% 计算距离
R_direct=(c * T_chirp * f_direct) / (2 * B);
R_zoom=(c * T_chirp * max_freq) / (2 * B);
toc;

disp(['运行时间: ',num2str(toc)]);
fprintf('=============== 测距结果对比 ===============\n');
fprintf('直接FFT法:   峰值频率 = %.5f Hz, 距离 = %.5f m\n', f_direct, R_direct);
fprintf('Zoom-FFT法:   峰值频率 = %.5f Hz, 距离 = %.5f m\n', max_freq, R_zoom);

% 绘制Zoom-FFT细化频谱
figure(1);
plot(f_zoom, X_zoom/max(X_zoom));
grid on;
xlabel('Frequency (Hz)'), ylabel('Amplitude (dB)');
