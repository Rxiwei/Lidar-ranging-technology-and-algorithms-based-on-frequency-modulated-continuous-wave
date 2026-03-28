OutputPort1 = InputPort1;   %创建输出端口数据结构
%获取样本数目
SampledNumber=length(InputPort1.Sampled);
%测距距离
R=200; 
%衰减参数
attenution = sqrt(exp(-0.01178*R));
% 设置目标信噪比（-10 dB）
SNR_dB = 0;

if(SampledNumber>0)
    for i=1:SampledNumber
        % 衰减
        OutputPort1.Sampled(i).Signal(1,:)=InputPort1.Sampled(i).Signal(1,:)*attenution;
        % 添加高斯白噪声
        % 提取当前信号（假设为复数光场复包络）
        signal = OutputPort1.Sampled(i).Signal(1,:);
        % 计算信号功率（复数信号功率公式）
        Ps = mean(abs(signal).^2);
        % 根据SNR计算噪声功率（SNR = Ps/Pn）
        Pn = Ps / (10^(SNR_dB/10));
        % 生成高斯白噪声（总功率为Pn）
        Noise = sqrt(Pn/2) * (randn(size(signal)) + 1j*randn(size(signal)));  % 正确复数噪声
        % 添加瑞利分布杂波（示例）
        clutter_power = Ps * 10^(-5); 
        clutter = sqrt(clutter_power/2) * (randn(size(signal)) + 1j*randn(size(signal)));
        % 将噪声叠加到信号
        OutputPort1.Sampled(i).Signal(1,:) = signal + Noise + clutter;
    end
end