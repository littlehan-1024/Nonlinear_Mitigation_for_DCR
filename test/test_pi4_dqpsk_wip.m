%% π/4-DQPSK信号非线性失真抑制测试
clear; close all; clc;
addpath(genpath(pwd));

%% 1. 基本参数设置
fs = 25e6;          % 采样率 25 MHz
fc = 2.6e6;         % 载波频率 2.6 MHz
Rs = 1e6;           % 符号率 1 MHz
Rb = 2e6;           % 比特率 2 Mbps (DQPSK每个符号2比特)
rolloff = 0.35;     % 滚降系数
N_symbols = 500;    % 符号数
sps = fs / Rs;      % 每符号采样点数

% 非线性参数
alpha = struct();
alpha.a1 = 5.62;
alpha.a2 = -(84351 + 1j*74391);
alpha.a3 = 3.16;
alpha.a4 = -1588.7;

% I/Q不平衡参数
iq_params.gm = 0.99;
iq_params.phi_m = deg2rad(3.6);

% NLMS算法参数
M = 5;
mu = [1, 1, 0.01, 1, 1];
alpha_nlms = [1e-9, 1e-8, 1e-4, 1e-9, 1e-8];

% 滤波器参数（根据信号带宽设置）
bandwidth = Rs * (1 + rolloff);  % 信号带宽
f_low = fc - bandwidth/2 - 0.2e6;
f_high = fc + bandwidth/2 + 0.2e6;
f1 = f_low;
f2 = f_high;

%% 2. 生成π/4-DQPSK信号
fprintf('生成π/4-DQPSK信号...\n');

% 2.1 生成随机比特序列
N_bits = N_symbols * 2;  % DQPSK每个符号2比特
data_bits = randi([0, 1], 1, N_bits);

% 2.2 π/4-DQPSK星座映射（差分编码）
% 星座点：±45°, ±135°
constellation = exp(1j * pi/4 * [1, 3, -3, -1]);  % 45°, 135°, -135°, -45°
symbol_map = [0 0; 0 1; 1 1; 1 0];  % 格雷码映射

% 比特到符号映射
symbol_indices = zeros(1, N_symbols);
for i = 1:N_symbols
    bits = data_bits(2*i-1:2*i);
    for j = 1:4
        if isequal(bits, symbol_map(j, :))
            symbol_indices(i) = j;
            break;
        end
    end
end
initial_symbols = constellation(symbol_indices);

% 2.3 差分编码
diff_symbols = zeros(1, N_symbols);
prev_symbol = 1;  % 初始相位
for i = 1:N_symbols
    diff_symbols(i) = prev_symbol * initial_symbols(i);
    prev_symbol = diff_symbols(i);
end

% 2.4 脉冲成形（升余弦滤波）
span = 10;  % 滤波器跨度
rrc_filter = rcosdesign(rolloff, span, sps, 'sqrt');

% 过采样
upsampled_symbols = zeros(1, N_symbols * sps);
upsampled_symbols(1:sps:end) = diff_symbols;

% 脉冲成形
shaped_signal = conv(upsampled_symbols, rrc_filter, 'same');

% 2.5 上变频到载波频率
t = (0:length(shaped_signal)-1) / fs;
pi4_dqpsk_signal = real(shaped_signal) .* cos(2*pi*fc*t) - ...
                   imag(shaped_signal) .* sin(2*pi*fc*t);

% 显示原始信号时域波形
figure(1);
subplot(3,2,1);
plot(t(1:min(500, length(t)))*1e6, pi4_dqpsk_signal(1:min(500, length(t))));
xlabel('时间 (μs)'); ylabel('幅度'); title('π/4-DQPSK时域波形');
grid on;

subplot(3,2,2);
plot(real(shaped_signal(1:min(500, length(shaped_signal)))), ...
     imag(shaped_signal(1:min(500, length(shaped_signal)))), '.');
xlabel('同相分量'); ylabel('正交分量'); title('基带信号星座图');
axis equal; grid on;

%% 3. 计算原始信号频谱
N_fft = 8192;
f_axis = (-N_fft/2:N_fft/2-1) * fs / N_fft / 1e6;

subplot(3,2,[3,4]);
Y_original = fftshift(fft(pi4_dqpsk_signal, N_fft));
plot(f_axis, 20*log10(abs(Y_original)+1e-12), 'b', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('原始π/4-DQPSK信号频谱');
xlim([-15, 15]); ylim([-100, 50]); grid on;
hold on;
plot([fc/1e6, fc/1e6], [-100, 50], 'r--');
plot([-fc/1e6, -fc/1e6], [-100, 50], 'r--');
legend('信号', '载波频率', 'Location', 'best');

%% 4. 非线性失真建模
fprintf('模拟非线性失真...\n');
[y_received, y_t, y_IQ] = nonlinear_distortion_model(pi4_dqpsk_signal, alpha, iq_params, ...
    'SNR_dB', 61, 'P_in_dBm', -30);

% 显示失真后星座图
subplot(3,2,5);
% 下变频
t_rec = (0:length(y_received)-1) / fs;
I_rec = y_received .* cos(2*pi*fc*t_rec) * 2;
Q_rec = -y_received .* sin(2*pi*fc*t_rec) * 2;

% 低通滤波
[b_lp, a_lp] = butter(6, (Rs * (1+rolloff))/(fs/2));
I_filtered = filter(b_lp, a_lp, I_rec);
Q_filtered = filter(b_lp, a_lp, Q_rec);

% 抽取
decimation_factor = floor(sps/2);
I_decimated = I_filtered(1:decimation_factor:end);
Q_decimated = Q_filtered(1:decimation_factor:end);

plot(I_decimated(101:200), Q_decimated(101:200), '.');
xlabel('同相分量'); ylabel('正交分量'); title('失真后星座图（采样后）');
axis equal; grid on;

%% 5. 显示失真频谱
subplot(3,2,6);
Y_distorted = fftshift(fft(y_received, N_fft));
plot(f_axis, 20*log10(abs(Y_distorted)+1e-12), 'r', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('非线性失真后频谱');
xlim([-15, 15]); ylim([-120, 50]); grid on;
hold on;
plot([fc/1e6, fc/1e6], [-120, 50], 'k--');
plot([-fc/1e6, -fc/1e6], [-120, 50], 'k--');
legend('失真信号', '载波频率', 'Location', 'best');

%% 6. 非线性失真抑制
fprintf('进行非线性失真抑制...\n');
[d, x_hat, y_final, w, e, avg_suppression, suppression_results] = ...
    nonlinear_cancellation(y_received, fs, f_low, f_high, M, mu, alpha_nlms, 0, 0, ...
    'signal_type', 'modulated', 'fc', fc, 'N_fft', N_fft);

%% 7. 显示抑制结果
figure(2);
set(gcf, 'Position', [100, 100, 1200, 800]);

% 7.1 抑制前后频谱对比
subplot(3,2,[1,2]);
hold on; grid on;
Y_before = fftshift(fft(y_received, N_fft));
Y_after = fftshift(fft(y_final, N_fft));

h1 = plot(f_axis, 20*log10(abs(Y_before)+1e-12), 'b-', 'LineWidth', 1.5, 'DisplayName', '抑制前');
h2 = plot(f_axis, 20*log10(abs(Y_after)+1e-12), 'r-', 'LineWidth', 2.0, 'DisplayName', '抑制后');
plot([fc/1e6, fc/1e6], [-120, 50], 'k--', 'LineWidth', 1.0);
plot([-fc/1e6, -fc/1e6], [-120, 50], 'k--', 'LineWidth', 1.0);

xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); 
title('π/4-DQPSK信号抑制前后频谱对比');
xlim([-15, 15]); ylim([-120, 50]);
legend('Location', 'best');

% 7.2 局部放大
subplot(3,2,3);
hold on; grid on;
plot(f_axis, 20*log10(abs(Y_before)+1e-12), 'b-', 'LineWidth', 1.0);
plot(f_axis, 20*log10(abs(Y_after)+1e-12), 'r-', 'LineWidth', 1.5);
plot([fc/1e6, fc/1e6], [-120, 50], 'k--');

xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); 
title('正频率部分局部放大');
xlim([1, 5]); ylim([-100, 50]);

% 7.3 抑制后星座图
subplot(3,2,4);
% 对抑制后信号下变频
I_final = y_final .* cos(2*pi*fc*t_rec) * 2;
Q_final = -y_final .* sin(2*pi*fc*t_rec) * 2;

I_final_filtered = filter(b_lp, a_lp, I_final);
Q_final_filtered = filter(b_lp, a_lp, Q_final);

I_final_decimated = I_final_filtered(1:decimation_factor:end);
Q_final_decimated = Q_final_filtered(1:decimation_factor:end);

plot(I_final_decimated(101:200), Q_final_decimated(101:200), 'r.', 'MarkerSize', 8);
xlabel('同相分量'); ylabel('正交分量'); title('抑制后星座图');
axis equal; grid on;

% 7.4 误差收敛曲线
subplot(3,2,5);
e_power = 10*log10(abs(e).^2 + 1e-12);
plot(1:length(e_power), e_power, 'b-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('误差功率 (dB)'); title('NLMS算法误差收敛曲线');
grid on;

% 7.5 滤波器系数收敛
subplot(3,2,6);
w_norm = zeros(size(w, 3), 1);
for n = 1:size(w, 3)
    w_norm(n) = norm(w(:,:,n), 'fro');
end
plot(1:length(w_norm), 10*log10(w_norm + 1e-12), 'g-', 'LineWidth', 1.5);
xlabel('迭代次数'); ylabel('滤波器系数范数 (dB)'); title('滤波器系数收敛');
grid on;

%% 8. 加窗细节分析
fprintf('进行加窗细节分析...\n');
figure(3);
set(gcf, 'Position', [100, 100, 1200, 400]);

N_window = 2048;
start_idx = 1000;

y_before_segment = y_received(start_idx:start_idx+N_window-1);
y_after_segment = y_final(start_idx:start_idx+N_window-1);

hann_window = hann(N_window)';
y_before_windowed = y_before_segment .* hann_window;
y_after_windowed = y_after_segment .* hann_window;

N_fft_local = 2048;
f_axis_local = (-N_fft_local/2:N_fft_local/2-1)*fs/N_fft_local/1e6;

Y_before_windowed_fft = fftshift(fft(y_before_windowed, N_fft_local));
Y_after_windowed_fft = fftshift(fft(y_after_windowed, N_fft_local));

hold on;
h1 = plot(f_axis_local, 20*log10(abs(Y_before_windowed_fft)+1e-12), ...
    'b-', 'LineWidth', 1.5, 'DisplayName', '抑制前');
h2 = plot(f_axis_local, 20*log10(abs(Y_after_windowed_fft)+1e-12), ...
    'r-', 'LineWidth', 2.0, 'DisplayName', '抑制后');

plot([fc/1e6, fc/1e6], [-150, 50], 'k--', 'LineWidth', 1.0);
plot([-fc/1e6, -fc/1e6], [-150, 50], 'k--', 'LineWidth', 1.0);
plot([f_low/1e6, f_low/1e6], [-150, 50], 'g:', 'LineWidth', 1.0);
plot([f_high/1e6, f_high/1e6], [-150, 50], 'g:', 'LineWidth', 1.0);

text(fc/1e6 + 0.2, 30, sprintf('fc=%.1fMHz', fc/1e6), 'FontSize', 10, 'Color', 'k');
text(f_low/1e6 - 0.5, -100, sprintf('%.1fMHz', f_low/1e6), 'FontSize', 9, 'Color', 'g');
text(f_high/1e6 + 0.2, -100, sprintf('%.1fMHz', f_high/1e6), 'FontSize', 9, 'Color', 'g');

xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
title('π/4-DQPSK信号抑制前后细节对比（加窗分析）');
xlim([-15, 15]); ylim([-150, 50]); grid on;
legend('Location', 'best', 'FontSize', 10); box on;

%% 9. 性能评估
fprintf('\n==================== π/4-DQPSK信号测试报告 ====================\n');
fprintf('信号参数:\n');
fprintf('   采样率: %.1f MHz\n', fs/1e6);
fprintf('   载波频率: %.1f MHz\n', fc/1e6);
fprintf('   符号率: %.1f Msps\n', Rs/1e6);
fprintf('   比特率: %.1f Mbps\n', Rb/1e6);
fprintf('   滚降系数: %.2f\n', rolloff);
fprintf('   信号带宽: %.1f MHz\n', bandwidth/1e6);
fprintf('   滤波器通带: %.1f-%.1f MHz\n\n', f_low/1e6, f_high/1e6);

fprintf('非线性失真抑制性能:\n');
fprintf('   平均带外抑制比: %.1f dB\n', avg_suppression);
if isscalar(suppression_results)
    fprintf('   带内信号功率比: %.1f dB\n', suppression_results);
elseif length(suppression_results) >= 2
    fprintf('   带内信号质量改善: %.1f dB\n', suppression_results(1));
    fprintf('   带外抑制比: %.1f dB\n', suppression_results(2));
end
fprintf('\n');

fprintf('NLMS算法配置:\n');
fprintf('   滤波器阶数: M = %d\n', M);
fprintf('   步长参数: [%.2f, %.2f, %.4f, %.2f, %.2f]\n', mu);
fprintf('   正则化参数: [%.1e, %.1e, %.1e, %.1e, %.1e]\n', alpha_nlms);
fprintf('===============================================================\n');

%% 10. 计算并显示EVM（误差矢量幅度）
fprintf('\n计算误差矢量幅度(EVM)...\n');

% 对抑制后信号进行解调
I_final_decimated_aligned = I_final_decimated(span+1:end-span);
Q_final_decimated_aligned = Q_final_decimated(span+1:end-span);

% 判决
received_symbols = I_final_decimated_aligned + 1j*Q_final_decimated_aligned;
detected_symbols = zeros(1, length(received_symbols));
for i = 1:length(received_symbols)
    % 找到最近的星座点
    distances = abs(received_symbols(i) - constellation);
    [~, idx] = min(distances);
    detected_symbols(i) = constellation(idx);
end

% 计算EVM
evm_rms = sqrt(mean(abs(received_symbols - detected_symbols).^2)) / ...
          sqrt(mean(abs(detected_symbols).^2)) * 100;

fprintf('   抑制后信号EVM: %.2f%%\n', evm_rms);

% 对原始信号也计算EVM（作为对比）
I_original = I_decimated(span+1:end-span);
Q_original = Q_decimated(span+1:end-span);
original_symbols = I_original + 1j*Q_original;

detected_original = zeros(1, length(original_symbols));
for i = 1:length(original_symbols)
    distances = abs(original_symbols(i) - constellation);
    [~, idx] = min(distances);
    detected_original(i) = constellation(idx);
end

evm_original = sqrt(mean(abs(original_symbols - detected_original).^2)) / ...
               sqrt(mean(abs(detected_original).^2)) * 100;
fprintf('   原始信号EVM: %.2f%%\n', evm_original);
fprintf('   EVM改善: %.2f%%\n', evm_original - evm_rms);

%% 11. 保存结果
save('pi4_dqpsk_results.mat', 'y_received', 'y_final', 'avg_suppression', ...
    'suppression_results', 'fc', 'Rs', 'Rb', 'fs', 'M', 'evm_rms', 'evm_original');
fprintf('\n结果已保存至 pi4_dqpsk_results.mat\n');
fprintf('测试完成！\n');