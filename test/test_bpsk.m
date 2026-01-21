%% BPSK信号非线性失真抑制测试
clear; close all; clc;
addpath(genpath(pwd));

%% 1. 基本参数设置
fs = 25e6;          % 采样率
fc = 2.6e6;         % 载波频率
Rb = 788e3;         % 比特率
rolloff = 0.5;      % 滚降系数
N_bits = 1000;      % 比特数
f1 = fc - Rb/2;
f2 = fc + Rb/2;

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
f_low = fc - Rb/2 - 0.1e6;
f_high = fc + Rb/2 + 0.1e6;

%% 2. 生成BPSK信号
fprintf('生成BPSK信号...\n');
data_bits = randi([0, 1], 1, N_bits);
bpsk_symbols = 2*data_bits - 1;  % 0->-1, 1->+1

samples_per_symbol = N_bits * floor(fs / Rb);
upsampled_signal = zeros(1, samples_per_symbol);
for i = 1:N_bits
    upsampled_signal((((i-1)*samples_per_symbol)/N_bits) + 1) = bpsk_symbols(i);
end

span = 10;
rcos_filter = rcosdesign(rolloff, span, samples_per_symbol/N_bits, 'normal');
baseband_signal = conv(upsampled_signal, rcos_filter, 'same');

t = (0:length(baseband_signal)-1) / fs;
carrier = cos(2*pi*fc*t);
bpsk_signal = baseband_signal .* carrier;

N_fft = min(8192, length(bpsk_signal));
f_axis = (-N_fft/2:N_fft/2-1)*fs/N_fft/1e6;

figure(1);
plot(f_axis, 20*log10(abs(fftshift(fft(bpsk_signal, N_fft))))+1e-12);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('原始bpsk信号频谱');
xlim([-15, 15]); ylim([-100, 50]); grid on;
hold on;
plot([f1/1e6, f1/1e6], [-100, 50], 'r--');
plot([f2/1e6, f2/1e6], [-100, 50], 'r--');
legend('信号', 'f1', 'f2');

%% 3. 非线性失真建模
fprintf('模拟非线性失真...\n');
[y_received, y_t, y_IQ] = nonlinear_distortion_model(bpsk_signal, alpha, iq_params, ...
    'SNR_dB', 61, 'P_in_dBm', -30);

%% 4. 非线性失真抑制
fprintf('进行非线性失真抑制...\n');
[d, x_hat, y_final, w, e, avg_suppression, suppression_results] = ...
    nonlinear_cancellation(y_received, fs, f_low, f_high, M, mu, alpha_nlms, 0, 0, ...
    'signal_type', 'bpsk', 'fc', fc, 'N_fft', 8192);

%% 5. 显示结果

figure(2);
subplot(2,1,1);
Y_before_fft = fftshift(fft(y_received, N_fft));
plot(f_axis, 20*log10(abs(Y_before_fft)+1e-12), 'b', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('BPSK信号抑制前频谱');
xlim([-15, 15]); ylim([-120, 50]); grid on;
hold on;
plot([fc/1e6, fc/1e6], [-120, 50], 'r--');
legend('信号', '中心频率');

subplot(2,1,2);
Y_final_fft = fftshift(fft(y_final, N_fft));
plot(f_axis, 20*log10(abs(Y_final_fft)+1e-12), 'g', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('BPSK信号抑制后频谱');
xlim([-15, 15]); ylim([-120, 50]); grid on;
hold on;
plot([fc/1e6, fc/1e6], [-120, 50], 'r--');
legend('抑制后信号', '中心频率');

%% 6. 对比图
figure(3);
plot(f_axis, 20*log10(abs(Y_before_fft)+1e-12), 'b', 'LineWidth', 1.5);
hold on;
plot(f_axis, 20*log10(abs(Y_final_fft)+1e-12), 'g', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('BPSK信号抑制前后频谱对比');
xlim([-15, 15]); ylim([-120, 50]); grid on;
legend('抑制前', '抑制后', 'Location', 'best');

%% 7. 加窗细节分析
fprintf('进行加窗细节分析...\n');
figure('Position', [100, 100, 1200, 400]);

N_window = 2500;
start_idx = 1;

y_before_segment = y_received(start_idx:start_idx+N_window-1);
y_after_segment = y_final(start_idx:start_idx+N_window-1);

hann_window = hann(N_window)';
y_before_windowed = y_before_segment .* hann_window;
y_after_windowed = y_after_segment .* hann_window;

N_fft_local = 8192;
f_axis_local = (-N_fft_local/2:N_fft_local/2-1)*fs/N_fft_local/1e6;

Y_before_windowed_fft = fftshift(fft(y_before_windowed, N_fft_local));
Y_after_windowed_fft = fftshift(fft(y_after_windowed, N_fft_local));

hold on;
h1 = plot(f_axis_local, 20*log10(abs(Y_before_windowed_fft)+1e-12), ...
    'b-', 'LineWidth', 1.5, 'DisplayName', '抑制前');
h2 = plot(f_axis_local, 20*log10(abs(Y_after_windowed_fft)+1e-12), ...
    'r-', 'LineWidth', 2.0, 'DisplayName', '抑制后');

plot([f1/1e6, f1/1e6], [-100, 50], 'k--', 'LineWidth', 0.8);
plot([f2/1e6, f2/1e6], [-100, 50], 'k--', 'LineWidth', 0.8);
text(f1/1e6+0.1, 30, 'f1', 'FontSize', 9, 'Color', 'k');
text(f2/1e6+0.1, 30, 'f2', 'FontSize', 9, 'Color', 'k');

xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
title('抑制前后细节对比（加窗分析）');
xlim([-15, 15]); ylim([-200, 50]); grid on;
legend('Location', 'best', 'FontSize', 9); box on;

%% 8. 输出性能报告
fprintf('\n==================== BPSK信号测试报告 ====================\n');
fprintf('载波频率: %.1f MHz\n', fc/1e6);
fprintf('比特率: %.1f kbps\n', Rb/1e3);
fprintf('滚降系数: %.2f\n', rolloff);
fprintf('带外抑制比: %.1f dB\n', avg_suppression);
fprintf('带内信号质量改善: %.1f dB\n', suppression_results(1));
fprintf('=====================================================\n');

%% 9. 保存结果
save('bpsk_results.mat', 'y_received', 'y_final', 'avg_suppression', ...
    'suppression_results', 'fc', 'Rb', 'fs', 'M');
fprintf('结果已保存至 bpsk_results.mat\n');