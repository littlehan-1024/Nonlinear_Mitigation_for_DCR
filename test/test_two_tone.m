%% 双音信号非线性失真抑制测试
clear; close all; clc;
addpath(genpath(pwd));

%% 1. 基本参数设置
fs = 25e6;          % 采样率
T = 1e-4;           % 信号时长
t = 0:1/fs:T-1/fs;
N = length(t);

% 双音信号参数
f1 = 2.3e6;
f2 = 2.9e6;
A1 = 0.3;
A2 = 0.3;
P_in_dBm = -30;
P_in = 10^((P_in_dBm-30)/10);

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

% 滤波器参数
f_low = 2.2e6;
f_high = 3.0e6;

%% 2. 生成双音信号
fprintf('生成双音信号...\n');
x1 = A1 * exp(1j*2*pi*f1*t);
x2 = A2 * exp(1j*2*pi*f2*t);
x = (x1 + x2) / sqrt(2);

% 显示原始信号频谱
figure(1);
N_fft = min(8192, N);
f_axis = (-N_fft/2:N_fft/2-1)*fs/N_fft/1e6;
plot(f_axis, 20*log10(abs(fftshift(fft(x, N_fft))))+1e-12);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('原始双音信号频谱');
xlim([-15, 15]); ylim([-100, 50]); grid on;
hold on;
plot([f1/1e6, f1/1e6], [-100, 50], 'r--');
plot([f2/1e6, f2/1e6], [-100, 50], 'r--');
legend('信号', 'f1', 'f2');

%% 3. 非线性失真建模
fprintf('模拟非线性失真...\n');
[y_received, y_t, y_IQ] = nonlinear_distortion_model(x, alpha, iq_params, ...
    'SNR_dB', 61, 'P_in_dBm', -30);

% 显示失真后信号
figure(2);
subplot(2,1,1);
plot(f_axis, 20*log10(abs(fftshift(fft(y_received, N_fft))))+1e-12);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('非线性失真后信号频谱');
xlim([-15, 15]); ylim([-120, 50]); grid on;
hold on;

% 标记关键频率
f_imd1 = 2*f1 - f2;
f_imd2 = 2*f2 - f1;
f_mirror1 = -f1;
f_mirror2 = -f2;

plot([f1/1e6, f1/1e6], [-120, 50], 'r--');
plot([f2/1e6, f2/1e6], [-120, 50], 'r--');
plot([f_mirror1/1e6, f_mirror1/1e6], [-120, 50], 'm--');
plot([f_imd1/1e6, f_imd1/1e6], [-120, 50], 'g--');
legend('信号', 'f1/f2', '镜像', 'IMD3');

%% 4. 非线性失真抑制
fprintf('进行非线性失真抑制...\n');
[d, x_hat, y_final, w, e, avg_suppression, suppression_results] = ...
    nonlinear_cancellation(y_received, fs, f_low, f_high, M, mu, alpha_nlms, f1, f2, ...
    'signal_type', 'two_tone', 'N_fft', N_fft);

%% 5. 显示抑制结果
figure(2);
subplot(2,1,2);
Y_final_fft = fftshift(fft(y_final, N_fft));
plot(f_axis, 20*log10(abs(Y_final_fft)+1e-12), 'g', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('抑制后信号频谱');
xlim([-15, 15]); ylim([-120, 50]); grid on;
hold on;
plot([f1/1e6, f1/1e6], [-120, 50], 'r--');
plot([f2/1e6, f2/1e6], [-120, 50], 'r--');
legend('抑制后信号', 'f1', 'f2');

%% 6. 对比图
figure(3);
Y_before_fft = fftshift(fft(y_received, N_fft));
plot(f_axis, 20*log10(abs(Y_before_fft)+1e-12), 'b', 'LineWidth', 1.5);
hold on;
plot(f_axis, 20*log10(abs(Y_final_fft)+1e-12), 'g', 'LineWidth', 1.5);
xlabel('频率 (MHz)'); ylabel('幅度 (dB)'); title('抑制前后频谱对比');
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
plot([f_mirror1/1e6, f_mirror1/1e6], [-100, 50], 'm:', 'LineWidth', 1.0);
plot([f_imd1/1e6, f_imd1/1e6], [-100, 50], 'g:', 'LineWidth', 1.0);

text(f1/1e6+0.1, 30, 'f1', 'FontSize', 9, 'Color', 'k');
text(f2/1e6+0.1, 30, 'f2', 'FontSize', 9, 'Color', 'k');
text(f_mirror1/1e6+0.1, -60, '镜像', 'FontSize', 9, 'Color', 'm');
text(f_imd1/1e6+0.1, -70, 'IMD3', 'FontSize', 9, 'Color', 'g');

xlabel('频率 (MHz)'); ylabel('幅度 (dB)');
title('抑制前后细节对比（加窗分析）');
xlim([-15, 15]); ylim([-200, 50]); grid on;
legend('Location', 'best', 'FontSize', 9); box on;

%% 8. 输出性能报告
fprintf('\n==================== 双音信号测试报告 ====================\n');
fprintf('基于并行自适应结构的非线性失真抑制性能:\n\n');

fprintf('1. 总体性能:\n');
fprintf('   平均非线性失真抑制比: %.1f dB\n', avg_suppression);
fprintf('   与论文结果(32.6 dB)差异: %.1f dB\n\n', avg_suppression - 32.6);

fprintf('2. 各分量抑制性能:\n');
band_names = {'原始信号', '镜像', 'IMD', '三阶谐波'};
for i = 1:length(band_names)
    fprintf('   %s: %.1f dB\n', band_names{i}, suppression_results(i));
end
fprintf('\n');

fprintf('3. 算法配置:\n');
fprintf('   采样率: %.1f MHz\n', fs/1e6);
fprintf('   滤波器长度: %d\n', M);
fprintf('   信号频率: f1=%.1f MHz, f2=%.1f MHz\n', f1/1e6, f2/1e6);
fprintf('   滤波器频带: %.1f-%.1f MHz\n', f_low/1e6, f_high/1e6);
fprintf('=====================================================\n');

%% 9. 保存结果
save('two_tone_results.mat', 'y_received', 'y_final', 'avg_suppression', ...
    'suppression_results', 'f1', 'f2', 'fs', 'M');
fprintf('结果已保存至 two_tone_results.mat\n');