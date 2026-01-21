function [d, x_hat, y_final, w, e, avg_suppression, suppression_results] = ...
         nonlinear_cancellation(y_received, fs, f_low, f_high, M, mu, alpha_nlms, f1, f2, varargin)
    % 非线性失真抑制算法主函数
    % 输入参数:
    %   y_received: 接收到的含失真信号
    %   fs: 采样率 (Hz)
    %   f_low: 带通滤波器下截止频率 (Hz)
    %   f_high: 带通滤波器上截止频率 (Hz)
    %   M: 滤波器阶数
    %   mu: NLMS步长向量 [1×5]
    %   alpha_nlms: NLMS正则化参数向量 [1×5]
    %   f1, f2: 信号频率 (Hz)，用于性能评估
    %   varargin: 可选参数，包括'signal_type'
    % 输出参数:
    %   d: 失真参考信号
    %   x_hat: 估计的干净信号
    %   y_final: 抑制后的信号
    %   w: 滤波器系数
    %   e: 误差信号
    %   avg_suppression: 平均抑制比 (dB)
    %   suppression_results: 各频带抑制结果
    
    % 参数解析
    p = inputParser;
    addParameter(p, 'signal_type', 'two_tone', @ischar);  % 默认双音信号
    addParameter(p, 'fc', 0, @isnumeric);  % 载波频率（用于调制信号）
    addParameter(p, 'N_fft', 8192, @isnumeric);  % FFT点数
    parse(p, varargin{:});
    params = p.Results;
    
    N = length(y_received);
    N_fft = params.N_fft;
    
    f_axis = (-N_fft/2:N_fft/2-1)*fs/N_fft/1e6;
    
    if strcmp(params.signal_type, 'two_tone')
        f_harmonic_pos = 3 * (f1+f2)/2;
        f_harmonic_neg = -3 * (f1+f2)/2;
        f_imd1 = 2*f1 - f2;
        f_imd2 = 2*f2 - f1;
        f_mirror1 = -f1;
        f_mirror2 = -f2;
    end
    
    %% 滤波器设计
    Fpass = [f_low, f_high];
    
    b_cs = fir1(150, Fpass/(fs/2), 'stop', hanning(151));
    
    b_cp = fir1(150, Fpass/(fs/2), 'bandpass', hanning(151));
    B_cp = hilbert(b_cp);
    
    [b_rs, a_rs] = butter(8, Fpass*2/fs, 'stop');
    
    x_hat = filter(B_cp, 1, y_received);
    
    d = filter(b_cs, 1, y_received);
    
    %% 构建非线性失真项
    term1 = conj(x_hat);
    
    z = x_hat .* abs(x_hat).^2;
    z_rs = filter(b_rs, a_rs, z);
    term2_direct = z_rs;
    term2_conj = conj(z);
    
    x_cube = x_hat.^3;
    x_cube_real_rs = filter(b_rs, a_rs, real(x_cube));
    x_cube_imag_rs = filter(b_rs, a_rs, imag(x_cube));
    term3_real = x_cube_real_rs;
    term3_imag = x_cube_imag_rs;
    
    %% NLMS并行自适应抑制
    xx = [term1; term2_direct; term2_conj; term3_real; term3_imag];
    
    [e, y, w] = NLMS(d, xx, M, alpha_nlms, mu);
    
    if isrow(y)
        y_mitigated = y;
    else 
        y_mitigated = y';
    end
    
    y_final = y_mitigated + x_hat;
    
    %% 性能评估
    Y_BB_fft = fftshift(fft(y_received, N_fft));
    Y_final_fft = fftshift(fft(y_final, N_fft));
    
    if strcmp(params.signal_type, 'two_tone')
        key_bands = {
            '原始信号', [f1-0.1e6, f1+0.1e6; f2-0.1e6, f2+0.1e6];
            '镜像', [f_mirror1-0.1e6, f_mirror1+0.1e6; f_mirror2-0.1e6, f_mirror2+0.1e6];
            'IMD', [f_imd1-0.1e6, f_imd1+0.1e6; f_imd2-0.1e6, f_imd2+0.1e6];
            '三阶谐波', [f_harmonic_pos-0.2e6, f_harmonic_pos+0.2e6];
        };
        
        suppression_results = zeros(length(key_bands), 1);
        
        for i = 1:length(key_bands)
            band = key_bands{i,2};
            P_before = 0; P_after = 0;
            
            if size(band,1) == 1
                idx = find(f_axis >= band(1,1) & f_axis <= band(1,2));
                P_before = sum(abs(Y_BB_fft(idx)).^2);
                P_after = sum(abs(Y_final_fft(idx)).^2);
            else
                for b = 1:size(band,1)
                    idx = find(f_axis >= band(b,1) & f_axis <= band(b,2));
                    P_before = P_before + sum(abs(Y_BB_fft(idx)).^2);
                    P_after = P_after + sum(abs(Y_final_fft(idx)).^2);
                end
            end
            
            if P_after > 0
                suppression_results(i) = 20*log10(P_before/P_after);
            else
                suppression_results(i) = 40;
            end
        end
        
        original_bands = [f1-0.3e6, f1+0.3e6; f2-0.3e6, f2+0.3e6];
        P_distortion_before = 0;
        P_distortion_after = 0;
        
        for k = 1:N_fft
            f_k = f_axis(k);
            in_original = false;
            for b = 1:size(original_bands,1)
                if f_k >= original_bands(b,1) && f_k <= original_bands(b,2)
                    in_original = true;
                    break;
                end
            end
            if ~in_original
                P_distortion_before = P_distortion_before + abs(Y_BB_fft(k))^2;
                P_distortion_after = P_distortion_after + abs(Y_final_fft(k))^2;
            end
        end
        
        if P_distortion_after > 0
            avg_suppression = 20*log10(P_distortion_before/P_distortion_after);
        else
            avg_suppression = Inf;
        end
        
    else
        % 对于调制信号，计算带外抑制比
        f_center = params.fc;
        f_bandwidth = 0.5e6;  % 假设500kHz带宽
        
        idx_inband = find(f_axis >= (f_center - f_bandwidth) & f_axis <= (f_center + f_bandwidth));
        P_inband_before = sum(abs(Y_BB_fft(idx_inband)).^2);
        P_inband_after = sum(abs(Y_final_fft(idx_inband)).^2);
        
        idx_outband = find(f_axis < (f_center - f_bandwidth) | f_axis > (f_center + f_bandwidth));
        P_outband_before = sum(abs(Y_BB_fft(idx_outband)).^2);
        P_outband_after = sum(abs(Y_final_fft(idx_outband)).^2);
        
        if P_outband_after > 0
            avg_suppression = 20*log10(P_outband_before/P_outband_after);
        else
            avg_suppression = Inf;
        end
        
        suppression_results = [P_inband_before/P_inband_after; avg_suppression];
    end
end

%% NLMS算法实现
function [e,yy,w] = NLMS(d, x, M, alpha, mu)
    d_length = length(d);
    if (d_length <= M)
        print('error: 信号长度小于滤波器阶数！');
        return; 
    end
    if (d_length ~= length(x(1,:)))  
        print('error: 输入信号和参考信号长度不同！');
        return; 
    end

    xx1 = zeros(M,1);
    xx2 = zeros(M,1);
    xx3 = zeros(M,1);
    xx4 = zeros(M,1);
    xx5 = zeros(M,1);
    w1 = zeros(M,1);
    w2 = zeros(M,1);
    w3 = zeros(M,1);
    w4 = zeros(M,1);
    w5 = zeros(M,1);
    y = zeros(1, d_length*M);
    yy = zeros(1, d_length);
    e = zeros(1, d_length);
    w = zeros(M, d_length);
    
    for n = 1:d_length
        for i = 1:M
            xx1 = [xx1(2:M,1);x(1,n)];
            xx2 = [xx2(2:M,1);x(2,n)];
            xx3 = [xx3(2:M,1);x(3,n)];
            xx4 = [xx4(2:M,1);x(4,n)];
            xx5 = [xx5(2:M,1);x(5,n)];
            xx = [xx1;xx2;xx3;xx4;xx5];
            ww = [w1;w2;w3;w4;w5];
            nromed1 = norm(xx1);
            nromed2 = norm(xx2);
            nromed3 = norm(xx3);
            nromed4 = norm(xx4);
            nromed5 = norm(xx5);
            mu_eff1 = mu(1) / (alpha(1) + nromed1);
            mu_eff2 = mu(2) / (alpha(2) + nromed2);
            mu_eff3 = mu(3) / (alpha(3) + nromed3);
            mu_eff4 = mu(4) / (alpha(4) + nromed4);
            mu_eff5 = mu(5) / (alpha(5) + nromed5);
            mu_eff = [mu_eff1;mu_eff2;mu_eff3;mu_eff4;mu_eff5];
            for k = 1:M
                y(k) = ww(k)' * xx(k);
                yy(n) = yy(n) + y(k);
                e(n) = d(n) - yy(n);
            end

            ww(i) = w1(i) + mu_eff(i) * e(n) * xx(i);
            w(i,n) = ww(i);
        end
    end
end