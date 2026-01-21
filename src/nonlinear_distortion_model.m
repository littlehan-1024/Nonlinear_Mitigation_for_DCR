function [y_received, y_t, y_IQ] = nonlinear_distortion_model(x, alpha_params, iq_params, varargin)
    % 非线性失真模型
    % 输入参数:
    %   x: 输入信号
    %   alpha_params: 非线性参数结构体 [a1, a2, a3, a4]
    %   iq_params: I/Q不平衡参数 [gm, phi_m]
    %   varargin: 可选参数
    % 输出参数:
    %   y_received: 接收信号（含所有失真）
    %   y_t: 射频非线性输出
    %   y_IQ: I/Q不平衡输出
    
    % 参数解析
    p = inputParser;
    addParameter(p, 'SNR_dB', 61, @isnumeric);
    addParameter(p, 'P_in_dBm', -30, @isnumeric);
    parse(p, varargin{:});
    params = p.Results;
    
    %% 1. 添加噪声
    if params.SNR_dB > 0
        signal_power = mean(abs(x).^2);
        noise_power = signal_power / (10^(params.SNR_dB/10));
        noise = sqrt(noise_power/2) * (randn(size(x)) + 1j*randn(size(x)));
        x_with_noise = x + noise;
    else
        x_with_noise = x;
    end
    
    %% 2. 功率归一化
    if params.P_in_dBm ~= 0
        P_in = 10^((params.P_in_dBm-30)/10);
        x_power = mean(abs(x_with_noise).^2);
        x_normalized = x_with_noise * sqrt(P_in/x_power);
    else
        x_normalized = x_with_noise;
    end
    
    %% 3. 射频非线性
    x_t = x_normalized;
    y_t = alpha_params.a1 * x_t + 3 * alpha_params.a2 * ((abs(x_t) .^ 2) .* x_t);
    
    %% 4. I/Q不平衡
    gm = iq_params.gm;
    phi_m = iq_params.phi_m;
    k1 = (1 + gm*exp(-1j*phi_m))/2;
    k2 = (1 - gm*exp(1j*phi_m))/2;
    y_IQ = k1 * y_t + k2 * conj(y_t);
    
    %% 5. 基带非线性
    y_I = real(y_IQ);
    y_Q = imag(y_IQ);
    y_I_BB = alpha_params.a3 * y_I + alpha_params.a4 * (y_I.^3);
    y_Q_BB = alpha_params.a3 * y_Q + alpha_params.a4 * (y_Q.^3);
    y_received = y_I_BB + 1j*y_Q_BB;
end