# 宽带射频直采接收机非线性失真抑制系统

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-orange)](https://www.mathworks.com/products/matlab.html)

## 📖 概述

本开源项目实现了一个基于并行自适应滤波的宽带直接变频接收机（DCR）非线性失真抑制算法。该算法能够有效抑制通信系统中由射频非线性、I/Q不平衡和基带非线性引起的各种失真分量，包括互调失真（IMD）、镜像干扰和谐波失真。

## ✨ 特性

- **并行自适应结构**：采用5个并行的NLMS滤波器同时处理多种非线性失真机制
- **多信号类型支持**：支持双音信号、BPSK信号和π/4-DQPSK信号
- **模块化设计**：算法核心与测试脚本分离，便于集成和扩展
- **全面性能评估**：提供频谱分析、抑制比计算、EVM评估等多种性能指标
- **学术研究友好**：完整复现论文Fig.7(a)结果，便于算法对比研究

## 🏗️ 算法结构

### 核心思想
非线性失真抑制通过构建失真模型并进行自适应抵消实现：

```
通过并行自适应滤波器对5项参数进行估计，恢复原始信号
```

### 三级非线性模型
1. **射频非线性**：三阶互调产物生成
2. **I/Q不平衡**：增益和相位不平衡导致的镜像干扰
3. **基带非线性**：基带电路的三阶非线性

### 并行自适应结构
```
输入信号 → 信号分解 → 5个并行失真项 → NLMS自适应滤波 → 信号重建
```
5个并行失真项：
1. **Term1**: 一阶项（处理I/Q不平衡）
2. **Term2_direct**: 三阶非线性直接项
3. **Term2_conj**: 三阶非线性共轭项  
4. **Term3_real**: 立方项实部
5. **Term3_imag**: 立方项虚部

## 📊 适配范围

### 支持的信号类型
- ✅ 双音测试信号（用于算法验证）
- ✅ BPSK调制信号
- ✅ π/4-DQPSK调制信号_work in progress
- 🔄 可扩展支持其他数字调制格式

### 系统参数范围
- **采样率**: 10-100 MHz
- **信号带宽**: 最高5 MHz
- **非线性抑制比**: 最高35+ dB
- **I/Q不平衡**: 增益不平衡≤3dB，相位不平衡≤10°

### 应用场景
- 宽带软件定义无线电（SDR）
- 5G/Wi-Fi接收机非线性补偿
- 卫星通信系统
- 测试与测量设备
- 学术研究和算法开发

## 🚀 使用方法

### 环境要求
- MATLAB R2020b或更高版本
- Signal Processing Toolbox

### 快速开始

1. **克隆仓库**
```bash
git clone https://github.com/yourusername/dcr-nonlinear-cancellation.git
cd dcr-nonlinear-cancellation
```

2. **运行双音信号测试**
```matlab
% 测试基本功能
test_two_tone;
```

3. **运行调制信号测试**
```matlab
% 测试BPSK信号
test_bpsk;

% 测试π/4-DQPSK信号  
test_pi4_dqpsk_wip;
```

### 自定义使用

#### 基本使用流程
```matlab
% 1. 准备输入信号
fs = 25e6;  % 采样率
t = 0:1/fs:1e-3;
f1 = 2.3e6; f2 = 2.9e6;
x = 0.3*exp(1j*2*pi*f1*t) + 0.3*exp(1j*2*pi*f2*t);

% 2. 添加非线性失真
alpha = struct('a1', 5.62, 'a2', -(84351+1j*74391), 'a3', 3.16, 'a4', -1588.7);
iq_params = struct('gm', 0.99, 'phi_m', deg2rad(3.6));
y_received = nonlinear_distortion_model(x, fs, alpha, iq_params);

% 3. 应用失真抑制算法
f_low = 2.2e6; f_high = 3.0e6;
M = 5; mu = [1, 1, 0.01, 1, 1]; alpha_nlms = [1e-9, 1e-8, 1e-4, 1e-9, 1e-8];

[d, x_hat, y_final, w, e, avg_suppression] = ...
    nonlinear_cancellation(y_received, fs, f_low, f_high, M, mu, alpha_nlms, f1, f2);
```

#### 参数说明

**非线性失真模型参数**：
```matlab
alpha.a1 = 5.62;      % 线性增益
alpha.a2 = - (84351 + 1j*74391);  % 三阶非线性系数
alpha.a3 = 3.16;      % 基带线性增益
alpha.a4 = -1588.7;   % 基带三阶非线性系数

iq_params.gm = 0.99;           % I/Q增益不平衡
iq_params.phi_m = deg2rad(3.6); % I/Q相位不平衡
```

**NLMS算法参数**：
```matlab
M = 5;  % 滤波器阶数
mu = [1, 1, 0.01, 1, 1];  % 各分支步长
alpha_nlms = [1e-9, 1e-8, 1e-4, 1e-9, 1e-8];  % 正则化参数
```

### 性能评估

每个测试脚本都会生成：
1. **频谱对比图**：抑制前后频谱对比
2. **细节分析图**：加窗FFT细节对比
3. **性能报告**：各频带抑制比和总体性能
4. **数据文件**：保存关键结果到.mat文件

## 📁 文件结构

```
dcr-nonlinear-cancellation/
├── README.md                  # 本文档
├── src/
│   ├── nonlinear_cancellation.m   # 主抑制算法函数
│   ├── nonlinear_distortion_model.m # 非线性失真模型
├── test/                # 测试结果
│   ├── test_two_tone.m           # 双音信号测试脚本
│   ├── test_bpsk.m              # BPSK信号测试脚本
│   ├── test_pi4_dqpsk.m         # π/4-DQPSK信号测试脚本
└── dco/                # 参考文献
```

## 📈 预期性能

| 信号类型 | 平均抑制比 | 带内EVM改善 | 计算复杂度 |
|---------|-----------|-------------|-----------|
| 双音信号 | 32-36 dB | - | O(N×M) |
| BPSK信号 | 25-30 dB | 15-20% | O(N×M) |
| π/4-DQPSK | 22-28 dB | 10-15% | O(N×M) |

*注：实际性能取决于具体系统参数*

## 🔧 扩展与定制

### 添加新的信号类型
1. 创建新的测试脚本（如`test_qpsk.m`）
2. 遵循`nonlinear_distortion_model`接口
3. 调用`nonlinear_cancellation`函数

### 调整算法参数
- 修改滤波器阶数`M`以平衡性能与复杂度
- 调整步长`mu`以优化收敛速度
- 更改滤波器频带`[f_low, f_high]`适应不同信号

### 集成到现有系统
```matlab
% 在现有接收机处理链中集成
function processed_signal = receiver_chain(rx_signal, fs)
    % 其他处理步骤...
    
    % 非线性失真抑制
    [~, ~, processed_signal] = nonlinear_cancellation(rx_signal, fs, ...);
    
    % 后续处理步骤...
end
```

## 📚 参考文献

1. [原始论文] "Wideband DCR Nonlinearity Cancellation Using Parallel Adaptive Structure"
2. J. Tsimbinos, "Compensation of Nonlinear Distortion in OFDM Systems"
3. S. A. Bassam et al., "I/Q Imbalance Compensation in Wideband Receivers"

## 📄 许可证

本项目采用 **GNU General Public License v3.0** 许可证。

```
Copyright (C) 2024 Your Name

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
```

## 🤝 贡献

欢迎提交Issue和Pull Request！

1. Fork本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 打开Pull Request

## 📧 联系

如有问题或建议，请通过以下方式联系：

- 作者: daihong.han@outlook.com
- GitHub Issues

**⚠️ 免责声明**: 本软件按"原样"提供，不附带任何明示或暗示的保证。使用者需自行承担风险。
