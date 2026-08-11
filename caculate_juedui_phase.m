% =========================================================================
% 格雷码 + 四步相移 → 绝对相位解包裹
% 参考文献: Sansoni et al., Zhang et al., OpenCV structured_light 模块
%
% 原理:
%   格雷码 → 整数周期 k (0~127)
%   四步相移 → 包裹相位 φ ∈ [-π, π]
%   绝对相位: Φ = φ + 2π * k_corrected  (连续, 无锯齿)
%   亚像素列号: col = k_corrected + (φ + π) / (2π)
%
% 前提: 正弦条纹周期 = 1 个格雷码列宽度
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
% 投影仪分辨率（用于计算物理列号, 设为0则只输出归一化相位）
projectorWidth = 0;  % 投影仪水平分辨率（像素），如1920。0=跳过物理列号计算

% 条纹周期占几个格雷码列（通常=1）
fringePeriodsPerGrayCol = 1;

%% ======================== 加载前置结果 ========================
fprintf('加载格雷码结果...\n');
grayData = load('graycode_result.mat');
rayMap    = grayData.rayMap;      % 格雷码列号 0~127
validMask = grayData.validMask;   % 有效像素掩码

fprintf('加载包裹相位结果...\n');
phaseData = load('wrapped_phase_result.mat');
wrappedPhase = phaseData.wrappedPhase;  % [-π, π]
reliability  = phaseData.reliability;

[imgH, imgW] = size(wrappedPhase);
fprintf('图像尺寸: %d × %d\n', imgW, imgH);

%% ======================== 边界矫正 ========================
% 问题: 格雷码判决边界与相位 wrap 边界在实际系统中存在亚像素偏移
% 解决: 利用包裹相位值纠正格雷码列号
% 
%   φ ∈ [0, π)   → 相位在周期"前半段"，若格雷码已跳到下一列，需减1
%   φ ∈ [-π, 0)  → 相位在周期"后半段"，格雷码列号正确
%
% 参考: Sansoni et al. "Three-dimensional vision based on a combination
%        of Gray-code and phase-shift light projection"

fprintf('正在边界矫正...\n');

k = rayMap;  % 初始周期号 = 格雷码列号

% 将相位映射到 [0, 2π) 便于判断
phi_0_2pi = wrappedPhase;
phi_0_2pi(phi_0_2pi < 0) = phi_0_2pi(phi_0_2pi < 0) + 2 * pi;

% 矫正: 相位 < π (前半个周期) → 若格雷码已跳变, 周期号减1
%       避免少数边界像素处格雷码和相位"错位一位"
maxGray = max(rayMap(validMask));
correctionMask = (phi_0_2pi < pi) & (k > 0);
k(correctionMask) = k(correctionMask) - 1;
k = max(min(k, maxGray), 0);

%% ======================== 绝对相位计算 ========================
% Φ = φ + 2π * k  (k 为矫正后的整数周期号)
fprintf('计算绝对相位...\n');
absPhase = wrappedPhase + 2 * pi * k;

% 无效像素置 NaN
absPhase(~validMask) = NaN;
k(~validMask) = NaN;

%% ======================== 亚像素投影仪列号 ========================
% 归一化相位 [0, 1) 映射到一个格雷码列内的亚像素位置
phi_norm = (wrappedPhase + pi) / (2 * pi);  % [-π, π] → [0, 1)

% 亚像素列号 = 整数列号 + 列内小数偏移
absColumn = k + phi_norm;
absColumn(~validMask) = NaN;

%% ======================== 物理投影仪列号（可选） ========================
if projectorWidth > 0
    colsPerGrayCode = projectorWidth / (maxGray + 1);
    fprintf('投影仪宽度: %d px, 每格雷码列: %.1f px\n', projectorWidth, colsPerGrayCode);
    physicalColumn = absColumn * colsPerGrayCode;
    physicalColumn(~validMask) = NaN;
else
    physicalColumn = [];
end

fprintf('解包裹完成。\n');
validPixels = sum(validMask(:));
fprintf('有效像素: %d\n', validPixels);
fprintf('绝对相位范围: [%.4f, %.4f] rad\n', ...
    min(absPhase(validMask)), max(absPhase(validMask)));

%% ======================== 结果可视化 ========================
figure('Name', '绝对相位解包裹', 'Position', [50, 50, 1500, 800]);

% (1) 格雷码列号
subplot(2, 4, 1);
imagesc(rayMap);
axis image; colorbar; colormap(gca, jet);
title(sprintf('格雷码列号 (整数, 0~%d)', maxGray));
caxis([0, maxGray]);

% (2) 包裹相位
subplot(2, 4, 2);
imagesc(wrappedPhase);
axis image; colorbar; colormap(gca, hsv);
title('包裹相位 [-\pi, \pi]');
caxis([-pi, pi]);

% (3) 矫正后周期号 k
subplot(2, 4, 3);
imagesc(k);
axis image; colorbar; colormap(gca, jet);
title('矫正后周期号 k');
caxis([0, maxGray]);

% (4) 绝对相位
subplot(2, 4, 4);
imagesc(absPhase);
axis image; colorbar; colormap(gca, hsv);
title('绝对相位 (连续)');

% (5) 亚像素列号
subplot(2, 4, 5);
imagesc(absColumn);
axis image; colorbar; colormap(gca, jet);
title(sprintf('亚像素列号 [0, %.1f]', maxGray + 1));

% (6) 可靠性图
subplot(2, 4, 6);
imagesc(reliability);
axis image; colorbar;
title('可靠性 (B/A)');

% (7) 中间行剖面对比: 包裹 vs 绝对
subplot(2, 4, 7);
midRow = round(imgH / 2);
plot(1:imgW, wrappedPhase(midRow, :), 'b-', 'LineWidth', 0.5); hold on;
plot(1:imgW, absPhase(midRow, :), 'r-', 'LineWidth', 1.5);
xlabel('列'); ylabel('相位 (rad)');
title(sprintf('第%d行: 包裹(蓝) vs 绝对(红)', midRow));
legend('包裹相位', '绝对相位', 'Location', 'best');
grid on;

% (8) 矫正量统计
subplot(2, 4, 8);
correctionMap = rayMap - k;  % 正值=减了, 0=没变
correctionMap(~validMask) = NaN;
histogram(correctionMap(validMask), 'BinEdges', -1.5:1:1.5);
xlabel('矫正量 (格雷码列 - 矫正后列)');
ylabel('像素数');
title('边界矫正统计');

sgtitle('格雷码 + 四步相移 → 绝对相位解包裹');

%% ======================== 局部放大对比 ========================
figure('Name', '局部放大 — 包裹 vs 绝对', 'Position', [100, 100, 1200, 400]);

% 取中间行, 显示前 200 列看清细节
showRange = 1:min(200, imgW);

subplot(1, 2, 1);
yyaxis left;
plot(showRange, wrappedPhase(midRow, showRange), 'b-', 'LineWidth', 1);
ylabel('包裹相位 (rad)'); ylim([-pi, pi]);
yyaxis right;
stairs(showRange, rayMap(midRow, showRange), 'k-', 'LineWidth', 1);
ylabel('格雷码列号');
xlabel('列'); grid on;
title(sprintf('包裹相位 + 格雷码 (第%d行前200列)', midRow));
legend('包裹相位', '格雷码', 'Location', 'best');

subplot(1, 2, 2);
plot(showRange, absPhase(midRow, showRange), 'r-', 'LineWidth', 1.5);
xlabel('列'); ylabel('绝对相位 (rad)');
title(sprintf('解包裹后绝对相位 (第%d行前200列)', midRow));
grid on;

%% ======================== 保存结果 ========================
if projectorWidth > 0
    save('absolute_phase_result.mat', ...
        'absPhase', 'absColumn', 'k', 'validMask', 'reliability', 'physicalColumn');
else
    save('absolute_phase_result.mat', ...
        'absPhase', 'absColumn', 'k', 'validMask', 'reliability');
end
fprintf('\n结果已保存至 absolute_phase_result.mat\n');
fprintf('主要输出变量:\n');
fprintf('  absPhase  - 绝对相位 (连续, rad)\n');
fprintf('  absColumn - 亚像素投影仪列号\n');
fprintf('  k         - 矫正后整数周期号\n');
