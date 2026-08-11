% =========================================================================
% 四步相移法 — 包裹相位计算
% 参考文献: 四步相移算法标准公式，每次移相 π/2 (90°)
%   I1 = A + B*cos(φ)           (0°)
%   I2 = A + B*cos(φ + π/2)     (90°)
%   I3 = A + B*cos(φ + π)       (180°)
%   I4 = A + B*cos(φ + 3π/2)    (270°)
%   包裹相位: φ = atan2(I4 - I2, I1 - I3)  ∈ [-π, π]
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
imageDir  = 'phase_images/';   % 存放4张相移图片的文件夹
phaseFiles = {'phase_0.png', ...    % 相移 0°   (I1)
              'phase_90.png', ...   % 相移 90°  (I2)
              'phase_180.png', ...  % 相移 180° (I3)
              'phase_270.png'};     % 相移 270° (I4)

%% ======================== 加载图像 ========================
fprintf('正在加载图像...\n');
I = zeros(0, 0, 4);  % 预分配
for k = 1:4
    img = imread(fullfile(imageDir, phaseFiles{k}));
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    if k == 1
        [imgH, imgW] = size(img);
        I = zeros(imgH, imgW, 4);
    end
    I(:, :, k) = double(img);
end
fprintf('图像加载完成。图像尺寸: %d × %d\n', imgW, imgH);

I1 = I(:, :, 1);
I2 = I(:, :, 2);
I3 = I(:, :, 3);
I4 = I(:, :, 4);

%% ======================== 计算包裹相位 ========================
% φ = atan2(I4 - I2, I1 - I3)
% atan2 自动处理四个象限，输出范围 [-π, π]
fprintf('正在计算包裹相位...\n');
wrappedPhase = atan2(I4 - I2, I1 - I3);

%% ======================== 辅助信息 ========================
% 平均光强（背景项 A）
avgIntensity = (I1 + I2 + I3 + I4) / 4;

% 调制度（条纹对比度 B）
% 推导: I4-I2 = 2B*sin(φ), I1-I3 = 2B*cos(φ)
%       B = sqrt((I4-I2)² + (I1-I3)²) / 2
modulation = sqrt((I4 - I2).^2 + (I1 - I3).^2) / 2;

% 可靠性图: 调制度 / 平均光强
% 值越大说明该像素的条纹信号越强，相位越可靠
% 用于后续相位解包裹或点云计算时筛选有效点
reliability = modulation ./ (avgIntensity + 1e-8);

fprintf('相位计算完成。\n');
fprintf('相位范围: [%.4f, %.4f] (弧度)  ≈ [%.1f°, %.1f°]\n', ...
    min(wrappedPhase(:)), max(wrappedPhase(:)), ...
    rad2deg(min(wrappedPhase(:))), rad2deg(max(wrappedPhase(:))));

%% ======================== 结果可视化 ========================
figure('Name', '四步相移法 — 包裹相位', 'Position', [50, 50, 1400, 700]);

% 子图1~4：四张原始相移图
for k = 1:4
    subplot(2, 4, k);
    imshow(uint8(I(:, :, k)));
    title(sprintf('I_%d (相移 %d°)', k, (k-1)*90));
end

% 子图5：包裹相位（伪彩色）
subplot(2, 4, 5);
imagesc(wrappedPhase);
axis image; colorbar; colormap(gca, hsv);
title('包裹相位 \phi [-\pi, \pi]');
caxis([-pi, pi]);

% 子图6：包裹相位（灰度）
subplot(2, 4, 6);
imagesc(wrappedPhase);
axis image; colorbar; colormap(gca, gray);
title('包裹相位 (灰度)');
caxis([-pi, pi]);

% 子图7：调制度
subplot(2, 4, 7);
imagesc(modulation);
axis image; colorbar;
title('调制度 B');

% 子图8：可靠性图
subplot(2, 4, 8);
imagesc(reliability);
axis image; colorbar;
title('可靠性 (B/A)');

sgtitle('四步相移法 — 包裹相位计算结果');

%% ======================== 显示某一行相位剖面 ========================
figure('Name', '相位剖面', 'Position', [100, 100, 900, 400]);

% 取中间行
midRow = round(imgH / 2);

subplot(1, 2, 1);
plot(1:imgW, wrappedPhase(midRow, :), 'b-', 'LineWidth', 1);
xlabel('列 (像素)'); ylabel('相位 (弧度)');
title(sprintf('第 %d 行相位剖面', midRow));
grid on;
ylim([-pi, pi]);

subplot(1, 2, 2);
plot(1:imgW, reliability(midRow, :), 'r-', 'LineWidth', 1);
xlabel('列 (像素)'); ylabel('可靠性 B/A');
title(sprintf('第 %d 行可靠性剖面', midRow));
grid on;

%% ======================== 保存结果 ========================
save('wrapped_phase_result.mat', ...
    'wrappedPhase', 'avgIntensity', 'modulation', 'reliability');
fprintf('\n结果已保存至 wrapped_phase_result.mat\n');
