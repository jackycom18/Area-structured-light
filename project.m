% =========================================================================
% 多项式拟合法 — 逐像素 列号→深度 关系式
%
% 原理 (达飞鹏,盖绍彦; 刘顺涛 et al.):
%   对每个相机像素 (u,v)，在其 13 个标定板深度下收集 (col, Zc) 对
%   拟合: Zc = a0 + a1*col + a2*col² + a3*col³
%   之后扫物体时: 像素(u,v) → 相位得 col → 代入公式得 Zc
%   再由 Xc = (u-cx)*Zc/fx,  Yc = (v-cy)*Zc/fy 得完整 3D 坐标
%
% 前置: 对 13 个标定板位姿，各运行一次 relationship_mapping.m
%       结果分别放在 data/pos_01/ ~ data/pos_13/
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
nPositions = 13;
dataBase   = 'data/pos_%02d';
polyOrder  = 3;              % 多项式阶数 (3阶推荐, 论文标准)
minSamples = 5;              % 最少需要几个位置的数据才拟合

% 加载相机内参
load(fullfile(sprintf(dataBase,1), 'pixel_to_projector_mapping.mat'), 'K');
[imgH, imgW] = deal(0, 0);  % 稍后从数据确定

%% ======================== 收集数据: 逐像素 (col, Zc) ========================
fprintf('===== 收集 %d 组映射数据 =====\n', nPositions);

% 用 cell 数组暂存每个像素的数据
% cell{u}(v) = [col_list, Zc_list] 的 struct，或用稀疏方式
% 由于像素量大，用稀疏累积+分组方式处理

allU = []; allV = []; allCol = []; allZc = [];  % 所有位置的原始数据

for i = 1:nPositions
    folder = sprintf(dataBase, i);
    file   = fullfile(folder, 'pixel_to_projector_mapping.mat');
    
    if ~exist(file, 'file')
        warning('未找到 %s，跳过', file);
        continue;
    end
    
    data = load(file);
    
    % 标定板局部坐标 → 相机坐标系
    Xc = (data.Rmat * data.worldXYZ' + data.t_c)';
    
    % 记录每个点的 (u, v, col, Zc)
    allU   = [allU;   data.camUV(:,1)];
    allV   = [allV;   data.camUV(:,2)];
    allCol = [allCol; data.projCols];
    allZc  = [allZc;  Xc(:,3)];          % Zc = 相机坐标系深度
    
    fprintf('  位姿 %2d: %5d 个点, Zc∈[%.0f,%.0f]\n', ...
        i, size(data.worldXYZ,1), min(Xc(:,3)), max(Xc(:,3)));
    
end

% 获取图像尺寸（从最后一个有效数据推导）
imgH = max(allV) + 1;
imgW = max(allU) + 1;
fprintf('图像尺寸: %d × %d\n', imgW, imgH);
fprintf('总计: %d 个 (u,v,col,Zc) 数据点\n', length(allU));

%% ======================== 逐像素多项式拟合 ========================
fprintf('\n===== 逐像素多项式拟合 Zc = f(col) =====\n');

% 系数矩阵: coeffs(rows, cols, order+1)
% coeffs(v,u,:) = [a0, a1, a2, a3]
coeffs = NaN(imgH, imgW, polyOrder + 1);
fitR2  = NaN(imgH, imgW);          % R² 拟合优度
fitRMSE = NaN(imgH, imgW);         % RMS 误差
nData   = zeros(imgH, imgW);       % 每像素有效数据点数

% 统计量
nFitted  = 0;  % 成功拟合的像素数
nSkipped = 0;  % 数据不足跳过的像素数

% 逐像素处理
% 先按行分组加速（同一行像素可能来自不同位置的数据）
for v = 1:imgH
    if mod(v, 100) == 0
        fprintf('  处理行 %d / %d ...\n', v, imgH);
    end
    
    rowMask = (allV == v);
    if sum(rowMask) < minSamples
        continue;
    end
    
    rowU   = allU(rowMask);
    rowCol = allCol(rowMask);
    rowZc  = allZc(rowMask);
    
    for u = unique(rowU)'
        colMask = (rowU == u);
        cols = rowCol(colMask);
        zcs  = rowZc(colMask);
        
        nData(v, u) = length(cols);
        
        if nData(v, u) < minSamples
            nSkipped = nSkipped + 1;
            continue;
        end
        
        % 多项式拟合: Zc = a0 + a1*col + a2*col² + a3*col³
        try
            [p, S] = polyfit(cols, zcs, polyOrder);
            coeffs(v, u, :) = p;  % p = [a3, a2, a1, a0] (polyfit 降幂排列)
            
            % 拟合优度 R²
            zc_pred = polyval(p, cols);
            SSres = sum((zcs - zc_pred).^2);
            SStot = sum((zcs - mean(zcs)).^2);
            if SStot > 0
                fitR2(v, u) = 1 - SSres / SStot;
            end
            
            % RMS 误差
            fitRMSE(v, u) = sqrt(SSres / length(cols));
            
            nFitted = nFitted + 1;
        catch
            nSkipped = nSkipped + 1;
        end
    end
end

fprintf('拟合完成: %d 个像素成功, %d 个跳过\n', nFitted, nSkipped);
fprintf('拟合 R² 中位数: %.4f\n', median(fitR2(~isnan(fitR2))));
fprintf('拟合 RMSE 中位数: %.4f mm\n', median(fitRMSE(~isnan(fitRMSE))));

%% ======================== 验证: 用拟合公式反算 3D 坐标 ========================
fprintf('\n===== 拟合验证 =====\n');

% 随机抽取 500 个已拟合像素验证
validPixels = find(~isnan(fitR2));
nSamples = min(500, length(validPixels));
sampleIdx = validPixels(randperm(length(validPixels), nSamples));

sampleErr = zeros(nSamples, 1);

for si = 1:nSamples
    idx = sampleIdx(si);
    [v, u] = ind2sub([imgH, imgW], idx);
    
    p = squeeze(coeffs(v, u, :))';  % [a3, a2, a1, a0]
    
    % 找到该像素的所有实测数据
    mask = (allU == u) & (allV == v);
    cols = allCol(mask);
    zcs_measured  = allZc(mask);
    zcs_predicted = polyval(p, cols);
    
    sampleErr(si) = mean(abs(zcs_measured - zcs_predicted));
end

fprintf('随机 %d 像素验证:\n', nSamples);
fprintf('  平均深度误差: %.4f mm\n', mean(sampleErr));
fprintf('  最大深度误差: %.4f mm\n', max(sampleErr));

%% ======================== 可视化 ========================
figure('Name', '多项式拟合: col→Zc', 'Position', [50, 50, 1500, 700]);

% (1) 每像素数据点数
subplot(2, 3, 1);
imagesc(nData);
axis image; colorbar; colormap(gca, jet);
title(sprintf('每像素有效数据点数 (需≥%d)', minSamples));

% (2) 拟合 R²
subplot(2, 3, 2);
imagesc(fitR2);
axis image; colorbar; caxis([0.9, 1]);
title('拟合 R²');

% (3) 拟合 RMSE
subplot(2, 3, 3);
imagesc(fitRMSE);
axis image; colorbar;
title('拟合 RMSE (mm)');

% (4) 随机5个像素的拟合曲线
subplot(2, 3, 4);
colors = lines(5);
for si = 1:5
    idx = sampleIdx(si);
    [v, u] = ind2sub([imgH, imgW], idx);
    p = squeeze(coeffs(v, u, :))';
    
    mask = (allU == u) & (allV == v);
    cols = allCol(mask);
    zcs  = allZc(mask);
    
    [colsSorted, sortIdx] = sort(cols);
    colFit = linspace(min(cols), max(cols), 100);
    zcFit  = polyval(p, colFit);
    
    plot(cols, zcs, 'o', 'Color', colors(si,:), 'MarkerSize', 3); hold on;
    plot(colFit, zcFit, '-', 'Color', colors(si,:), 'LineWidth', 1.5);
end
xlabel('投影仪列号'); ylabel('Zc (mm)');
title('5个随机像素拟合曲线');
grid on;

% (5) 某个像素的多位置 Zc-col 散点
subplot(2, 3, 5);
% 取一个像素,显示其 13 个位置的数据
demoV = round(imgH/2);
demoU = round(imgW/2);
demoP = squeeze(coeffs(demoV, demoU, :))';
mask = (allU == demoU) & (allV == demoV);
scatter(allCol(mask), allZc(mask), 40, 'b', 'filled'); hold on;
colRange = linspace(min(allCol(mask)), max(allCol(mask)), 100);
plot(colRange, polyval(demoP, colRange), 'r-', 'LineWidth', 2);
xlabel('投影仪列号'); ylabel('Zc (mm)');
title(sprintf('中心像素 (%d,%d): Zc=%.2f+%.4f*col+...', demoU, demoV, ...
    demoP(end), demoP(end-1)));
grid on;

% (6) 深度误差分布
subplot(2, 3, 6);
histogram(fitRMSE(~isnan(fitRMSE)), 50);
xlabel('RMSE (mm)'); ylabel('频数');
title(sprintf('深度拟合 RMSE 分布 (中位数=%.4f mm)', median(fitRMSE(~isnan(fitRMSE)))));

sgtitle('多项式拟合法: 逐像素 列号→深度 关系式');

%% ======================== 保存 ========================
save('per_pixel_depth_model.mat', 'coeffs', 'polyOrder', 'fitR2', 'fitRMSE', ...
    'nData', 'minSamples', 'K', 'imgH', 'imgW');
fprintf('\n逐像素深度模型已保存至 per_pixel_depth_model.mat\n');
fprintf('变量:\n');
fprintf('  coeffs(v,u,:) = [a%d,...,a0]  -- 多项式系数 (降幂)\n', polyOrder);
fprintf('  K             -- 相机内参\n');
fprintf('\n使用时:\n');
fprintf('  1. 对像素 (u,v)，解码得到投影仪列号 col\n');
fprintf('  2. p = squeeze(coeffs(v,u,:));\n');
fprintf('  3. Zc = polyval(p, col);\n');
fprintf('  4. Xc = (u-K(1,3))*Zc/K(1,1);  Yc = (v-K(2,3))*Zc/K(2,2);\n');
