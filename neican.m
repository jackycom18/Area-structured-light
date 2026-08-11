% =========================================================================
% 张正友标定法 — 仅求解内参矩阵 K
% 参考文献: Z. Zhang, "A Flexible New Technique for Camera Calibration", PAMI 2000
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
% 棋盘格参数
squareSize = 25;        % 棋盘格方格边长（mm，按实际修改）
numRows    = 8;         % 棋盘格内角点行数
numCols    = 11;        % 棋盘格内角点列数

% 图像路径（支持 jpg/png/bmp）
imageDir  = 'calib_images/';   % 标定图片所在文件夹
imageExt  = '*.jpg';            % 图片格式

%% ======================== 准备世界坐标 ========================
% 棋盘格在世界坐标系 Z=0 平面上
[Xw, Yw] = meshgrid(0:numCols-1, 0:numRows-1);
worldPts = [Xw(:) * squareSize, Yw(:) * squareSize, zeros(numRows*numCols, 1)];

%% ======================== 读取图像并检测角点 ========================
imageFiles = dir(fullfile(imageDir, imageExt));
if isempty(imageFiles)
    error('在 %s 中未找到 %s 文件', imageDir, imageExt);
end

fprintf('共找到 %d 张图片\n', length(imageFiles));

validIdx = [];
allImagePts = {};   % 每张图检测到的角点

for i = 1:length(imageFiles)
    imgPath = fullfile(imageDir, imageFiles(i).name);
    img = imread(imgPath);
    
    % 检测棋盘格角点
    [imagePts, boardSize] = detectCheckerboardPoints(img);
    
    % 验证是否检测到正确数量的角点
    if size(imagePts,1) == numRows * numCols
        validIdx = [validIdx, i];
        allImagePts{end+1} = imagePts;
        fprintf('  [%2d] %s  ✓ 检测到 %d 个角点\n', i, imageFiles(i).name, size(imagePts,1));
    else
        fprintf('  [%2d] %s  ✗ 角点数不足 (检测到 %d, 期望 %d)\n', ...
            i, imageFiles(i).name, size(imagePts,1), numRows*numCols);
    end
end

numViews = length(allImagePts);
fprintf('\n有效视角数: %d\n', numViews);

if numViews < 3
    error('至少需要 3 个有效视角才能求解内参');
end

%% ======================== 步骤1：计算各视角的单应性矩阵 H ========================
fprintf('\n计算单应性矩阵...\n');

H_cell = cell(1, numViews);
worldPts2D = worldPts(:, 1:2);  % 只取 X,Y（Z=0）

for v = 1:numViews
    H_cell{v} = computeHomography(worldPts2D, allImagePts{v});
end

%% ======================== 步骤2：闭式解求解内参矩阵 K ========================
fprintf('求解内参矩阵...\n');

K = solveIntrinsics(H_cell);

fprintf('\n========== 内参矩阵 K ==========\n');
fprintf('  fx = %.4f  (焦距 x 方向，像素)\n', K(1,1));
fprintf('  fy = %.4f  (焦距 y 方向，像素)\n', K(2,2));
fprintf('  cx = %.4f  (主点 x 坐标，像素)\n', K(1,3));
fprintf('  cy = %.4f  (主点 y 坐标，像素)\n', K(2,3));
fprintf('  s  = %.4f  (倾斜因子)\n', K(1,2));
fprintf('\n内参矩阵:\n');
disp(K);

%% ======================== 步骤3：评估重投影误差 ========================
fprintf('计算重投影误差...\n');

totalError = 0;
totalPts = 0;
errors_per_view = zeros(1, numViews);

for v = 1:numViews
    % 用当前单应性（等价于外参已知）计算重投影
    H = H_cell{v};
    projPts = (K * (K \ [worldPts2D, ones(size(worldPts2D,1),1)]'))';  % 近似
    
    % 更精确的重投影：用 H 直接映射
    worldHomog = [worldPts2D, ones(size(worldPts2D,1),1)];
    projHomog = (H * worldHomog')';
    projPts = projHomog(:, 1:2) ./ projHomog(:, 3);
    
    errors = sqrt(sum((allImagePts{v} - projPts).^2, 2));
    errors_per_view(v) = mean(errors);
    totalError = totalError + sum(errors);
    totalPts = totalPts + length(errors);
end

meanError = totalError / totalPts;

fprintf('\n========== 重投影误差 ==========\n');
fprintf('平均误差: %.4f 像素\n', meanError);
fprintf('各视角误差:\n');
for v = 1:numViews
    fprintf('  视角 %2d: %.4f 像素\n', v, errors_per_view(v));
end

%% ======================== 保存内参 ========================
save('camera_intrinsics.mat', 'K');
fprintf('\n内参矩阵 K 已保存至 camera_intrinsics.mat\n');

% =========================================================================
%                         辅助函数
% =========================================================================

%% 函数：计算单应性矩阵 H (DLT 方法)
function H = computeHomography(worldPts, imagePts)
    % 使用直接线性变换 (DLT) 计算单应性矩阵
    % 输入：
    %   worldPts - Nx2 世界坐标（Z=0 平面）
    %   imagePts - Nx2 图像坐标
    % 输出：
    %   H - 3x3 单应性矩阵
    
    N = size(worldPts, 1);
    
    % 构建线性方程组 A*h = 0
    A = zeros(2*N, 9);
    for i = 1:N
        X = worldPts(i, 1);
        Y = worldPts(i, 2);
        u = imagePts(i, 1);
        v = imagePts(i, 2);
        
        A(2*i-1, :) = [X, Y, 1, 0, 0, 0, -u*X, -u*Y, -u];
        A(2*i,   :) = [0, 0, 0, X, Y, 1, -v*X, -v*Y, -v];
    end
    
    % SVD 分解，取最小奇异值对应的向量
    [~, ~, V] = svd(A);
    h = V(:, 9);
    
    % 重塑为 3x3 矩阵
    H = reshape(h, [3, 3])';
    
    % 归一化（使 H(3,3) ≈ 1）
    H = H / H(3,3);
end

%% 函数：闭式解求解内参矩阵 K
function K = solveIntrinsics(H_cell)
    % 使用张正友方法求解内参矩阵
    % 利用旋转矩阵列向量正交约束：h1'*B*h2 = 0, h1'*B*h1 = h2'*B*h2
    % 其中 B = K^(-T) * K^(-1)
    
    numViews = length(H_cell);
    
    % 构建关于 b 的线性方程组 V * b = 0
    % b = [B11, B12, B22, B13, B23, B33]'
    V = zeros(2*numViews, 6);
    
    for v = 1:numViews
        H = H_cell{v};
        h1 = H(:, 1);
        h2 = H(:, 2);
        
        % v12 = [h1*h2 关于 b 的系数]
        v12 = [h1(1)*h2(1), ...
               h1(1)*h2(2) + h1(2)*h2(1), ...
               h1(2)*h2(2), ...
               h1(3)*h2(1) + h1(1)*h2(3), ...
               h1(3)*h2(2) + h1(2)*h2(3), ...
               h1(3)*h2(3)];
        
        % v11 - v22 = [h1*h1 - h2*h2 关于 b 的系数]
        v11_v22 = [h1(1)^2 - h2(1)^2, ...
                   2*(h1(1)*h1(2) - h2(1)*h2(2)), ...
                   h1(2)^2 - h2(2)^2, ...
                   2*(h1(3)*h1(1) - h2(3)*h2(1)), ...
                   2*(h1(3)*h1(2) - h2(3)*h2(2)), ...
                   h1(3)^2 - h2(3)^2];
        
        V(2*v-1, :) = v12;
        V(2*v,   :) = v11_v22;
    end
    
    % SVD 求解 V * b = 0
    [~, ~, Vt] = svd(V);
    b = Vt(:, 6);
    
    % 从 b 提取 B 矩阵元素
    B11 = b(1); B12 = b(2); B22 = b(3);
    B13 = b(4); B23 = b(5); B33 = b(6);
    
    % 从 B 计算内参
    v0 = (B12*B13 - B11*B23) / (B11*B22 - B12^2);
    lambda = B33 - (B13^2 + v0*(B12*B13 - B11*B23)) / B11;
    alpha = sqrt(lambda / B11);
    beta  = sqrt(lambda * B11 / (B11*B22 - B12^2));
    gamma = -B12 * alpha^2 * beta / lambda;
    u0 = gamma * v0 / beta - B13 * alpha^2 / lambda;
    
    % 构建内参矩阵
    K = [alpha, gamma, u0;
         0,     beta,  v0;
         0,     0,     1];
end