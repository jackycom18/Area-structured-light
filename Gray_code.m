% =========================================================================
% 格雷码解码 — 相机像素 → 投影仪光线列号
% 参考文献: Inokuchi et al. (1984), OpenCV structured_light 模块
% =========================================================================
clear; close all; clc;

%% ======================== 参数设置 ========================
% 图像文件夹路径
imageDir   = 'graycode_images/';   % 存放9张图像的文件夹
whiteFile  = 'white.png';          % 全白投影图像文件名
blackFile  = 'black.png';          % 全黑投影图像文件名
grayFiles  = {'pattern_01.png', ...   % 7张格雷码图案，第1张对应最高位(bit 6)
              'pattern_02.png', ...   % 第2张对应 bit 5
              'pattern_03.png', ...   % ...
              'pattern_04.png', ...
              'pattern_05.png', ...
              'pattern_06.png', ...
              'pattern_07.png'};      % 第7张对应最低位(bit 0)

% 阴影判定阈值：白图与黑图灰度差低于此值的像素视为无效（阴影/无信号）
% 论文建议取 5~15（0-255 灰度范围），按实际图像对比度调整
shadowThreshold = 10;

% 格雷码位数
nBits = 7;  % 2^7 = 128 列

%% ======================== 加载图像 ========================
fprintf('正在加载图像...\n');
whiteImg = imread(fullfile(imageDir, whiteFile));
blackImg = imread(fullfile(imageDir, blackFile));

% 统一转为灰度图（若输入为RGB）
if size(whiteImg, 3) == 3
    whiteImg = rgb2gray(whiteImg);
end
if size(blackImg, 3) == 3
    blackImg = rgb2gray(blackImg);
end

% 转换为 double 便于计算
whiteImg = double(whiteImg);
blackImg = double(blackImg);

[imgH, imgW] = size(whiteImg);

% 加载7张格雷码图案
grayPatterns = zeros(imgH, imgW, nBits);
for i = 1:nBits
    img = imread(fullfile(imageDir, grayFiles{i}));
    if size(img, 3) == 3
        img = rgb2gray(img);
    end
    grayPatterns(:, :, i) = double(img);
end
fprintf('图像加载完成。图像尺寸: %d × %d\n', imgW, imgH);

%% ======================== 逐像素阈值计算 ========================
% threshold(x,y) = (I_white(x,y) + I_black(x,y)) / 2
fprintf('正在计算逐像素阈值...\n');
threshold = (whiteImg + blackImg) / 2;

%% ======================== 有效像素掩码 ========================
% 白-黑 差值太小说明该像素没有收到投影光（阴影、遮挡或超出投影范围）
diffImg = whiteImg - blackImg;
validMask = diffImg > shadowThreshold;
fprintf('有效像素占比: %.1f%%\n', 100 * sum(validMask(:)) / numel(validMask));

%% ======================== 二值化 + 格雷码解码 ========================
% 对每张格雷码图：I > threshold → 1, 否则 → 0
% 得到每个像素的 nBits 位格雷码: [G_{nBits-1}, ..., G_0]
%   格雷码 → 二进制: B_{nBits-1} = G_{nBits-1}
%                     B_i = G_i XOR B_{i+1}  (i 从高位-1 递减到 0)
%   二进制 → 十进制 → 投影仪列号

fprintf('正在解码...\n');

% 初始化格雷码矩阵和射线图
grayCode = false(imgH, imgW, nBits);   % 逻辑型，节省内存
rayMap = NaN(imgH, imgW);               % 输出：投影仪列号

for i = 1:nBits
    grayCode(:, :, i) = grayPatterns(:, :, i) > threshold;
end

% 格雷码 → 二进制 → 十进制（向量化运算，比逐像素循环快）
binaryBits = false(imgH, imgW, nBits);
binaryBits(:, :, nBits) = grayCode(:, :, nBits);  % 最高位不变

for i = nBits-1 : -1 : 1
    % B_i = G_i XOR B_{i+1}
    binaryBits(:, :, i) = xor(grayCode(:, :, i), binaryBits(:, :, i+1));
end

% 二进制位 → 十进制列号
% 权重: bit_i 对应 2^(nBits-i)，即 bit1(最高位)=2^6=64, ..., bit7(最低位)=2^0=1
rayMap = zeros(imgH, imgW);
for i = 1:nBits
    rayMap = rayMap + double(binaryBits(:, :, i)) * 2^(nBits - i);
end

% 无效像素设为 NaN
rayMap(~validMask) = NaN;

fprintf('解码完成。\n');

%% ======================== 结果可视化 ========================
figure('Name', '格雷码解码结果', 'Position', [100, 100, 1200, 500]);

% 子图1：全白图像
subplot(2, 4, 1);
imshow(uint8(whiteImg));
title('全白投影图像');

% 子图2：全黑图像
subplot(2, 4, 2);
imshow(uint8(blackImg));
title('全黑投影图像');

% 子图3：阈值图
subplot(2, 4, 3);
imshow(uint8(threshold));
title('逐像素阈值 (White+Black)/2');

% 子图4：有效像素掩码
subplot(2, 4, 4);
imshow(validMask);
title(sprintf('有效像素 (白-黑 > %d)', shadowThreshold));

% 子图5~7：3张代表性格雷码二值化结果
sampleBits = [1, 4, 7];  % 展示 bit6(最粗), bit3, bit0(最细)
for idx = 1:3
    subplot(2, 4, 4 + idx);
    imshow(grayCode(:, :, sampleBits(idx)));
    title(sprintf('格雷码位 %d (bit %d)', sampleBits(idx), nBits - sampleBits(idx)));
end

% 子图8：解码结果（伪彩色显示列号）
subplot(2, 4, 8);
% 有效区域用伪彩色，无效区域显示为黑色
rayMapDisplay = rayMap;
rayMapDisplay(isnan(rayMapDisplay)) = -1;  % 临时标记
imagesc(rayMapDisplay);
colormap(gca, [0 0 0; jet(128)]);  % 黑色=无效, 彩色=有效列号
caxis([0 127]);
colorbar;
title(sprintf('投影仪列号 (0~%d)', 2^nBits - 1));

sgtitle('Gray Code 解码结果');

%% ======================== 输出统计 ========================
fprintf('\n============ 解码统计 ============\n');
fprintf('图像尺寸: %d × %d\n', imgW, imgH);
fprintf('格雷码位数: %d (编码 %d 列)\n', nBits, 2^nBits);
fprintf('有效像素数: %d / %d (%.1f%%)\n', ...
    sum(validMask(:)), numel(validMask), 100*sum(validMask(:))/numel(validMask));
fprintf('列号范围: %d ~ %d\n', min(rayMap(validMask)), max(rayMap(validMask)));

% 显示解码结果（数值矩阵，仅前20行20列作为示例）
fprintf('\n投影仪列号矩阵示例 (前20×20像素):\n');
disp(rayMap(1:min(20, imgH), 1:min(20, imgW)));

%% ======================== 保存结果 ========================
save('graycode_result.mat', 'rayMap', 'validMask', 'threshold', 'grayCode');
fprintf('\n结果已保存至 graycode_result.mat\n');
