% Divide images into patches

% Original code: NoRMCorre/construct_grid.m && NoRMCOrre/mat2cell_ov.m
%     https://github.com/flatironinstitute/NoRMCorre

% Constructed by ELi, 20230104
% add sidelobe, 20230207
function [P, PPara] = mat2patch(I, patchN, ovFactor, min_patch_size, sidelobe)
%% check inputs
if nargin < 1 || isempty(I); error('Error: no input'); end
if ndims(I) > 3; error('Error: input size not supported'); end
[d1_sidelobe,d2_sidelobe,d3_sidelobe] = size(I);

if nargin < 2 || isempty(patchN); patchN = [1,1,1]; end
patchN = patchN(:)';
if length(patchN) == 1; patchN = [patchN, patchN, 1]; end
if length(patchN) == 2; patchN = [patchN, 1]; end
patchN = patchN(1:3);

if nargin < 3 || isempty(ovFactor); ovFactor = [0,0,0]; end
ovFactor = ovFactor(:)';
if length(ovFactor) == 1; ovFactor = [ovFactor, ovFactor, 0]; end
if length(ovFactor) == 2; ovFactor = [ovFactor, 0]; end
ovFactor = ovFactor(1:3);

if nargin < 4 || isempty(min_patch_size); min_patch_size = [1,1,1]; end
min_patch_size = min_patch_size(:)';
if length(min_patch_size) == 1; min_patch_size = [min_patch_size, min_patch_size, 1]; end
if length(min_patch_size) == 2; min_patch_size = [min_patch_size, 1]; end
min_patch_size = min_patch_size(1:3);

if nargin < 5 || isempty(sidelobe); sidelobe = [0,0,0]; end
sidelobe = sidelobe(:)';
if length(sidelobe) == 1; sidelobe = [sidelobe,sidelobe,0]; end
if length(sidelobe) == 2; sidelobe = [sidelobe,0]; end
sidelobe = sidelobe(1:3);
I = I(sidelobe(1)+1:d1_sidelobe-sidelobe(1),...
      sidelobe(2)+1:d2_sidelobe-sidelobe(2),...
      sidelobe(3)+1:d3_sidelobe-sidelobe(3));
[d1,d2,d3] = size(I);

%% construct_grid.m
grid_size = ceil([d1,d2,d3] ./ patchN);
mot_uf = [1,1,1];
[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size);

%% mat2cell_ov.m
overlap = ceil(ovFactor .* grid_size);
sz = [d1,d2,d3];
[P,xx_s_ov,xx_f_ov,yy_s_ov,yy_f_ov,zz_s_ov,zz_f_ov] = mat2cell_ov(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz);

%% construct outPara
PPara.patchN = patchN;
PPara.ovFactor = ovFactor;
PPara.min_patch_size = min_patch_size;
PPara.grid_size = grid_size;
PPara.mot_uf = mot_uf;
PPara.overlap = overlap;
PPara.sidelobe = sidelobe;
PPara.xx_s = xx_s;
PPara.xx_f = xx_f;
PPara.yy_s = yy_s;
PPara.yy_f = yy_f;
PPara.zz_s = zz_s;
PPara.zz_f = zz_f;
PPara.xx_s_sidelobe = xx_s + sidelobe(1);
PPara.xx_f_sidelobe = xx_f + sidelobe(1);
PPara.yy_s_sidelobe = yy_s + sidelobe(2);
PPara.yy_f_sidelobe = yy_f + sidelobe(2);
PPara.zz_s_sidelobe = zz_s + sidelobe(3);
PPara.zz_f_sidelobe = zz_f + sidelobe(3);
PPara.xx_s_ov = xx_s_ov;
PPara.xx_f_ov = xx_f_ov;
PPara.yy_s_ov = yy_s_ov;
PPara.yy_f_ov = yy_f_ov;
PPara.zz_s_ov = zz_s_ov;
PPara.zz_f_ov = zz_f_ov;
PPara.xx_s_ov_sidelobe = xx_s_ov + sidelobe(1);
PPara.xx_f_ov_sidelobe = xx_f_ov + sidelobe(1);
PPara.yy_s_ov_sidelobe = yy_s_ov + sidelobe(2);
PPara.yy_f_ov_sidelobe = yy_f_ov + sidelobe(2);
PPara.zz_s_ov_sidelobe = zz_s_ov + sidelobe(3);
PPara.zz_f_ov_sidelobe = zz_f_ov + sidelobe(3);
PPara.sz = sz;
PPara.sz_sidelobe = [d1_sidelobe,d2_sidelobe,d3_sidelobe];
end

%% ------------------------------------------------------------------------
%% NoRMCorre/construct_grid.m
function [xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size)

xx_s = 1:grid_size(1):d1;
yy_s = 1:grid_size(2):d2;
zz_s = 1:grid_size(3):d3;

xx_f = [xx_s(2:end)-1,d1];
yy_f = [yy_s(2:end)-1,d2];
zz_f = [zz_s(2:end)-1,d3];

if xx_f(end)-xx_s(end) + 1 < min_patch_size(1) && length(xx_s) > 1; xx_s(end) = []; xx_f(end-1) = []; end
if yy_f(end)-yy_s(end) + 1 < min_patch_size(2) && length(yy_s) > 1; yy_s(end) = []; yy_f(end-1) = []; end
if zz_f(end)-zz_s(end) + 1 < min_patch_size(3) && length(zz_s) > 1; zz_s(end) = []; zz_f(end-1) = []; end

grid_size_us = floor(grid_size./mot_uf);
if mot_uf(1) > 1
    xx_us = 1:grid_size_us(1):d1;
    xx_uf = [xx_us(2:end)-1,d1];
else
    xx_us = xx_s; xx_uf = xx_f;
end
if mot_uf(2) > 1
    yy_us = 1:grid_size_us(2):d2;
    yy_uf = [yy_us(2:end)-1,d2];
else
    yy_us = yy_s; yy_uf = yy_f;
end
if mot_uf(3) > 1
    zz_us = 1:grid_size_us(3):d3;
    zz_uf = [zz_us(2:end)-1,d3];
else
    zz_us = zz_s; zz_uf = zz_f;
end
end

%% ------------------------------------------------------------------------
%% NoRMCOrre/mat2cell_ov.m
function [I,xx_s_ov,xx_f_ov,yy_s_ov,yy_f_ov,zz_s_ov,zz_f_ov] = mat2cell_ov(X,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)

% converts a matrix into a cell array with overlapping elements
% INPUTS:
% X:            Input matrix
% grid_size:    size of each element without overlap
% overlap:      amount of overlap
% sz:           spatial size of X

% OUTPUT:
% I:            output cell array

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

I = cell(length(xx_s),length(yy_s),length(zz_s));
nd = length(sz);
if nd == 2; sz(3) = 1; end
for i = 1:length(xx_s)
    xx_s_ov(i) = max(xx_s(i)-overlap(1),1);
    xx_f_ov(i) = min(xx_f(i)+overlap(1),sz(1));
end
for j = 1:length(yy_s)
    yy_s_ov(j) = max(yy_s(j)-overlap(2),1);
    yy_f_ov(j) = min(yy_f(j)+overlap(2),sz(2));
end
for k = 1:length(zz_s)
    zz_s_ov(k) = max(zz_s(k)-overlap(3),1);
    zz_f_ov(k) = min(zz_f(k)+overlap(3),sz(3));
end

for i = 1:length(xx_s)
    for j = 1:length(yy_s)
        for k = 1:length(zz_s)            
            if nd == 2
                I{i,j} = X(xx_s_ov(i):xx_f_ov(i),yy_s_ov(j):yy_f_ov(j),:);
            else
                I{i,j,k} = X(xx_s_ov(i):xx_f_ov(i),yy_s_ov(j):yy_f_ov(j),zz_s_ov(k):zz_f_ov(k),:);
            end
        end
    end
end
end