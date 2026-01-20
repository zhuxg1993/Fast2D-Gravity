%% ------------------------------------------------------------------------ %%
%% ------------------------------------------------------------------------ %%
% note:Fz, Fx, Fxz, Fzz are kernel function matrices, and Gz, Gx, Gxz, Gzz are the numerical values ??of the gravity components.
% mod is the density difference model
% nodex and nodez are the grid node coordinate files, respectively.
%% set paraments
mod = load('mod.txt');
nodex = load('nodex.txt');
nodez = load('nodez.txt');
%% ------------------------------------------------------------------------ %%
%% ------------------------------------------------------------------------ %%
%% Polygon calculation formula£¬calculate Gkernel (Fx,Fz,Fxz,Fzz) and gravity value (Gx,Gz,Gxz,Gzz)
%% mothed1£ºUnoptimized Algorithm
[Fz,Gz] = cal_Fz_mothed_1(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_mothed_1(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_mothed_1(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_mothed_1(mod,gobs,nodex,nodez);
%% mothed2£ºOnly the calculation formulas for the vertical and horizontal sides were optimized.
[Fz,Gz] = cal_Fz_mothed_2(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_mothed_2(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_mothed_2(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_mothed_2(mod,gobs,nodex,nodez);
%% mothed3£ºThe edge indexing was optimized based on mothed2, meaning each edge is calculated only once.
[Fz,Gz] = cal_Fz_mothed_3(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_mothed_3(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_mothed_3(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_mothed_3(mod,gobs,nodex,nodez);
%% mothed4£ºThe node indexing was optimized based on mothed3, meaning that the arctan and ln terms for each node are calculated only once.
[Fz,Gz] = cal_Fz_mothed_4(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_mothed_4(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_mothed_4(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_mothed_4(mod,gobs,nodex,nodez);
%% ------------------------------------------------------------------------ %%
%% ------------------------------------------------------------------------ %%
%% Rectangle calculation formula£¬calculate Gkernel (Fx,Fz,Fxz,Fzz) and gravity value (Gx,Gz,Gxz,Gzz)
%% mothed1£ºUnoptimized Algorithm
[Fz,Gz] = cal_Fz_rectangle_mothed_1(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_rectangle_mothed_1(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_rectangle_mothed_1(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_rectangle_mothed_1(mod,gobs,nodex,nodez);
%% mothed2£ºOptimization Algorithm
[Fz,Gz] = cal_Fz_rectangle_mothed_2(mod,gobs,nodex,nodez);
[Fx,Gx] = cal_Fx_rectangle_mothed_2(mod,gobs,nodex,nodez);
[Fxz,Gxz] = cal_Fxz_rectangle_mothed_2(mod,gobs,nodex,nodez);
[Fzz,Gzz] = cal_Fzz_rectangle_mothed_2(mod,gobs,nodex,nodez);
%% ------------------------------------------------------------------------ %%
%% ------------------------------------------------------------------------ %%
