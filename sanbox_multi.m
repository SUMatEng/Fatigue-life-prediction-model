clear all; close all; clc;

tic
%options=optimset('TolFun',10,'TolX',1);%,'Display','iter','PlotFcns',@optimplotfval);
%options=optimset('Display','iter','PlotFcns',@optimplotfval);
options=optimset('TolFun',200,'TolX',1); %TolX is the tolerance on stress, Tol fun is the tol on cycles
[stress,err]=fminsearch(@Nasgro_many_01,400,options);
toc

%%
% clc; close all;
% stress=[];err=[];n=cell(1,1);
% n{1}=@Nasgro_many_01;n{2}=@Nasgro_many_02;n{3}=@Nasgro_many_03;n{4}=@Nasgro_many_04;n{5}=@Nasgro_many_05;
% n{6}=@Nasgro_many_06;n{7}=@Nasgro_many_07;n{8}=@Nasgro_many_08;n{9}=@Nasgro_many_09;n{10}=@Nasgro_many_10;
% n{11}=@Nasgro_many_11;n{12}=@Nasgro_many_12;n{13}=@Nasgro_many_13;n{14}=@Nasgro_many_14;n{15}=@Nasgro_many_15;
% n{16}=@Nasgro_many_16;n{17}=@Nasgro_many_17;n{18}=@Nasgro_many_18;n{19}=@Nasgro_many_19;n{20}=@Nasgro_many_20;
% n{21}=@Nasgro_many_21;n{22}=@Nasgro_many_22;n{23}=@Nasgro_many_23;n{24}=@Nasgro_many_24;n{25}=@Nasgro_many_25;
% n{26}=@Nasgro_many_26;n{27}=@Nasgro_many_27;n{28}=@Nasgro_many_28;n{29}=@Nasgro_many_29;n{30}=@Nasgro_many_30;
% options=optimset('TolFun',200,'TolX',1);
% p = gcp;
% tic
% parfor i=1:30
%     [stress(i),err(i)]=fminsearch(n{i},400,options);
% end
% toc