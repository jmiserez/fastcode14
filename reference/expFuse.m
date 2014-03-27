% This script is supposed to reside in ./reference
% 
% Start this script with the following parameters:
% octave expFuse.m <sourcePath> <refPath> <prefix> <M> <Wc> <Ws> <We>

pkg load image

arg_list = argv();
sourcePath = [arg_list{1}];
refPath = [arg_list{2}];
prefix = arg_list{3};
M = arg_list{4};
Wc = arg_list{5};
Ws = arg_list{6};
We = arg_list{7};

outFile = [refPath '/' prefix '-' M '-' Wc '-' Ws '-' We '.tif' ];

sz = size(imread([sourcePath '/' prefix '.0.tif']));
r = floor(sz(1));
c = floor(sz(2));
I = zeros(r,c,3,str2num(M));

for i = 0:str2num(M)
	filename = [sourcePath '/' prefix '.' sprintf("%i", i) '.tif'];
	im = double(imread(filename))/255;
	if (size(im,1) ~= sz(1) || size(im,2) ~= sz(2))
		error('images must all have the same size');
	end

	I(:,:,:,i+1) = im;
end

R = exposure_fusion(I, [str2double(Wc) str2double(Ws) str2double(We)]);
imwrite(R, outFile)
clear all
