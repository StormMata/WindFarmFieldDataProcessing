function [data] = LoadFullData()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Field Data Power Modeling/functions/figplots');
addpath('/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Field Data Power Modeling/functions');

folderCoord = '/Users/stormmata/Library/Mobile Documents/com~apple~CloudDocs/Courses/Thesis/Field Data Power Modeling/data/';

data = load(strcat(folderCoord, 'AllData_NP_corrected.mat')); 
data = data.AllData;

end