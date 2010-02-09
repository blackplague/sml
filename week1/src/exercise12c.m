function [ output_args ] = exercise12c( correct_mean, estimated_mean, data_points )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

hold on
scatterplot(data_points, 'x', 'MarkerFaceColor', 'g')

scatterplot(correct_mean, 'o', 'MarkerFaceColor', 'b')

scatterplot(estimated_mean, 'o', 'MarkerFaceColor', 'r')

hold off

end

