% Function to calculate logistic function value at a given x
function y = logistic_function(x, design_pc, p)
    y = 1 ./ (1 + exp(-design_pc * p));
end