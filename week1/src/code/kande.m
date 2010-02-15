function k = kande()

Q8()
k = 0;  
end

function k = Q8()
image = imread('kande1', 'pnm');
testset = image(264:328,150:330,:);
size(testset)
image(testset)
k = 0;
end