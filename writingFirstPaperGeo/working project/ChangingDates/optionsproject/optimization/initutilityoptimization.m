function initutilityoptimization()
%INITMATLAB Make sure tha path is correct
global iuo3
if iuo3==true
else
    addpath '..\utilityoptimization';
end
iuo3 = true;
end

