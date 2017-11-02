function initMosek()
global initMosekComplete;
if initMosekComplete==true;
else
    addpath 'c:\Program Files\mosek\8\toolbox\r2014a'
end
initMosekComplete = true;
end

