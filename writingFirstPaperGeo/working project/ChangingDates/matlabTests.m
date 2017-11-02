%
%   Error handling for the MATLAB build
%

try 

	stars = '\n**********************\n%s\n**********************\n\n';
	work = pwd();
		
	cd (work);
	cd 'optionsproject\\utilityoptimization';
	fprintf(stars,pwd());
	run('testutilityoptimization.m');
	
	cd (work);
	cd 'optionsproject\\optimization';
	fprintf(stars,pwd());
	run('testoptimization.m');

	
	fprintf(stars,'All tests passed');
	exit(0);
catch e
    fprintf(stars,'An exception occurred');
    disp(e);
	exit(1);
end	