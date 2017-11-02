function saveFigure( num, file, size )
%SAVEFIGURE Summary of this function goes here
%   Detailed explanation goes here
f = figure(num);
set( f, 'PaperSize', size);
set( f, 'PaperPosition', horzcat([0 0], size));
saveas(f, file );

end

