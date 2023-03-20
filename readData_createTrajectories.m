clear all
format long e

% basePath = uigetdir;
basePath = '/Users/matthiasmonden/Downloads/deXtreme_pilot_001/001-noBX-Alch/6_Alchemist_20220125165037';
% save_dir = uigetdir;
save_dir = '/Users/matthiasmonden/Downloads/deXtreme_pilot_001/saves';

atEnd = false;
filenamesDone = {};
index = 3;
figNumber = 1;

rotationMatrix_toYZPlane = zeros(3, 3);
rotationMatrix = zeros(3, 3);
xPoints = {};
yPoints = {};
zPoints = {};


while ~atEnd
	fileList = dir(basePath);

	if length(fileList) < index
		atEnd = true;
		continue
	end

	file = fileList(index).name;

	%	Check if file is already done OR if it is the baseline
	[filepath, name, ext] = fileparts(file);

	if ext ~= ".csv"
		index = index + 1;
		continue
	end

	[x, y, z, orgX, orgY, orgZ] = getdata(append(basePath, "/", file));

	index = index + 1;
end

fig_dist = plotPosition(x, y, z, orgX, orgY, orgZ);
saveas(fig_dist, fullfile(save_dir, file(1:length(file)-4)), 'fig')

function [xData, yData, zData, orgX, orgY, orgZ] = getdata(csvPath)
	opts = detectImportOptions(csvPath);
	M = readtable(csvPath, opts);

	% Calculate size of the imported data table.
	[n, m] = size(M);

	% Extract data from table to array/matrix.
	dataMatrix = M{5:n-2, 1:m-1};

	startX = M{1, 2};
	startY = M{1, 3};
	startZ = M{1, 4};

	endX = M{3, 2};
	endY = M{3, 3};
	endZ = M{3, 4};

	orgX = [startX, endX];
	orgY = [startY, endY];
	orgZ = [startZ, endZ];

	endpoint = [endX; endY; startZ];
	proj = ([norm([endX, endY, endZ]); 0; 0]' * endpoint/norm(endpoint)) * endpoint/norm(endpoint);
	angle_toXYPlane = acos(sqrt((endX-startX)^2 + (endY-startY)^2 + (startZ-startZ)^2) / sqrt((endX-startX)^2 + (endY-startY)^2 + (endZ-startZ)^2));
	rotationMatrix_toYZPlane(1, 1) = cos(angle_toXYPlane);
	rotationMatrix_toYZPlane(1, 3) = -sin(angle_toXYPlane);
	rotationMatrix_toYZPlane(2, 2) = 1;
	rotationMatrix_toYZPlane(3, 1) = sin(angle_toXYPlane);
	rotationMatrix_toYZPlane(3, 3) = cos(angle_toXYPlane);

	proj_end = (rotationMatrix_toYZPlane * [endX; endY; endZ]) - [startX; startY; startZ]

	endpoint = [endX; startY; startZ];
	proj = (proj_end' * endpoint/norm(endpoint)) * endpoint/norm(endpoint);
	angle = acos(sqrt((proj(1)-startX)^2 + (proj(2)-startY)^2 + (proj(3)-startZ)^2) / sqrt((proj_end(1)-startX)^2 + (proj_end(2)-startY)^2 + (proj_end(3)-startZ)^2));
	rotationMatrix(3, 3) = 1;
	rotationMatrix(1, 1) = cos(angle);
	rotationMatrix(1, 2) = -sin(angle);
	rotationMatrix(2, 1) = sin(angle);
	rotationMatrix(2, 2) = cos(angle);

	finalvector = rotationMatrix * proj_end;

	% Put data in column vectors for later use.
	timeData = dataMatrix(:, 1);

	xData = dataMatrix(:, 2);
	yData = dataMatrix(:, 3);
	zData = dataMatrix(:, 4);
end

function [fig] = plotPosition(x, y, z, orgX, orgY, orgZ)
	grid on
	hold on
	fig = plot3(x, y, z, '-o', 'Color', 'b');
	plot3(orgX, orgY, orgZ, '-o', 'Color', 'm');
	hold off
	title('Patients trajectory in space -- 3D')
	xlabel('x')
	ylabel('y')
	zlabel('z')
	legend('Patients averaged trajectory', 'Correct trajectory')
end