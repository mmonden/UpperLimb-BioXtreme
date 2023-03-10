clear all
format long e

basePath = uigetdir
%basePath = '/Users/matthiasmonden/Downloads/marketstand001/6_Market Stand_20220125164648_baseline';

% dirList = dir(basePath);
% direc = dirList(4).name	% Start vanaf 3 (. en .. zijn ook dirs)

atEnd = false;
filenamesDone = [];
index = 3;
figNumber = 1;

delete(fullfile(basePath, "*.fig"))
delete(fullfile(basePath, "*.xlsx"))

targetReached_arr = [];
reactionTime_arr = [];
initialDirectionAngle_arr = [];
initialDistanceRatio_arr = [];
initialSpeedRatio_arr = [];
speedMaximaCount_arr = [];
minMaxSpeed_arr = [];
momementTime_arr = [];
pathLengthRatio_arr = [];
maxSpeed_arr = [];
avgDeviation_arr = [];

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

	[time, x, y, z, orgX, orgY, orgZ, endX, endY, endZ, supination, flexion, abduction, deviation, averageDeviation] = getdata(append(basePath, "/", file));

	%	Here the NT (No Trial) parameter is calculated
	targetReached = hasReachedTarget(x, y, z, endX, endY, endZ);

	%	Here we calculate the velocities and accelerations TM (Total Movement) parameters
	[vx, vy, vz] = calculateAvgSpeed(time, x, y, z);
	[ax, ay, az] = calculateAvgAccel(time, vx, vy, vz);

	v = sqrt(vx.^2 + vy.^2 + vz.^2);
	a = sqrt(ax.^2 + ay.^2 + az.^2);

	maxSpeed = max(v) * 100;
	% pathLengthRatio is calculated at the FM parameters
	movementTime = time(length(time)) / 1000;

	fig_dist = plotPosition(time, x, y, z, orgX, orgY, orgZ, figNumber);
	saveas(fig_dist, fullfile(basePath, file(1:length(file)-4)), 'fig')

	%	Here we calculate the FM (First Movement) parameters
	[initialDirectionAngle, initialDistanceRatio, initialSpeedRatio, totalDistance] = getFMParameters(time, x, y, z, orgX, orgY, orgZ, v);
	initialSpeedRatio = initialSpeedRatio * 100;
	pathLengthRatio = totalDistance/sqrt((orgX(2) - orgX(1))^2 + (orgY(2) - orgY(1))^2 + (orgZ(2) - orgZ(1))^2);

	%	Now we calculate the CM (Corrective Movement) parameters
	speedMaximaCount = getSpeedMaxCnt(a);
	sp = getMinMaxSpeed(time, v);

	if sp ~= -1
		minMaxSpeed = sum(sp)/length(sp);
	else
		minMaxSpeed = -1;	% If minMaxSpeed == -1; this means that there was no minimum found after a maximum in the "CM test".
	end

	%	Here we calculate the reaction time (VR - Visual Reaction)
	reactionTime = getReactionTime(time);

	targetReached_arr = [targetReached_arr targetReached];
	reactionTime_arr = [reactionTime_arr reactionTime];
	initialDirectionAngle_arr = [initialDirectionAngle_arr initialDirectionAngle];
	initialDistanceRatio_arr = [initialDistanceRatio_arr initialDistanceRatio];
	initialSpeedRatio_arr = [initialSpeedRatio_arr initialSpeedRatio];
	speedMaximaCount_arr = [speedMaximaCount_arr speedMaximaCount];
	minMaxSpeed_arr = [minMaxSpeed_arr minMaxSpeed];
	momementTime_arr = [momementTime_arr movementTime];
	pathLengthRatio_arr = [pathLengthRatio_arr pathLengthRatio];
	maxSpeed_arr = [maxSpeed_arr maxSpeed];
	avgDeviation_arr = [avgDeviation_arr averageDeviation];

	filenamesDone = [filenamesDone; file];
	index = index + 1;
	figNumber = figNumber + 1;
end

[path, filename, ext] = fileparts(basePath);
ext = ".xlsx";
writematrix(targetReached_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "A2");
writematrix(reactionTime_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "B2");
writematrix(initialDirectionAngle_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "C2");
writematrix(initialDistanceRatio_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "D2");
writematrix(initialSpeedRatio_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "E2");
writematrix(speedMaximaCount_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "F2");
writematrix(minMaxSpeed_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "G2");
writematrix(momementTime_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "H2");
writematrix(pathLengthRatio_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "I2");
writematrix(maxSpeed_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "J2");
writematrix(avgDeviation_arr', fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "K2");

writematrix("targetReached", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "A1");
writematrix("reactionTime", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "B1");
writematrix("initialDirectionAngle", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "C1");
writematrix("initialDistanceRatio", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "D1");
writematrix("initialSpeedRatio", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "E1");
writematrix("speedMaximaCount", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "F1");
writematrix("minMaxSpeed", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "G1");
writematrix("momementTime", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "H1");
writematrix("pathLengthRatio", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "I1");
writematrix("maxSpeed", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "J1");
writematrix("avgDeviation", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', "K1");

avg = length(reactionTime_arr) + 2;

index = 1;
while index < length(minMaxSpeed_arr)
	if minMaxSpeed_arr(index) == -1
		minMaxSpeed_arr = [minMaxSpeed_arr(1:index-1), minMaxSpeed_arr(index+1:length(minMaxSpeed_arr))];
	else
		index = index + 1;
	end
end

writematrix("Average", fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range', strcat("L", string(avg)));
writematrix(sum(targetReached_arr)/length(targetReached_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("A", string(avg)));
writematrix(sum(reactionTime_arr)/length(reactionTime_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("B", string(avg)));
writematrix(sum(initialDirectionAngle_arr)/length(initialDirectionAngle_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("C", string(avg)));
writematrix(sum(initialDistanceRatio_arr)/length(initialDistanceRatio_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("D", string(avg)));
writematrix(sum(initialSpeedRatio_arr)/length(initialSpeedRatio_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("E", string(avg)));
writematrix(sum(speedMaximaCount_arr)/length(speedMaximaCount_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("F", string(avg)));
writematrix(sum(minMaxSpeed_arr)/length(minMaxSpeed_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("G", string(avg)));
writematrix(sum(momementTime_arr)/length(momementTime_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("H", string(avg)));
writematrix(sum(pathLengthRatio_arr)/length(pathLengthRatio_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("I", string(avg)));
writematrix(sum(maxSpeed_arr)/length(maxSpeed_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("J", string(avg)));
writematrix(sum(avgDeviation_arr)/length(avgDeviation_arr), fullfile(append(path, "/", filename), append(filename, ext)), 'Sheet', 1, 'Range',  strcat("K", string(avg)));

function [timeData, xData, yData, zData, orgX, orgY, orgZ, endX, endY, endZ, supination, flexion, abduction, deviation, averageDeviation] = getdata(csvPath)
	% Import raw csv file from robot.
	opts = detectImportOptions(csvPath);
	M = readtable(csvPath, opts);

	% Calculate size of the imported data table.
	[n, m] = size(M);

	% Extract data from table to array/matrix.
	dataMatrix = M{5:n-2, 1:m-1};
	% averageDeviation = M{n, m-1}

	startX = M{1, 2};
	startY = M{1, 3};
	startZ = M{1, 4};

	endX = M{3, 2};
	endY = M{3, 3};
	endZ = M{3, 4};

	orgX = [startX, endX];
	orgY = [startY, endY];
	orgZ = [startZ, endZ];

	% Put data in column vectors for later use.
	timeData = dataMatrix(:, 1);

	xData = dataMatrix(:, 2);
	yData = dataMatrix(:, 3);
	zData = dataMatrix(:, 4);

	supination = dataMatrix(:, 5);
	flexion = dataMatrix(:, 6);
	abduction = dataMatrix(:, 7);
	deviation = dataMatrix(:, 8);

	% deviation(isnan(deviation)) = 0;
	index = 1;
	while index < 5
		if isnan(deviation(index))
			deviation = deviation(index+1:length(deviation));
		else
			index = index + 1;
		end
	end

	averageDeviation = sum(deviation)/length(deviation);
end

function [speedX, speedY, speedZ] = calculateAvgSpeed(timeData, xData, yData, zData)
	% Get vector length.
	n = size(timeData, 1);

	% First point has no average speed.
	speedX(1) = 0;
	speedY(1) = 0;
	speedZ(1) = 0;

	% Calculate average speed between this and previous point.
	for t=2:n
		deltaT = (timeData(t) - timeData(t-1))/1000;

		speedX(t) = (xData(t) - xData(t-1))/deltaT;
		speedY(t) = (yData(t) - yData(t-1))/deltaT;
		speedZ(t) = (zData(t) - zData(t-1))/deltaT;
	end

	% Make column vectors again.
	speedX = speedX';
	speedY = speedY';
	speedZ = speedZ';
end

function [accX, accY, accZ] = calculateAvgAccel(timeData, xSpeed, ySpeed, zSpeed)
	% Get vector length.
	n = size(timeData, 1);

	% First point has no average speed.
	accX(1) = 0;
	accY(1) = 0;
	accZ(1) = 0;

	% Calculate average speed between this and previous point.
	for t=3:n
		deltaT = (timeData(t) - timeData(t-1))/1000;

		accX(t) = (xSpeed(t) - xSpeed(t-1))/deltaT;
		accY(t) = (ySpeed(t) - ySpeed(t-1))/deltaT;
		accZ(t) = (zSpeed(t) - zSpeed(t-1))/deltaT;
	end

	% Make column vectors again.
	accX = accX';
	accY = accY';
	accZ = accZ';
end

function [fig] = plotSpeed(time, x, y, z, num)
	figure(num)
	grid on
	hold on
	fig = plot3(x, y, z, '-o', 'Color', 'r');
	hold off
	title('Patients speed in space -- 3D')
	xlabel('v_x')
	ylabel('v_y')
	zlabel('v_z')
end

function [fig] = plotPosition(time, x, y, z, orgX, orgY, orgZ, num)
	figure(num);
	grid on
	hold on
	fig = plot3(x, y, z, '-o', 'Color', 'b');
	plot3(orgX, orgY, orgZ, '-o', 'Color', 'm');
	hold off
	title('Patients trajectory in space -- 3D')
	xlabel('x')
	ylabel('y')
	zlabel('z')
	legend('Patients trajectory', 'Correct trajectory')
end

function [spMaxCnt] = getSpeedMaxCnt(a)
	currentMinima = a(1);
	currentMaxima = a(1);
	minimaReached = false;

	maxima = [];

	for index=2:1:length(a)
		if (a(index) < currentMinima) & (minimaReached == false)
			currentMinima = a(index);
		elseif (a(index) > currentMaxima) & (minimaReached == true)
			currentMaxima = a(index);
		end
		
		if (a(index) > currentMinima) & (minimaReached == false) & (index ~= 1)
			currentMaxima = currentMinima;
			minimaReached = true;
		elseif (a(index) < currentMaxima) & (minimaReached == true) & (index ~= 1)
			maxima = [maxima currentMaxima];

			currentMinima = currentMaxima;
			minimaReached = false;
		end
	end

	spMaxCnt = length(maxima);
end

function [mmSpeed] = getMinMaxSpeed(time, v)
	maxSpeed = max(v);
	indexMaxSpeed = 0;

	mmSpeed = [];

	for index=1:1:length(time)
		if v(index) == maxSpeed
			indexMaxSpeed = index;
		end
	end

	[ttmin, min_, ttmax, max_] = getMin(indexMaxSpeed, time, v);

	if length(ttmin) ~= 0
		if length(ttmin) > length(ttmax)
			for index=1:1:length(ttmax)
				mmSpeed = [mmSpeed abs(min_(index) - max_(index)) abs(min_(index + 1) - max_(index))];
			end
		elseif length(ttmin) < length(ttmax)
			for index=1:1:length(ttmin)
				mmSpeed = [mmSpeed abs(max_(index) - min_(index)) abs(max_(index + 1) - min_(index))];
			end
		else
			for index=1:1:(length(ttmin) - 1)
				mmSpeed = [mmSpeed abs(max_(index) - min_(index)) abs(max_(index + 1) - min_(index))];
			end

			index = length(ttmin);
			mmSpeed = [mmSpeed abs(max_(index) - min_(index))];
		end
	else
		mmSpeed = -1;
	end
end

function [timeToMinima, minima, timeToMaxima, maxima] = getMin(startIndex, time, v)
	currentMinima = v(startIndex);
	currentMaxima = v(startIndex);
	minimaReached = false;

	minima = [];
	maxima = [];
	timeToMinima = [];
	timeToMaxima = [];	

	for index=startIndex+1:1:length(time)
		if (v(index) < currentMinima) & (minimaReached == false)
			currentMinima = v(index);
		elseif (v(index) > currentMaxima) & (minimaReached == true)
			currentMaxima = v(index);
		end
		
		if (v(index) > currentMinima) & (minimaReached == false) & (index ~= startIndex)
			minima = [minima currentMinima];
			timeToMinima = [timeToMinima time(index)];

			currentMaxima = currentMinima;
			minimaReached = true;
		elseif (v(index) < currentMaxima) & (minimaReached == true) & (index ~= startIndex)
			maxima = [maxima currentMaxima];
			timeToMaxima = [timeToMaxima time(index)];

			currentMinima = currentMaxima;
			minimaReached = false;
		end
	end
end

function [angle, distanceRatio, speedRatio, totalDistance] = getFMParameters(time, x, y, z, orgX, orgY, orgZ, velocity)
	%	Initial angle calculation
	firstTrajectory = [x(2) - x(1); y(2) - y(1); z(2) - z(1)];
	path = [orgX(2) - orgX(1); orgY(1); orgZ(2) - orgZ(1)];
	angle = radtodeg(acos(dot(firstTrajectory, path)/(norm(firstTrajectory)*norm(path))));

	%	Speed ratio calculation
	initialMaxSpeed = 0;
	[ttmin, min_, ttmax, max_] = getMin(1, time, velocity);

	if ttmax(1) < ttmin(1)
		initialMaxSpeed = max_(1);
	elseif ttmax(1) > ttmin(1)
		index = 1;
		initialMaxSpeed = 0;

		while initialMaxSpeed == 0
			initialMaxSpeed = velocity(index);
			index = index + 1;
		end
	end

	speedRatio = initialMaxSpeed/max(velocity);

	%	Distance ration calculation
	firstStageDistance = 0;
	totalDistance = 0;

	index = 2;
	while ttmin(1) ~= time(index)
		firstStageDistance = firstStageDistance + sqrt((x(index) - x(index-1))^2 + (y(index) - y(index-1))^2 + (z(index) - z(index-1))^2);
		index = index + 1;
	end

	for index=2:1:length(time)
		totalDistance = totalDistance + sqrt((x(index) - x(index-1))^2 + (y(index) - y(index-1))^2 + (z(index) - z(index-1))^2);
	end

	distanceRatio = firstStageDistance/totalDistance;
end

function [reached] = hasReachedTarget(x, y, z, endX, endY, endZ)
	if (x(length(x)) == endX) & (y(length(x)) == endY) & (z(length(x)) == endZ)
		reached = true;
	end

	reached = false;
end

function [reactTime] = getReactionTime(time)
	reactTime = time(2);
end

% deviation
% Overlap same trajectories of different days