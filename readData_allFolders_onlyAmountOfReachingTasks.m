clear all
format long e

basePath_dirs = uigetdir;
save_dir = uigetdir;

dirList = dir(basePath_dirs);
delete(fullfile(save_dir, "*.xlsx"))

amount_of_reaching_tasks = [];

atEnd_dirs = false;
foldernamesDone = [];
index_dirs = 4;

while ~atEnd_dirs
	if length(dirList) < index_dirs
		atEnd_dirs = true;
		continue
	end

	basePath = append(basePath_dirs, "/", dirList(index_dirs).name);

	[filepath, name, ext] = fileparts(basePath);

	if ext ~= ""
		index_dirs = index_dirs + 1;
		continue
	end

	atEnd = false;
	filenamesDone = [];
	index = 3;

	amount_of_files = 0;

	% delete(fullfile(basePath, "*.fig"))
	delete(fullfile(basePath, "*.xlsx"))

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

		filenamesDone = [filenamesDone, file];
		index = index + 1;
		amount_of_files = amount_of_files + 1;
	end

	amount_of_reaching_tasks = [amount_of_reaching_tasks, amount_of_files];

	foldernamesDone = [foldernamesDone, dirList(index_dirs).name];
	index_dirs = index_dirs + 1;
end

filename = "averages_over_all_games";
ext = ".xlsx";

writematrix("amountOfReachingTasks", fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "A1");
writematrix(amount_of_reaching_tasks', fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "A2");

avg_dirs = length(amount_of_reaching_tasks) + 3;
writematrix(sum(amount_of_reaching_tasks), fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', strcat("A", string(avg_dirs)));