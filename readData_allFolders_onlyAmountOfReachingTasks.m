clear all
format long e

basePath_dirs = uigetdir;
save_dir = uigetdir;

dirList = dir(basePath_dirs);

amount_of_reaching_tasks_MS = [];
amount_of_reaching_tasks_Alch = [];

atEnd_dirs = false;
foldernamesDone = {};
index_dirs = 3;

alch = "Alchemist"
ms = "Market Stand"

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
	index = 3;
	amount_of_files = 0;

	while ~atEnd
		fileList = dir(basePath);

		if length(fileList) < index
			atEnd = true;
			continue
		end

		file = fileList(index).name;
		[filepath, name, ext] = fileparts(file);

		if ext ~= ".csv"
			index = index + 1;
			continue
		end

		index = index + 1;
		amount_of_files = amount_of_files + 1;
	end

	if contains(dirList(index_dirs).name, ms)
		amount_of_reaching_tasks_MS = [amount_of_reaching_tasks_MS, amount_of_files];
		amount_of_reaching_tasks_Alch = [amount_of_reaching_tasks_Alch, 0];
	elseif contains(dirList(index_dirs).name, alch)
		amount_of_reaching_tasks_MS = [amount_of_reaching_tasks_MS, 0];
		amount_of_reaching_tasks_Alch = [amount_of_reaching_tasks_Alch, amount_of_files];
	end
	
	foldernamesDone{end+1} = dirList(index_dirs).name;
	index_dirs = index_dirs + 1;
end

filename = "movement_counts";
ext = ".xlsx";

writematrix("folderName", fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range',  "A1");
writematrix("amountOfReachingTasks_MS", fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "B1");
writematrix("amountOfReachingTasks_Alch", fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "C1");
writecell(foldernamesDone', fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range',  "A2");
writematrix(amount_of_reaching_tasks_MS', fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "B2");
writematrix(amount_of_reaching_tasks_Alch', fullfile(save_dir, append(filename, ext)), 'Sheet', 1, 'Range', "C2");