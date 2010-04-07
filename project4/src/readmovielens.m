function [X, timestamps, user, movie, genre]=readmovielens(pathstr, filename)
% [X, timestamps, user, movie, genre]=readmovielens(filename)
%
%  X - User-Item (Movie) matrix
%  timestamps - Matrix of timestamps for rankings
%  user - Cell arrays with user demograhic data
%  movie - Cell arrays with movie data
%  genre - Cell arrays with genre
%  pathstr - Path to MovieLens data files
%  filename - Name of data file, e.g. u.data, u1.base, etc.


% Load the data
data = dlmread([pathstr filename],'\t');

% Construct the User-Item matrix
X=zeros(943,1682); % User-Item (Movie) matrix. 
timestamps=zeros(943,1682); % Time stamps for rankings
for i=1:size(data,1)
    X(data(i,1), data(i,2)) = data(i,3);
    timestamps(data(i,1), data(i,2)) = data(i,4);
end

% Load user demographic data
fid=fopen([pathstr 'u.user']);
user = textscan(fid,'%n%n%s%s%s','Delimiter','|'); % Creates a struct
fclose(fid);

% Load Genre information
fid=fopen([pathstr 'u.genre']);
genre = textscan(fid,'%s%n','Delimiter','|'); % Creates a struct
fclose(fid);

% Load Movie information
fid=fopen([pathstr 'u.item']);
tmp = textscan(fid,'%n%s%s%s%s%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n','Delimiter','|'); % Creates a struct
moviegenrebit = zeros(size(tmp{6},1),19);
for i=1:19
   moviegenrebit(:,i) = tmp{5+i}; 
end

genrestr = cell(size(X,2),1);
for i=1:size(X,2)
    S=genre{1}(find(moviegenrebit(i,:)));
    str = '';
    for j=1:size(S)
        str = strcat(str, S{j}, ',');
    end
    genrestr{i}=str;
end
movie = {tmp{1},tmp{2},tmp{3},tmp{4},tmp{5}, moviegenrebit, genrestr};
fclose(fid);