function dirs = subdirectories(path)

    dirs = dir(path);
    dirs = string(setdiff({dirs([dirs(:).isdir]).name}, {'.','..'}));

end