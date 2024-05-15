function [] = fix_path()

% fix_path: adds to the matlab path the rest of the functions

    addpath("../")                                                                                       %Path to user defined functions
    addpath("../assembly/")
    addpath("../elements/")
    addpath("../errors/")
    addpath("../evaluation/")
    addpath("../matrices/")
    addpath("../mesh_files/")
    addpath("../preprocessing/")
    addpath("../quadrature/")
    addpath("../test/")
    addpath("../utility")
    %addpath("../blendenpik")

end