[fname, fpath, fltidx]  = uigetfile("*.h5","Choose file ot analyze","../exampleData");

% octave (broken)
%data = load(strcat(fpath,fname))

target = strcat(fpath,fname);

% matlab
h5disp(target)

info = h5info(target);