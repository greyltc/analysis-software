[fname, fpath, fltidx]  = uigetfile("*.h5","Choose file ot analyze","../exampleData");

% octave (broken)
%data = load(strcat(fpath,fname))

% matlab
h5disp(strcat(fpath,fname))