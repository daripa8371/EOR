function make_movie(moviename,dirname)

aviobj = avifile(moviename,'compression','None'); 
filename = dir(strcat(dirname, '/*.jpg'));  
for i = 1 : length(filename)
    frame = im2frame(imread(strcat(dirname, '/', filename(i).name)));
    aviobj = addframe(aviobj, frame);
end
aviobj = close(aviobj);



