function dat = read_fortran_box(fname, nx, ny, nz, dtype)
fid =    fopen(fname,'r');
dat = fread(fid,nx*ny*nz,dtype);
dat = reshape(dat,[nx ny nz]);
fclose(fid);