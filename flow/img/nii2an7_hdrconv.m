an7hdr.sizeof_hdr = niihdr.sizeof_hdr;
an7hdr.extents = niihdr.extents;
an7hdr.regular = niihdr.regular;
an7hdr.dims = niihdr.dim(1);
an7hdr.xdim = niihdr.dim(2);
an7hdr.ydim = niihdr.dim(3);
an7hdr.zdim = niihdr.dim(4);
an7hdr.tdim = niihdr.dim(5);
an7hdr.datatype = niihdr.datatype;
an7hdr.bits = niihdr.bitpix;
an7hdr.xsize = niihdr.pixdim(2);
an7hdr.ysize = niihdr.pixdim(3);
an7hdr.zsize = niihdr.pixdim(4);
an7hdr.glmax = niihdr.cal_max;
an7hdr.glmin = niihdr.cal_min;