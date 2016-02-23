#include "stdio.h"
#include "tiffio.h"

#ifdef __MPI
#include <mpi.h>
#endif

TIFF* tif;

tdata_t buf;
tdata_t tmpbuf;

// uint16 vDataType;
// typedef uint16 tdir_t;
uint16 tmpDataType;

tdir_t tdirt; /* directory index */

toff_t nextdir;
// TIFFErrorHandler oldhandler;


typedef struct {
  uint32 width;
  uint32 length;
  uint32 depth;
  uint16 datatype;
  uint16 bitsPerSample;
  uint16 compression;
  uint16 sampleFormat;
  uint16 samplesPerPixel;
  uint16 photometric;
  uint16 fillOrder;
  uint16 planarConfig;
  uint16 resolutionUnit;
  float  xResolution;
  float  yResolution;
} tiffinfo_ct;

tiffinfo_ct tiffinfo[4];

uint32 nb;
// counter for the tiff files
//you can read 4 different tiff files, but will use the
//the first one for writing the TIFF files

void libtiff_set_tiffnumber_(int *nbimage) {
  nb = (uint32) *nbimage;
}

void libtiff_open_tiff_(char *filename, int *infoC) {
//   oldhandler = TIFFSetWarningHandler(NULL);
  tif = TIFFOpen(filename, "r");
  if (!tif) {
    printf("Open File Failed!\n");
    *infoC = -1;
    return;
  }

/*
  int y=10;
  y=TIFFIsBigTIFF(tif);
  printf("File not open! %d",y);*/

//   TIFFSetWarningHandler(oldhandler);
//   TIFFGetField(tif, TIFFTAG_DATATYPE, &vDataType);
  buf = _TIFFmalloc(TIFFScanlineSize(tif));
}


void libtiff_open_write_tiff_(char *filename, int *infoC) {

  tif = TIFFOpen(filename, "a+");
  if (!tif) {
    printf("Open File Failed!\n");
    *infoC = -1;
    return;
  }

//   TIFFGetField(tif, TIFFTAG_DATATYPE, &vDataType);

}

void libtiff_close_tiff_(int *infoC) {

  if (buf != NULL) {
    _TIFFfree(buf);
    buf=NULL;
  }

  TIFFClose(tif);

  *infoC = 0;
  return;

}


void libtiff_read_tiff_info_(int *ngrid,int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  // This is where we read the header of the .TIF file and
  // obtain all required information
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,     &(tiffinfo[nb].width));
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH,    &(tiffinfo[nb].length));
  TIFFGetField(tif, TIFFTAG_IMAGEDEPTH,     &(tiffinfo[nb].depth));
  TIFFGetField(tif, TIFFTAG_DATATYPE,       &(tiffinfo[nb].datatype));
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,  &(tiffinfo[nb].bitsPerSample));
  TIFFGetField(tif, TIFFTAG_COMPRESSION,    &(tiffinfo[nb].compression));
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL,&(tiffinfo[nb].samplesPerPixel));
  TIFFGetField(tif, TIFFTAG_PHOTOMETRIC,    &(tiffinfo[nb].photometric));
  TIFFGetField(tif, TIFFTAG_FILLORDER,      &(tiffinfo[nb].fillOrder));
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG,   &(tiffinfo[nb].planarConfig));
  TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &(tiffinfo[nb].resolutionUnit));
  TIFFGetField(tif, TIFFTAG_XRESOLUTION,    &(tiffinfo[nb].xResolution));
  TIFFGetField(tif, TIFFTAG_YRESOLUTION,    &(tiffinfo[nb].yResolution));
  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,   &(tiffinfo[nb].sampleFormat));

  // force this one for now...
  if (tiffinfo[nb].sampleFormat <= 0) tiffinfo[nb].sampleFormat = 1;

  // Read the x and y values from the image and assign those to Ngrid
  // Be careful with the x/y assignment
  ngrid[0] = tiffinfo[nb].width;
  ngrid[1] = tiffinfo[nb].length;
  ngrid[2] = (int) TIFFNumberOfDirectories(tif);
}

void libtiff_read_tiff_info__(int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  // This is where we read the header of the .TIF file and
  // obtain all required information
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH,     &(tiffinfo[nb].width));
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH,    &(tiffinfo[nb].length));
  TIFFGetField(tif, TIFFTAG_IMAGEDEPTH,     &(tiffinfo[nb].depth));
  TIFFGetField(tif, TIFFTAG_DATATYPE,       &(tiffinfo[nb].datatype));
  TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE,  &(tiffinfo[nb].bitsPerSample));
  TIFFGetField(tif, TIFFTAG_COMPRESSION,    &(tiffinfo[nb].compression));
  TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL,&(tiffinfo[nb].samplesPerPixel));
  TIFFGetField(tif, TIFFTAG_PHOTOMETRIC,    &(tiffinfo[nb].photometric));
  TIFFGetField(tif, TIFFTAG_FILLORDER,      &(tiffinfo[nb].fillOrder));
  TIFFGetField(tif, TIFFTAG_PLANARCONFIG,   &(tiffinfo[nb].planarConfig));
  TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &(tiffinfo[nb].resolutionUnit));
  TIFFGetField(tif, TIFFTAG_XRESOLUTION,    &(tiffinfo[nb].xResolution));
  TIFFGetField(tif, TIFFTAG_YRESOLUTION,    &(tiffinfo[nb].yResolution));
  TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT,   &(tiffinfo[nb].sampleFormat));

  // force this one for now...
  if (tiffinfo[nb].sampleFormat <= 0) tiffinfo[nb].sampleFormat = 1;
}


#ifdef __MPI
void libtiff_read_tiff_info_bcast_(MPI_Fint *comm, int *infoC)
{
  MPI_Status status;
  MPI_Datatype tiffinfoMPI;
  MPI_Datatype type[3]={MPI_UNSIGNED,MPI_UNSIGNED_SHORT,MPI_FLOAT};
  int blocklen[3]={3,9,2};
  MPI_Aint disp[3];
  MPI_Aint base;
  MPI_Comm ccomm;
  int i;

  ccomm=MPI_Comm_f2c(*comm);

  MPI_Get_address(&tiffinfo[nb],disp);
  MPI_Get_address(&tiffinfo[nb].datatype,disp+1);
  MPI_Get_address(&tiffinfo[nb].xResolution,disp+2);

  base = disp[0];
  for (i=0; i<3; i++) disp[i] -= base;

  MPI_Type_create_struct(3,blocklen,disp,type,&tiffinfoMPI);
  MPI_Type_commit(&tiffinfoMPI);
  MPI_Bcast(&tiffinfo[nb],1,tiffinfoMPI,0,ccomm);
  MPI_Type_free(&tiffinfoMPI);
// //   MPI_Status status[14];
// //   MPI_Request request[14];
// //   MPI_Comm ccomm;
// //   ccomm=MPI_Comm_f2c(*comm);
// //   MPI_Ibcast(&tiffinfo[nb].width,1,MPI_UNSIGNED,0,ccomm,&request[0]);
// //   MPI_Ibcast(&tiffinfo[nb].length,1,MPI_UNSIGNED,0,ccomm,&request[1]);
// //   MPI_Ibcast(&tiffinfo[nb].depth,1,MPI_UNSIGNED,0,ccomm,&request[2]);
// //   MPI_Ibcast(&tiffinfo[nb].datatype,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[3]);
// //   MPI_Ibcast(&tiffinfo[nb].bitsPerSample,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[4]);
// //   MPI_Ibcast(&tiffinfo[nb].compression,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[5]);
// //   MPI_Ibcast(&tiffinfo[nb].sampleFormat,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[6]);
// //   MPI_Ibcast(&tiffinfo[nb].samplesPerPixel,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[7]);
// //   MPI_Ibcast(&tiffinfo[nb].photometric,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[8]);
// //   MPI_Ibcast(&tiffinfo[nb].fillOrder,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[9]);
// //   MPI_Ibcast(&tiffinfo[nb].planarConfig,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[10]);
// //   MPI_Ibcast(&tiffinfo[nb].resolutionUnit,1,MPI_UNSIGNED_SHORT,0,ccomm,&request[11]);
// //   MPI_Ibcast(&tiffinfo[nb].xResolution,1,MPI_FLOAT,0,ccomm,&request[12]);
// //   MPI_Ibcast(&tiffinfo[nb].yResolution,1,MPI_FLOAT,0,ccomm,&request[13]);
// //   MPI_Waitall(14, request, status);

}
#endif

void libtiff_write_tiff_header_2d_(int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  char software[] = "MOSAIC Group ppm_rc";

  if (tiffinfo[nb].compression == 0)
    tiffinfo[nb].compression = COMPRESSION_NONE;
  if (tiffinfo[nb].bitsPerSample == 0)
    tiffinfo[nb].bitsPerSample = 8;
  if (tiffinfo[nb].fillOrder == 0)
    tiffinfo[nb].fillOrder = FILLORDER_MSB2LSB;
  if (tiffinfo[nb].resolutionUnit == 0)
    tiffinfo[nb].resolutionUnit = RESUNIT_NONE;

  // This is where we read the header of the .TIFF file and
  // obtain all required information
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     tiffinfo[nb].width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    tiffinfo[nb].length);
  TIFFSetField(tif, TIFFTAG_IMAGEDEPTH,     tiffinfo[nb].depth);
  TIFFSetField(tif, TIFFTAG_DATATYPE,       tiffinfo[nb].datatype);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  tiffinfo[nb].bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    tiffinfo[nb].compression);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   tiffinfo[nb].sampleFormat);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,tiffinfo[nb].samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    tiffinfo[nb].photometric);
  TIFFSetField(tif, TIFFTAG_FILLORDER,      tiffinfo[nb].fillOrder);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   tiffinfo[nb].planarConfig);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, tiffinfo[nb].resolutionUnit);
  TIFFSetField(tif, TIFFTAG_XRESOLUTION,    tiffinfo[nb].xResolution);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION,    tiffinfo[nb].yResolution);

  TIFFSetField(tif, TIFFTAG_SOFTWARE,  software);

  _TIFFfree(buf);
  buf=NULL;

  buf = _TIFFmalloc(TIFFScanlineSize(tif));

}


void libtiff_write_tiff_header_3d_(int *infoC, int *page, int *npages)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  char software[] = "MOSAIC Group ppm_rc";

  if (tiffinfo[nb].compression == 0)
    tiffinfo[nb].compression = COMPRESSION_NONE;
  if (tiffinfo[nb].bitsPerSample == 0)
    tiffinfo[nb].bitsPerSample = 8;
  if (tiffinfo[nb].fillOrder == 0)
    tiffinfo[nb].fillOrder = FILLORDER_MSB2LSB;
  if (tiffinfo[nb].resolutionUnit == 0)
    tiffinfo[nb].resolutionUnit = RESUNIT_NONE;

  // This is where we read the header of the .TIF file and
  // obtain all required information
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     tiffinfo[nb].width);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    tiffinfo[nb].length);
  TIFFSetField(tif, TIFFTAG_IMAGEDEPTH,     tiffinfo[nb].depth);
  TIFFSetField(tif, TIFFTAG_DATATYPE,       tiffinfo[nb].datatype);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  tiffinfo[nb].bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    tiffinfo[nb].compression);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   tiffinfo[nb].sampleFormat);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,tiffinfo[nb].samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    tiffinfo[nb].photometric);
  TIFFSetField(tif, TIFFTAG_FILLORDER,      tiffinfo[nb].fillOrder);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   tiffinfo[nb].planarConfig);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, tiffinfo[nb].resolutionUnit);
  TIFFSetField(tif, TIFFTAG_XRESOLUTION,    tiffinfo[nb].xResolution);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION,    tiffinfo[nb].yResolution);

  /* We are writing single page of the multipage file */
  TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
  /* Set the page number */
  TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, npages);

  TIFFSetField(tif, TIFFTAG_SOFTWARE,  software);

  _TIFFfree(buf);
  buf=NULL;

  buf = _TIFFmalloc(TIFFScanlineSize(tif));

}

void libtiff_write_tiff_customheader_2d_(int *infoC,int *width,int *length)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 mywidth=(uint32)*width;
  uint32 mylength=(uint32)*length;

  char software[] = "MOSAIC Group ppm_rc";

  if (tiffinfo[nb].compression == 0)
    tiffinfo[nb].compression = COMPRESSION_NONE;
  if (tiffinfo[nb].bitsPerSample == 0)
    tiffinfo[nb].bitsPerSample = 8;
  if (tiffinfo[nb].fillOrder == 0)
    tiffinfo[nb].fillOrder = FILLORDER_MSB2LSB;
  if (tiffinfo[nb].resolutionUnit == 0)
    tiffinfo[nb].resolutionUnit = RESUNIT_NONE;

  // This is where we read the header of the .TIF file and
  // obtain all required information
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     mywidth);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    mylength);
  TIFFSetField(tif, TIFFTAG_IMAGEDEPTH,     tiffinfo[nb].depth);
  TIFFSetField(tif, TIFFTAG_DATATYPE,       tiffinfo[nb].datatype);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  tiffinfo[nb].bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    tiffinfo[nb].compression);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   tiffinfo[nb].sampleFormat);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,tiffinfo[nb].samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    tiffinfo[nb].photometric);
  TIFFSetField(tif, TIFFTAG_FILLORDER,      tiffinfo[nb].fillOrder);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   tiffinfo[nb].planarConfig);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, tiffinfo[nb].resolutionUnit);
  TIFFSetField(tif, TIFFTAG_XRESOLUTION,    tiffinfo[nb].xResolution);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION,    tiffinfo[nb].yResolution);

  TIFFSetField(tif, TIFFTAG_SOFTWARE,  software);

  _TIFFfree(buf);
  buf=NULL;

  buf = _TIFFmalloc(TIFFScanlineSize(tif));

}


void libtiff_write_tiff_customheader_3d_(int *infoC,int *page,int *npages,int * width,int *length)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 mywidth=(uint32) *width;
  uint32 mylength=(uint32) *length;

  char software[] = "MOSAIC Group ppm_rc";

  if (tiffinfo[nb].compression == 0)
    tiffinfo[nb].compression = COMPRESSION_NONE;
  if (tiffinfo[nb].bitsPerSample == 0)
    tiffinfo[nb].bitsPerSample = 8;
  if (tiffinfo[nb].fillOrder == 0)
    tiffinfo[nb].fillOrder = FILLORDER_MSB2LSB;
  if (tiffinfo[nb].resolutionUnit == 0)
    tiffinfo[nb].resolutionUnit = RESUNIT_NONE;

  // This is where we read the header of the .TIF file and
  // obtain all required information
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH,     mywidth);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH,    mylength);
  TIFFSetField(tif, TIFFTAG_IMAGEDEPTH,     tiffinfo[nb].depth);
  TIFFSetField(tif, TIFFTAG_DATATYPE,       tiffinfo[nb].datatype);
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE,  tiffinfo[nb].bitsPerSample);
  TIFFSetField(tif, TIFFTAG_COMPRESSION,    tiffinfo[nb].compression);
  TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,   tiffinfo[nb].sampleFormat);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL,tiffinfo[nb].samplesPerPixel);
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC,    tiffinfo[nb].photometric);
  TIFFSetField(tif, TIFFTAG_FILLORDER,      tiffinfo[nb].fillOrder);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG,   tiffinfo[nb].planarConfig);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, tiffinfo[nb].resolutionUnit);
  TIFFSetField(tif, TIFFTAG_XRESOLUTION,    tiffinfo[nb].xResolution);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION,    tiffinfo[nb].yResolution);

  /* We are writing single page of the multipage file */
  TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
  /* Set the page number */
  TIFFSetField(tif, TIFFTAG_PAGENUMBER, page, npages);

  TIFFSetField(tif, TIFFTAG_SOFTWARE, software);

  _TIFFfree(buf);
  buf=NULL;

  buf = _TIFFmalloc(TIFFScanlineSize(tif));

}


void libtiff_read_tiff_directory_(int *scandir,int *dirtn,int *xMax,int *yMax,int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  tsample_t vSample=1;

  tdir_t dirn = (tdir_t)*dirtn;

  uint32 myxMax=(uint32)*xMax;
  uint32 myyMax=(uint32)*yMax;

  uint32 x,y;

  TIFFSetDirectory(tif, dirn);

  switch (tiffinfo[nb].bitsPerSample) {
  case 8:
    for (y=0;y<myyMax;y++){
      TIFFReadScanline(tif, buf, y, vSample);
      for (x=0; x<myxMax; x++){
        scandir[x+y*myxMax] = *((unsigned char*)buf+x);
      }
    }
    break;

  case 16:
    for (y=0;y<myyMax;y++){
      TIFFReadScanline(tif, buf, y, vSample);
      for (x=0; x<myxMax; x++){
        scandir[x+y*myxMax] = *((unsigned short*)buf+x);
      }
    }
    break;

  case 32:
    for (y=0;y<myyMax;y++){
      TIFFReadScanline(tif, buf, y, vSample);
      for (x=0; x<myxMax; x++){
        scandir[x+y*myxMax] = *((unsigned int*)buf+x);
      }
    }
    break;
  }

}

void libtiff_read_tiff_scanline_(int *scanline,int *linenm,int *xMax,int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 x;
  uint32 row=(uint32)*linenm;

  tsample_t vSample=1;

  uint32 myxMax=(uint32)*xMax;

  TIFFReadScanline(tif, buf, row, vSample);

  switch (tiffinfo[nb].bitsPerSample) {
  case 8:
    for (x=0; x<myxMax; x++){
      scanline[x] = *((unsigned char*)buf+x);
    }
    break;

  case 16:
    for (x=0; x<myxMax; x++){
      scanline[x] = *((unsigned short*)buf+x);
    }
    break;

  case 32:
    for (x=0; x<myxMax; x++){
      scanline[x] = *((unsigned int*)buf+x);
    }
    break;

  }

}

void libtiff_read_tiff_scanline_int_(int *scanline,int *linenm,int *xMax,int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 x;
  uint32 row=(uint32)*linenm;

  tsample_t vSample=1;

  uint32 myxMax=(uint32)*xMax;

  TIFFReadScanline(tif, buf, row, vSample);

  switch (tiffinfo[nb].bitsPerSample) {
    case 8:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned char*)buf+x);
      }
      break;

    case 16:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned short*)buf+x);
      }
      break;

    case 32:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned int*)buf+x);
      }
      break;

  }

}

void libtiff_read_tiff_scanline_double_(double *scanline,int *linenm,int *xMax,int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 x;
  uint32 row=(uint32)*linenm;

  tsample_t vSample=1;

  uint32 myxMax=(uint32)*xMax;

  TIFFReadScanline(tif, buf, row, vSample);

  switch (tiffinfo[nb].bitsPerSample) {
    case 8:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned char*)buf+x);
      }
      break;

    case 16:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned short*)buf+x);
      }
      break;

    case 32:
      for (x=0; x<myxMax; x++){
        scanline[x] = *((unsigned int*)buf+x);
      }
      break;

  }

}

void libtiff_write_tiff_directory_(int* scandir, int *dirtn, int *xMax, int *yMax, int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  tsample_t vSample=1;

  tdir_t dirn = (tdir_t)*dirtn;

  uint32 myxMax=(uint32)*xMax;
  uint32 myyMax=(uint32)*yMax;

  uint32 x,y;

  TIFFSetDirectory(tif, dirn);

  switch (tiffinfo[nb].bitsPerSample) {
  case 8:
    for (y=0;y<myyMax;y++){
      for (x=0; x<myxMax; x++){
        *((unsigned char*)buf+x) = (unsigned char) scandir[x+y*myxMax];
      }
      TIFFWriteScanline(tif, buf, y, vSample);
    }
    break;

  case 16:
    for (y=0;y<myyMax;y++){
      for (x=0; x<myxMax; x++){
        *((unsigned short*)buf+x) = (unsigned short) scandir[x+y*myxMax];
      }
      TIFFWriteScanline(tif, buf, y, vSample);
    }
    break;

  case 32:
    for (y=0;y<myyMax;y++){
      for (x=0; x<myxMax; x++){
        *((unsigned int*)buf+x) = (unsigned int) scandir[x+y*myxMax];
      }
      TIFFWriteScanline(tif, buf, y, vSample);
    }
    break;

  }

  TIFFWriteDirectory(tif);

  _TIFFfree(buf);
  buf=NULL;

}

void libtiff_write_tiff_scanline_(int* scanline, int *linenm, int *xMax,  int *infoC)
{

  if (!tif) {
    printf("File not open!\n");
    *infoC = -1;
    return;
  }

  uint32 x;
  uint32 row=(uint32)*linenm;
  tsample_t vSample=1;
  uint32 myxMax=(uint32)*xMax;

  switch (tiffinfo[nb].bitsPerSample) {
  case 8:
    for (x=0; x<myxMax; x++) {
      *((unsigned char*)buf+x) = (unsigned char) scanline[x];
    }
    break;

  case 16:
    for (x=0; x<myxMax; x++) {
      *((unsigned short*)buf+x) = (unsigned short) scanline[x];
    }
    break;

  case 32:
    for (x=0; x<myxMax; x++) {
      *((unsigned int*)buf+x) = (unsigned int) scanline[x];
    }
    break;

  }

  //_TIFFmemcpy(wbuf,cbuf,myxMax);
  TIFFWriteScanline(tif, buf, row, vSample);

}


// // // // TIFF* TIFFOpenMPI(const char* name, const char* mode)
// // // // {
// // // //   static const char module[] = "TIFFOpenMPI";
// // // //   int m, fd;
// // // //   MPI_File fh;
// // // //   int MPImode;
// // // //
// // // //   m = _TIFFgetMode(mode, module);
// // // //   if (m == -1) return ((TIFF*)0);
// // // //
// // // //   switch (mode[0]) {
// // // //     case 'w':
// // // //       MPImode = MPI_MODE_WRONLY;
// // // //     default:
// // // //       TIFFErrorExt(0, module, "\"%s\": Bad mode", mode);
// // // //       break;
// // // //   }
// // // //
// // // //   fd = MPI_File_open(MPI_COMM_WORLD, name, MPImode, MPI_INFO_NULL, &fh);
// // // //   if (fd < 0) {
// // // //     printf( "Unable to open file \"temp\"\n" );
// // // //     fflush(stdout);
// // // //     return ((TIFF *)0);
// // // //   }
// // // //
// // // //   tif = TIFFFdOpenMPI((int)fd, name, mode);
// // // //   if (!tif) close(fd);
// // // //   return tif;
// // // // }

// // // // // TIFF* TIFFFdOpenMPI(int fd, const char* name, const char* mode)
// // // // // {
// // // // //   tif = TIFFClientOpenMPI(name, mode,
// // // // //                        (thandle_t) fd,
// // // // //                        _tiffReadProc, _tiffWriteProc,
// // // // //                        _tiffSeekProc, _tiffCloseProc, _tiffSizeProc,
// // // // //                        _tiffMapProc, _tiffUnmapProc);
// // // // //   if (tif)
// // // // //     tif->tif_fd = fd;
// // // // //   return (tif);
// // // // // }
// // // // //
// // // // //
// // // // // TIFF* TIFFClientOpenMPI(
// // // // //   const char* name, const char* mode,
// // // // //   thandle_t clientdata,
// // // // //   TIFFReadWriteProc readproc,
// // // // //   TIFFReadWriteProc writeproc,
// // // // //   TIFFSeekProc seekproc,
// // // // //   TIFFCloseProc closeproc,
// // // // //   TIFFSizeProc sizeproc,
// // // // //   TIFFMapFileProc mapproc,
// // // // //   TIFFUnmapFileProc unmapproc
// // // // // )
// // // // // {
// // // // //   static const char module[] = "TIFFClientOpenMPI";
// // // // //   int m;
// // // // //   const char* cp;
// // // // //
// // // // //   tif = (TIFF *)_TIFFmalloc((tmsize_t)(sizeof (TIFF) + strlen(name) + 1));
// // // // //   if (tif == NULL) {
// // // // //     TIFFErrorExt(clientdata, module, "%s: Out of memory (TIFF structure)", name);
// // // // //     goto bad2;
// // // // //   }
// // // // //   _TIFFmemset(tif, 0, sizeof (*tif));
// // // // //   tif->tif_name = (char *)tif + sizeof (TIFF);
// // // // //   strcpy(tif->tif_name, name);
// // // // //   tif->tif_mode = m &~ (O_CREAT|O_TRUNC);
// // // // //   tif->tif_curdir = (uint16) -1;          /* non-existent directory */
// // // // //   tif->tif_curoff = 0;
// // // // //   tif->tif_curstrip = (uint32) -1;        /* invalid strip */
// // // // //   tif->tif_row = (uint32) -1;             /* read/write pre-increment */
// // // // //   tif->tif_clientdata = clientdata;
// // // // //   if (!readproc || !writeproc || !seekproc || !closeproc || !sizeproc) {
// // // // //     TIFFErrorExt(clientdata, module,
// // // // //                  "One of the client procedures is NULL pointer.");
// // // // //     goto bad2;
// // // // //   }
// // // // //   tif->tif_readproc = readproc;
// // // // //   tif->tif_writeproc = writeproc;
// // // // //   tif->tif_seekproc = seekproc;
// // // // //   tif->tif_closeproc = closeproc;
// // // // //   tif->tif_sizeproc = sizeproc;
// // // // //   if (mapproc)
// // // // //     tif->tif_mapproc = mapproc;
// // // // //   else
// // // // //     tif->tif_mapproc = _tiffDummyMapProc;
// // // // //   if (unmapproc)
// // // // //     tif->tif_unmapproc = unmapproc;
// // // // //   else
// // // // //     tif->tif_unmapproc = _tiffDummyUnmapProc;
// // // // //   _TIFFSetDefaultCompressionState(tif);    /* setup default state */
// // // // //   /*
// // // // //    * Default is to return data MSB2LSB and enable the
// // // // //    * use of memory-mapped files and strip chopping when
// // // // //    * a file is opened read-only.
// // // // //    */
// // // // //   tif->tif_flags = FILLORDER_MSB2LSB;
// // // // //   if (m == O_RDONLY )
// // // // //     tif->tif_flags |= TIFF_MAPPED;
// // // // //
// // // // //   #ifdef STRIPCHOP_DEFAULT
// // // // //   if (m == O_RDONLY || m == O_RDWR)
// // // // //     tif->tif_flags |= STRIPCHOP_DEFAULT;
// // // // //   #endif
// // // // //
// // // // //   /*
// // // // //    * Process library-specific flags in the open mode string.
// // // // //    * The following flags may be used to control intrinsic library
// // // // //    * behaviour that may or may not be desirable (usually for
// // // // //    * compatibility with some application that claims to support
// // // // //    * TIFF but only supports some braindead idea of what the
// // // // //    * vendor thinks TIFF is):
// // // // //    *
// // // // //    * 'l' use little-endian byte order for creating a file
// // // // //    * 'b' use big-endian byte order for creating a file
// // // // //    * 'L' read/write information using LSB2MSB bit order
// // // // //    * 'B' read/write information using MSB2LSB bit order
// // // // //    * 'H' read/write information using host bit order
// // // // //    * 'M' enable use of memory-mapped files when supported
// // // // //    * 'm' disable use of memory-mapped files
// // // // //    * 'C' enable strip chopping support when reading
// // // // //    * 'c' disable strip chopping support
// // // // //    * 'h' read TIFF header only, do not load the first IFD
// // // // //    * '4' ClassicTIFF for creating a file (default)
// // // // //    * '8' BigTIFF for creating a file
// // // // //    *
// // // // //    * The use of the 'l' and 'b' flags is strongly discouraged.
// // // // //    * These flags are provided solely because numerous vendors,
// // // // //    * typically on the PC, do not correctly support TIFF; they
// // // // //    * only support the Intel little-endian byte order.  This
// // // // //    * support is not configured by default because it supports
// // // // //    * the violation of the TIFF spec that says that readers *MUST*
// // // // //    * support both byte orders.  It is strongly recommended that
// // // // //    * you not use this feature except to deal with busted apps
// // // // //    * that write invalid TIFF.  And even in those cases you should
// // // // //    * bang on the vendors to fix their software.
// // // // //    *
// // // // //    * The 'L', 'B', and 'H' flags are intended for applications
// // // // //    * that can optimize operations on data by using a particular
// // // // //    * bit order.  By default the library returns data in MSB2LSB
// // // // //    * bit order for compatibiltiy with older versions of this
// // // // //    * library.  Returning data in the bit order of the native cpu
// // // // //    * makes the most sense but also requires applications to check
// // // // //    * the value of the FillOrder tag; something they probably do
// // // // //    * not do right now.
// // // // //    *
// // // // //    * The 'M' and 'm' flags are provided because some virtual memory
// // // // //    * systems exhibit poor behaviour when large images are mapped.
// // // // //    * These options permit clients to control the use of memory-mapped
// // // // //    * files on a per-file basis.
// // // // //    *
// // // // //    * The 'C' and 'c' flags are provided because the library support
// // // // //    * for chopping up large strips into multiple smaller strips is not
// // // // //    * application-transparent and as such can cause problems.  The 'c'
// // // // //    * option permits applications that only want to look at the tags,
// // // // //    * for example, to get the unadulterated TIFF tag information.
// // // // //    */
// // // // //   for (cp = mode; *cp; cp++)
// // // // //     switch (*cp) {
// // // // //       case 'b':
// // // // //         #ifndef WORDS_BIGENDIAN
// // // // //         if (m&O_CREAT)
// // // // //           tif->tif_flags |= TIFF_SWAB;
// // // // //         #endif
// // // // //         break;
// // // // //       case 'l':
// // // // //         #ifdef WORDS_BIGENDIAN
// // // // //         if ((m&O_CREAT))
// // // // //           tif->tif_flags |= TIFF_SWAB;
// // // // //         #endif
// // // // //         break;
// // // // //       case 'B':
// // // // //         tif->tif_flags = (tif->tif_flags &~ TIFF_FILLORDER) |
// // // // //         FILLORDER_MSB2LSB;
// // // // //         break;
// // // // //       case 'L':
// // // // //         tif->tif_flags = (tif->tif_flags &~ TIFF_FILLORDER) |
// // // // //         FILLORDER_LSB2MSB;
// // // // //         break;
// // // // //       case 'H':
// // // // //         tif->tif_flags = (tif->tif_flags &~ TIFF_FILLORDER) |
// // // // //         HOST_FILLORDER;
// // // // //         break;
// // // // //       case 'M':
// // // // //         if (m == O_RDONLY)
// // // // //           tif->tif_flags |= TIFF_MAPPED;
// // // // //         break;
// // // // //       case 'm':
// // // // //         if (m == O_RDONLY)
// // // // //           tif->tif_flags &= ~TIFF_MAPPED;
// // // // //         break;
// // // // //       case 'C':
// // // // //         if (m == O_RDONLY)
// // // // //           tif->tif_flags |= TIFF_STRIPCHOP;
// // // // //         break;
// // // // //       case 'c':
// // // // //         if (m == O_RDONLY)
// // // // //           tif->tif_flags &= ~TIFF_STRIPCHOP;
// // // // //         break;
// // // // //       case 'h':
// // // // //         tif->tif_flags |= TIFF_HEADERONLY;
// // // // //         break;
// // // // //       case '8':
// // // // //         if (m&O_CREAT)
// // // // //           tif->tif_flags |= TIFF_BIGTIFF;
// // // // //         break;
// // // // //     }
// // // // //     /*
// // // // //      * Read in TIFF header.
// // // // //      */
// // // // //     if ((m & O_TRUNC) ||
// // // // //       !ReadOK(tif, &tif->tif_header, sizeof (TIFFHeaderClassic))) {
// // // // //       if (tif->tif_mode == O_RDONLY) {
// // // // //         TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                      "Cannot read TIFF header");
// // // // //         goto bad;
// // // // //       }
// // // // //       /*
// // // // //        * Setup header and write.
// // // // //        */
// // // // //       #ifdef WORDS_BIGENDIAN
// // // // //       tif->tif_header.common.tiff_magic = tif->tif_flags & TIFF_SWAB
// // // // //       ? TIFF_LITTLEENDIAN : TIFF_BIGENDIAN;
// // // // //       #else
// // // // //       tif->tif_header.common.tiff_magic = tif->tif_flags & TIFF_SWAB
// // // // //       ? TIFF_BIGENDIAN : TIFF_LITTLEENDIAN;
// // // // //       #endif
// // // // //       if (!(tif->tif_flags&TIFF_BIGTIFF))
// // // // //       {
// // // // //         tif->tif_header.common.tiff_version = TIFF_VERSION_CLASSIC;
// // // // //         tif->tif_header.classic.tiff_diroff = 0;
// // // // //         if (tif->tif_flags & TIFF_SWAB)
// // // // //           TIFFSwabShort(&tif->tif_header.common.tiff_version);
// // // // //         tif->tif_header_size = sizeof(TIFFHeaderClassic);
// // // // //       }
// // // // //       else
// // // // //       {
// // // // //         tif->tif_header.common.tiff_version = TIFF_VERSION_BIG;
// // // // //         tif->tif_header.big.tiff_offsetsize = 8;
// // // // //         tif->tif_header.big.tiff_unused = 0;
// // // // //         tif->tif_header.big.tiff_diroff = 0;
// // // // //         if (tif->tif_flags & TIFF_SWAB)
// // // // //         {
// // // // //           TIFFSwabShort(&tif->tif_header.common.tiff_version);
// // // // //           TIFFSwabShort(&tif->tif_header.big.tiff_offsetsize);
// // // // //         }
// // // // //         tif->tif_header_size = sizeof (TIFFHeaderBig);
// // // // //       }
// // // // //       /*
// // // // //        * The doc for "fopen" for some STD_C_LIBs says that if you
// // // // //        * open a file for modify ("+"), then you must fseek (or
// // // // //        * fflush?) between any freads and fwrites.  This is not
// // // // //        * necessary on most systems, but has been shown to be needed
// // // // //        * on Solaris.
// // // // //        */
// // // // //       TIFFSeekFile( tif, 0, SEEK_SET );
// // // // //       if (!WriteOK(tif, &tif->tif_header, (tmsize_t)(tif->tif_header_size))) {
// // // // //         TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                      "Error writing TIFF header");
// // // // //         goto bad;
// // // // //       }
// // // // //       /*
// // // // //        * Setup the byte order handling.
// // // // //        */
// // // // //       if (tif->tif_header.common.tiff_magic == TIFF_BIGENDIAN) {
// // // // //         #ifndef WORDS_BIGENDIAN
// // // // //         tif->tif_flags |= TIFF_SWAB;
// // // // //         #endif
// // // // //       } else {
// // // // //         #ifdef WORDS_BIGENDIAN
// // // // //         tif->tif_flags |= TIFF_SWAB;
// // // // //         #endif
// // // // //       }
// // // // //       /*
// // // // //        * Setup default directory.
// // // // //        */
// // // // //       if (!TIFFDefaultDirectory(tif))
// // // // //         goto bad;
// // // // //       tif->tif_diroff = 0;
// // // // //       tif->tif_dirlist = NULL;
// // // // //       tif->tif_dirlistsize = 0;
// // // // //       tif->tif_dirnumber = 0;
// // // // //       return (tif);
// // // // //       }
// // // // //       /*
// // // // //        * Setup the byte order handling.
// // // // //        */
// // // // //       if (tif->tif_header.common.tiff_magic != TIFF_BIGENDIAN &&
// // // // //         tif->tif_header.common.tiff_magic != TIFF_LITTLEENDIAN
// // // // //         #if MDI_SUPPORT
// // // // //         &&
// // // // //         #if HOST_BIGENDIAN
// // // // //         tif->tif_header.common.tiff_magic != MDI_BIGENDIAN
// // // // //         #else
// // // // //         tif->tif_header.common.tiff_magic != MDI_LITTLEENDIAN
// // // // //         #endif
// // // // //       ) {
// // // // //         TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                      "Not a TIFF or MDI file, bad magic number %d (0x%x)",
// // // // //                      #else
// // // // //         ) {
// // // // //           TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                        "Not a TIFF file, bad magic number %d (0x%x)",
// // // // //                      #endif
// // // // //                      tif->tif_header.common.tiff_magic,
// // // // //                      tif->tif_header.common.tiff_magic);
// // // // //         goto bad;
// // // // //         }
// // // // //         if (tif->tif_header.common.tiff_magic == TIFF_BIGENDIAN) {
// // // // //           #ifndef WORDS_BIGENDIAN
// // // // //           tif->tif_flags |= TIFF_SWAB;
// // // // //           #endif
// // // // //         } else {
// // // // //           #ifdef WORDS_BIGENDIAN
// // // // //           tif->tif_flags |= TIFF_SWAB;
// // // // //           #endif
// // // // //         }
// // // // //         if (tif->tif_flags & TIFF_SWAB)
// // // // //           TIFFSwabShort(&tif->tif_header.common.tiff_version);
// // // // //         if ((tif->tif_header.common.tiff_version != TIFF_VERSION_CLASSIC)&&
// // // // //           (tif->tif_header.common.tiff_version != TIFF_VERSION_BIG)) {
// // // // //           TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                        "Not a TIFF file, bad version number %d (0x%x)",
// // // // //                      tif->tif_header.common.tiff_version,
// // // // //                      tif->tif_header.common.tiff_version);
// // // // //         goto bad;
// // // // //           }
// // // // //           if (tif->tif_header.common.tiff_version == TIFF_VERSION_CLASSIC)
// // // // //           {
// // // // //             if (tif->tif_flags & TIFF_SWAB)
// // // // //               TIFFSwabLong(&tif->tif_header.classic.tiff_diroff);
// // // // //             tif->tif_header_size = sizeof(TIFFHeaderClassic);
// // // // //           }
// // // // //           else
// // // // //           {
// // // // //             if (!ReadOK(tif, ((uint8*)(&tif->tif_header) + sizeof(TIFFHeaderClassic)), (sizeof(TIFFHeaderBig)-sizeof(TIFFHeaderClassic))))
// // // // //             {
// // // // //               TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                            "Cannot read TIFF header");
// // // // //               goto bad;
// // // // //             }
// // // // //             if (tif->tif_flags & TIFF_SWAB)
// // // // //             {
// // // // //               TIFFSwabShort(&tif->tif_header.big.tiff_offsetsize);
// // // // //               TIFFSwabLong8(&tif->tif_header.big.tiff_diroff);
// // // // //             }
// // // // //             if (tif->tif_header.big.tiff_offsetsize != 8)
// // // // //             {
// // // // //               TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                            "Not a TIFF file, bad BigTIFF offsetsize %d (0x%x)",
// // // // //                            tif->tif_header.big.tiff_offsetsize,
// // // // //                            tif->tif_header.big.tiff_offsetsize);
// // // // //               goto bad;
// // // // //             }
// // // // //             if (tif->tif_header.big.tiff_unused != 0)
// // // // //             {
// // // // //               TIFFErrorExt(tif->tif_clientdata, name,
// // // // //                            "Not a TIFF file, bad BigTIFF unused %d (0x%x)",
// // // // //                            tif->tif_header.big.tiff_unused,
// // // // //                            tif->tif_header.big.tiff_unused);
// // // // //               goto bad;
// // // // //             }
// // // // //             tif->tif_header_size = sizeof(TIFFHeaderBig);
// // // // //             tif->tif_flags |= TIFF_BIGTIFF;
// // // // //           }
// // // // //           tif->tif_flags |= TIFF_MYBUFFER;
// // // // //           tif->tif_rawcp = tif->tif_rawdata = 0;
// // // // //           tif->tif_rawdatasize = 0;
// // // // //           tif->tif_rawdataoff = 0;
// // // // //           tif->tif_rawdataloaded = 0;
// // // // //
// // // // //           switch (mode[0]) {
// // // // //             case 'r':
// // // // //               if (!(tif->tif_flags&TIFF_BIGTIFF))
// // // // //                 tif->tif_nextdiroff = tif->tif_header.classic.tiff_diroff;
// // // // //               else
// // // // //                 tif->tif_nextdiroff = tif->tif_header.big.tiff_diroff;
// // // // //               /*
// // // // //                * Try to use a memory-mapped file if the client
// // // // //                * has not explicitly suppressed usage with the
// // // // //                * 'm' flag in the open mode (see above).
// // // // //                */
// // // // //               if (tif->tif_flags & TIFF_MAPPED)
// // // // //               {
// // // // //                 toff_t n;
// // // // //                 if (TIFFMapFileContents(tif,(void**)(&tif->tif_base),&n))
// // // // //                 {
// // // // //                   tif->tif_size=(tmsize_t)n;
// // // // //                   assert((toff_t)tif->tif_size==n);
// // // // //                 }
// // // // //                 else
// // // // //                   tif->tif_flags &= ~TIFF_MAPPED;
// // // // //               }
// // // // //               /*
// // // // //                * Sometimes we do not want to read the first directory (for example,
// // // // //                * it may be broken) and want to proceed to other directories. I this
// // // // //                * case we use the TIFF_HEADERONLY flag to open file and return
// // // // //                * immediately after reading TIFF header.
// // // // //                */
// // // // //               if (tif->tif_flags & TIFF_HEADERONLY)
// // // // //                 return (tif);
// // // // //
// // // // //               /*
// // // // //                * Setup initial directory.
// // // // //                */
// // // // //               if (TIFFReadDirectory(tif)) {
// // // // //                 tif->tif_rawcc = (tmsize_t)-1;
// // // // //                 tif->tif_flags |= TIFF_BUFFERSETUP;
// // // // //                 return (tif);
// // // // //               }
// // // // //               break;
// // // // //             case 'a':
// // // // //               /*
// // // // //                * New directories are automatically append
// // // // //                * to the end of the directory chain when they
// // // // //                * are written out (see TIFFWriteDirectory).
// // // // //                */
// // // // //               if (!TIFFDefaultDirectory(tif))
// // // // //                 goto bad;
// // // // //               return (tif);
// // // // //           }
// // // // //           bad:
// // // // //           tif->tif_mode = O_RDONLY;       /* XXX avoid flush */
// // // // //           TIFFCleanup(tif);
// // // // //           bad2:
// // // // //           return ((TIFF*)0);
// // // // //       }






