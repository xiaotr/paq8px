/*
  paq8px file compressor/archiver.

  Copyright (C) 2008 Matt Mahoney, Serge Osnach, Alexander Ratushnyak,
  Bill Pettis, Przemyslaw Skibinski, Matthew Fite, wowtiger, Andrew Paterson,
  Jan Ondrus, Andreas Morphis, Pavel L. Holoborodko, Kaido Orav, Simon Berger,
  Neill Corlett, Márcio Pais

  LICENSE

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details at
  Visit <http://www.gnu.org/copyleft/gpl.html>.
*/


#ifndef ENCODER_H_INCLUDED
#define ENCODER_H_INCLUDED

#include "misc.h"

typedef enum {COMPRESS, DECOMPRESS} Mode;

class Predictor {
  int pr;  // next prediction
public:
  Predictor();
  int p() const;
  void update();
};

class Encoder {
private:
  Predictor predictor;
  const Mode mode;       // Compress or decompress?
  FILE* archive;         // Compressed data file
  U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
  U32 x;                 // Decompress mode: last 4 input bytes of archive
  FILE *alt;             // decompress() source in COMPRESS mode
  float p1, p2;

  // Compress bit y or return decompressed bit
  int code(int i=0);

public:
  Encoder(Mode m, FILE* f, int lvl);
  Mode getMode() const {return mode;}
  long size() const {return ftell(archive);}  // length of archive so far
  void flush();  // call this when compression is finished
  void setFile(FILE* f) {alt=f;}
  void compress(int c); // Compress one byte
  int decompress(); // Decompress and return one byte

  void set_status_range(float perc1, float perc2) { p1=perc1; p2=perc2; }
  void print_status(int n, int size);
  void print_status();
};
#endif // #ifndef ENCODER_H_INCLUDED