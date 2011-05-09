/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2010 Jonas Behr, Regina Bohnert, Gunnar Raetsch
*   Copyright (C) 2010 Max Planck Society
*/


#ifndef __READ_H__
#define __READ_H__

#include <stdint.h>
#include <stdio.h>
#include <vector>
  using std::vector;


class CRead {   
 public:
  /** constructor
   */
  CRead();
  ~CRead();
  
  vector<int> block_starts;
  vector<int> block_lengths;
  char* read_id;
  char* sam_line;
  int start_pos;
  char * strand;
  int matches;
  int mismatches;
  int multiple_alignment_index;
  
  void get_coverage(int p_start_pos, int p_end_pos, uint32_t* coverage);
  void get_reads_sparse(int p_start_pos, int p_end_pos, double* reads, uint32_t & reads_c, uint32_t row_idx);
  void get_introns(vector<int>* introns);
  void get_acc_splice_sites(vector<int>* acc_pos);
  void get_don_splice_sites(vector<int>* acc_pos);
  int max_intron_len();
  int min_exon_len();
  bool operator==(const CRead& read) const;
  void print();
  void set_strand(char s);
  int get_mismatches();
};
#endif
