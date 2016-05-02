/*
 * cigargen.cc
 *
 *  Created on: Aug 6, 2014
 *      Author: Ivan Sovic
 */

//  The MIT License (MIT)
//
//  Copyright (c) 2014 Ivan Sovic
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.



#include "cigargen.h"

int32_t GenerateCigar(char *query, uint32_t query_length, char *reference, uint32_t reference_length, std::string *ret_cigar, uint32_t *ret_alignment_length, std::string *ret_alignment) {
  if (query == NULL || reference == NULL || query_length == 0 || reference_length == 0)
    return -1;

  // Preallocate the memory for the DP matrix and traceback, and initialize all values to 0.
  std::vector<std::vector<int32_t> > dp_matrix((query_length + 1), std::vector<int32_t>((reference_length + 1), 0));
  std::vector<std::vector<uint8_t> > dp_traceback((query_length + 1), std::vector<uint8_t>((reference_length + 1), 0));
  // Temporary variables for storing the options for the current step. (Could have been done differently, but increases readability.)
  int32_t up=0, left=0, diagonal=0;

  // Initialize the first row and column according to Needleman-Wunsch.
  for (uint32_t i=0; i<=query_length; i++)
    dp_matrix[i][0] = i * PENALTY_INSERTION;
  for (uint32_t i=0; i<=reference_length; i++)
    dp_matrix[0][i] = i * PENALTY_DELETION;

  int32_t max_score_index = 0;

  for (uint32_t i=1; i<=query_length; i++) {
    for (uint32_t j=1; j<=reference_length; j++) {
      // Calculate the options for the current value of the DP matrix.
      up = dp_matrix[i-1][j] + PENALTY_INSERTION;
      left = dp_matrix[i][j-1] + PENALTY_DELETION;
      diagonal = dp_matrix[i-1][j-1] + ((query[i-1] == reference[j-1]) ? PENALTY_MATCH : PENALTY_MISMATCH);

      // Find the maximum of the three values.
      dp_matrix[i][j] = diagonal;
      dp_traceback[i][j] = ((query[i-1] == reference[j-1])?0:1);
      if (up > dp_matrix[i][j]) {
        dp_matrix[i][j] = up;
        dp_traceback[i][j] = 2;
      }
      if (left > dp_matrix[i][j]) {
        dp_matrix[i][j] = left;
        dp_traceback[i][j] = 3;
      }

      // If this is the last row of the matrix, we need to look up for the maximum score. This will be our starting point for the
      // CIGAR traceback calculation.
      if (i == query_length && dp_matrix[query_length][j] > dp_matrix[query_length][max_score_index])
        max_score_index = j;
    }
  }

  // Initialize the traceback column to the index of the highest score, and the traceback row to the last row of the matrix.
  int32_t current_traceback_row = query_length, current_traceback_column = max_score_index;

  std::stringstream cigar_stream;
  std::string cigar = "";
  std::string alignment = "";
  uint32_t edit_distance = 0;

  // Traceback: Find the CIGAR string.
  uint32_t iterations=0, num_insertions=0, type=0, last_type=0, num_same_types=0;
  while (current_traceback_column > 0 && current_traceback_row > 0) {
    type = dp_traceback[current_traceback_row][current_traceback_column];

    // Count the number of same transitions in the DP matrix (whether if it was to the left, right or up.)
    if (iterations == 0) {
      // If this is the first traceback step, initialize the last type for counting.
      last_type = type;
      num_same_types = 1;
    }
    else {
      if (type == last_type) {
        // If types of transition are the same, just count them.
        num_same_types += 1;
      }
      else {
        // If the type has changed, stream the CIGAR count, and reset the type to the current one.
        cigar_stream << TypeToSymbol(last_type) << num_same_types;
        last_type = type;
        num_same_types = 1;
      }
    }

    if (type != TYPE_MATCH)
      edit_distance += 1;

    // Accumulate the characters for the alignment, but only if this was asked for (that's why there is
    // the NULL testing). We insert a special character '+' in this implementation, which represents that
    // there was an insertion on the reference (one character was removed from the query).
    // Character '-' represents the deletion in the reference (or an insertion on the query).
    if (ret_alignment != NULL)
        alignment +=  ((type == TYPE_DELETION) ? '-' :
                      ((type == TYPE_MATCH) ? query[current_traceback_row - 1] :
                      ((type == TYPE_INSERTION) ? '+' :
                      'x')));

    if (type == TYPE_INSERTION)
      num_insertions += 1;

    // Row changes in all cases except when a deletion occurs (type == 3).
    current_traceback_row -= (type != TYPE_DELETION) ? 1 : 0;
    // Column changes in all cases except when an insertion occurs (type == 2).
    current_traceback_column -= (type != TYPE_INSERTION) ? 1 : 0;

    // Count the traceback length.
    iterations += 1;
  }

  // Stream the remaining count to the CIGAR stream.
  cigar_stream << TypeToSymbol(last_type) << num_same_types;
  cigar = cigar_stream.str();
  std::reverse(cigar.begin(), cigar.end());
  std::reverse(alignment.begin(), alignment.end());

  if (ret_cigar != NULL)
    *ret_cigar = cigar;
  else
    Verbose(stdout, query, query_length, reference, reference_length, cigar, dp_matrix, dp_traceback);

  // Return the final alignment length. We need to subtract the number of insertions on the reference because the final
  // alignment should be shorter than the traceback for this amount.
  if (ret_alignment_length != NULL)
    *ret_alignment_length = (iterations - num_insertions);

  // Return the final alignment if required. The alignment.size() value is larger than ret_alignment_length, because
  // query deletions have not been completely removed, but changed with a special character.
  if (ret_alignment != NULL)
    *ret_alignment = alignment;

  return edit_distance;
}

void Verbose(FILE *fp, char *query, uint32_t query_length, char *reference, uint32_t reference_length, std::string &cigar, std::vector<std::vector<int32_t> > &dp_matrix, std::vector<std::vector<uint8_t> > &dp_traceback) {
  fprintf (fp, "%3c", ' ');
  for (uint32_t i=0; i<=reference_length; i++) {
    if (i > 0)
      fprintf (fp, "%3c", reference[i-1]);
    else
      fprintf (fp, "%3c", '-');
  }
  fprintf (fp, "\n");

  for (uint32_t i=0; i<=query_length; i++) {
    if (i > 0)
      fprintf (fp, "%3c", query[i-1]);
    else
      fprintf (fp, "%3c", '-');

    for (uint32_t j=0; j<=reference_length; j++) {
      fprintf (fp, "%3d", dp_matrix[i][j]);
    }
    fprintf (fp, "\n");
  }

  fprintf (fp, "\n");
  fprintf (fp, "Traceback:\n");

  fprintf (fp, "%3c", ' ');
  for (uint32_t i=0; i<=reference_length; i++) {
    if (i > 0)
      fprintf (fp, "%3c", reference[i-1]);
    else
      fprintf (fp, "%3c", '-');
  }
  fprintf (fp, "\n");

  for (uint32_t i=0; i<=query_length; i++) {
    if (i > 0)
      fprintf (fp, "%3c", query[i-1]);
    else
      fprintf (fp, "%3c", '-');

    for (uint32_t j=0; j<=reference_length; j++) {
      fprintf (fp, "%3d", dp_traceback[i][j]);
    }
    fprintf (fp, "\n");
  }

  fprintf (fp, "\n");

  fprintf (fp, "CIGAR: %s\n", cigar.c_str());
  fprintf (fp, "\n");
}
