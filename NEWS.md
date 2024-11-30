Release 3.9-r259 (30 November 2024)
-----------------------------------

Notable changes:

 * New feature: with option -p, the build command optionally runs partial
   suffix array construction and BWT merge at the same time.

(3.9 30 November 2024, r259)



Release 3.8-r249 (18 October 2024)
----------------------------------

Notable changes:

 * Improvement: a faster algorithm to locate a subset of positions in a suffix
   array interval. Given m highly similar genomes, the expected time is O(s/m)
   where s is the suffix array sample rate.

 * New feature: added option `mem -p` to optionally output a semi-random subset
   of SMEM positions.

 * New feature: added auxiliary script `rb3tools.js` for generating mappability
   filter and for simple SNP calling.

(3.8: 18 October 2024, r249)



Release 3.7-r226 (17 September 2024)
------------------------------------

In this release, the `sw` command outputs the `cs` tag and can report all
end-to-end hits in a compact format with option `--all-e2e`.

(3.7: 17 September 2024, r226)



Release 3.6-r217 (10 September 2024)
------------------------------------

Notable changes:

 * Bugfix: fixed a rare assertion failure.

 * Improvement: the `stat` command now works with the FMR format (#1).

 * Breaking: the `hapdiv` command now gives the counts of alleles up to edit
   distance five.

(3.6: 10 September 2024, r217)



Release 3.5-r203 (31 August 2024)
---------------------------------

Notable changes:

 * New feature: added the end-to-end mode (`-e`) to `sw`. It outputs local
   haplotypes similar to the entire query string.

 * New feature: `hapdiv` command to estimate local haplotype diversity. It
   applies the end-to-end mode to slide 101-mers and reports the number of
   haplotypes within certain edit distance threshold.

 * Improvement: ~30% faster BWT construction for large datasets.

 * Breaking: in `sw`, changed the default scoring to the BLASTN scoring

(3.5: 31 August 2024, r203)



Release 3.4-r167 (20 August 2024)
---------------------------------

Notable changes:

 * Improvement: in `sw`, removed out-of-band cells earlier. This speeds up
   BWA-SW by 20%.

 * Improvement: `sw` now computes positions in each thread. This improves the
   multi-threading performance for short queries.

 * Bugfix: with `mem`, ambiguous bases in query caused segmentation fault.
   Ropebwt3 now converts ambiguous bases to "A". This does not affect `sw`.

 * Breaking: removed greedy MEM

(3.4: 20 August 2024, r167)



Release 3.3-r149 (6 August 2024)
--------------------------------

Notable changes in the `sw` command:

 * New feature: option to try SW only when there is a long MEM. This is not
   enabled by default.

 * New feature: option to output unmapped reads in PAF

 * Bugfix: backtracking the F state could be wrong in corner cases

 * Bugfix: coordinates on the reverse strand were not flipped in PAF

 * Breaking: don't output the reference sequence in the rs tag by default

Other new features:

 * New feature: added the stat command to report the number of runs. Only
   working for the FMD format for now.

 * New feature: added `--min-gap` to the `mem` command to output regions not
   covered by long MEMs.

(3.3: 6 August 2024, r149)



Release 3.2-r137 (23 July 2024)
-------------------------------

This release implemented several critical features for sequence search:

 * New feature: Travis Gagie's algorithm for finding long MEMs. It is faster
   and now the default algorithm for MEM finding.

 * New feature: BWA-SW for local alignment. This algorithm allows mismatches
   and short gaps.

 * New feature: sampled suffix array for obtaining mapping positions.

 * Breaking: renamed `match` to `mem`.

(3.2: 23 July 2024, r137)



Release 3.1-r77 (15 June 2024)
------------------------------

Notable changes:

 * New feature: the `match` command supports multi-threading.

 * Improvement: generating suffix arrays (SA) with libsais16x64, a new addition
   to libsais by @IlyaGrebnov. libsais16x64 is faster and uses less memory for
   multi-string SA construction than libsais64. This reduces the peak memory of
   ropebwt3 for human genomes.

(3.1: 15 June 2024, r77)



Release 3.0-r67 (11 June 2024)
------------------------------

First release.
