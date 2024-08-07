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
