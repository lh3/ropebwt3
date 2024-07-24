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
