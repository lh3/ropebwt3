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
