Tools
=====

All scripts are located in ./workflow/scripts.

Step 1
------
* ``count_junctions.py`` - extracts splice junctions from alignment file (BAM) with their read counts.

* ``count_sites.py`` - extracts splice sites with only continuous reads support from alignment file (BAM).

Step 2
------

* ``aggregate_junctions.py`` - aggregates read counts for junctions from the first step (J1).
  Adds strand information to junctions and filters out junctions
  if their corresponding introns are too short or too long.

* ``aggregate_sites.py`` - aggregates read counts for sites from the first step (S1/PS1).

Step 3
------

* ``annotate_junctions.py`` - retrieves splice site sequences for aggregated junctions from the second step (J2).
  Adds their annotation status.

Step 4
------

* ``choose_strand.py`` - chooses correct strand for each junction from the third step (J3).

Step 6
------

* ``filter.py`` - filters junctions or sites from previous steps by various criteria.

Step 7
------

* ``compute_rates.py``

Auxiliary
---------

* ``describe_replicates.py``

* ``merge_junctions.py``


