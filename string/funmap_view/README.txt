The STRING view (the five files listed in SHA256SUMS, 66 MB) is rebuilt by  python code/build_string_from_shipped_nb.py
from string/string_top50_neighborhoods.csv.gz, string/string_edges.csv.gz and the frozen GO tables of funmap/; it is also
shipped in the Zenodo archive of this release.  Every neighbourhood-family script reads it through GOCO_DATA=string/funmap_view.
