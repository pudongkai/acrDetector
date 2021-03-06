#### AcrDetector

- Requirement

  In order to run this pragram, follow executables and python libs are required:

  `python3.*`
  `bioPython`
  `numpy`
  `json`
  `pandas`
  `joblib`
  
- Prepare

  To detect Acr, genomic island and prophage in genome should be detected. You could predict  genomic island by IslandViewer(http://www.pathogenomics.sfu.ca/islandviewer/) and prophage by  PHASTER(https://phaster.ca/). And GenBank annotation file is need.

  ```
  eg:
  	NZ_CP011374.1.gbff (GenBank annotation file)
  	NZ_CP011374.1.tsv (IslandViewer result file)
  	NZ_CP011374.1.phage (PHASTER result file)
  ```

- Usage

  ```shell
  cd acrDetector/
  python3 main.py basePath/NZ_CP011374.1
  ```
  **Note:** The basePath should be absolute Path.

- Result

  Result file will save in the base path. 
