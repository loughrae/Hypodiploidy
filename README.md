# Hypodiploidy
Scripts for paper on karyotypic evolution of hypodiploid tumours

**1. Processing Karyotypes **

All samples in the Mitelman Database on 25/09/2023 with morphology 'acute lymphoblastic leukaemia/lymphoblastic lymphoma' were downloaded via Google BigQuery. In the script 'run_cytoconverter_mitALL.R', long-form karyotypes were extracted and run through CytoConverter, which produces a table of copy number events (gains and losses). This takes ~30 minutes.



