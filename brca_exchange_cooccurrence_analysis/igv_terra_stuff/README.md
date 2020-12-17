## Running instructions

1) Make a new workspace in Terra
2) Download the .ipynb from [here](https://github.com/cmarkello/Sandbox/blob/master/brca_exchange_cooccurrence_analysis/igv_terra_stuff/IGV_WDL_Terra_setup.ipynb).
3) Upload `IGV_WDL_Terra_setup.ipynb` to the `NOTEBOOKS` tab in the Terra workspace
4) Create a table.tsv tab-delimited file like the example provided [here](https://github.com/cmarkello/Sandbox/blob/master/brca_exchange_cooccurrence_analysis/igv_terra_stuff/table.tsv) and add your data in the following way for each row:
{TOPMED_SAMPLE_ID}\t{VARIANT_POSITION}\t{VARIANT_RANGE}
where:
  - TOPMED_SAMPLE_ID: The TOPMed ID of the sample. Starts with an NWD prefix
  - VARIANT_POSITION: The position and definition of the variant of interest. Formatted as {VARIANT_CHR}:{VARIANT_POS}:{VARIANT_REF}:{VARIANT_ALT}
  - VARIANT_RANGE: The range of the variant to render in igv. Formatted as {VARIANT_CHR}:{POSITION_START}-{POSITION_END}

5) Upload your modified table.tsv to the Files directory found at the bottom of the `DATA` tab directory listing in the Terra workspace.
6) Upload the igv batch script formatting script from [here](https://github.com/cmarkello/Sandbox/blob/master/brca_exchange_cooccurrence_analysis/igv_terra_stuff/make_igv_batchfile.py) to the Files directory found at the bottom of the `DATA` tab directory listing in the Terra workspace.
7) Run through the `IGV_WDL_Terra_setup.ipynb` notebook in the `NOTEBOOKS` tab in the Terra workspace.
8) Go to the `WORKFLOWS` tab in the Terra workspace and add a new workflow, searching for `cmarkello/Sandbox` in the Dockstore repo.
9) Go back to the `WORKFLOWS` tab in the Terra workspace, select the `Sandbox` workflow. Then select the Run workflow with inputs defined by file paths button on the upper-left side of the workflow setup window. Then fill in the inputs list as the following:
- `batch_script` : Click on the folder icon on the right and select the file `make_igv_batchfile.py` from the dropdown window.
- `input_table` : Click on the folder icon on the right and select the file `igv_table.tsv` from the dropdown window.
- `output_folder_prefix` : "test_igv_output"

10) Click the `SAVE` button on the right side.
11) Click `RUN ANALYSIS` button.
12) Go to the `JOB HISTORY` tab in the Terra workspace and wait for the submitted job to complete.
13) Once completed it should have a `Status` of `Done` with a green checkmark if it ran to completion, otherwise there'll be a red icon.
14) Click on the `Submission job` for details.
15) Click on the folder icon on the right end of the table. That should redirect you to a google storage bucket directory.
16) Go to the `call-gather_shards` directory.
17) Download the file `test_igv_output.tar.gz` , unzip it and it should contain the screenshot files.
