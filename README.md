# ShinyApp-Interactive-Complex-Heatmap
Interactive tool that visualises Spatial Omics Data in R

### Data: ###
In the cells.csv dataset, the columns represent the genes/markers and rows represent the patientID. The intersecting cell between the row and column is the gene expression value corresponding to that specific gene and patient. 

###  Approach: ###
* To derive the average expression value for each marker, for each patient and normalise it using Standard and Z-Score normalisation.
* Visualise the Original, Standard and Z-Score normalised data in Complex Heatmaps in R.
* Annotate the heatmaps with patientID and indication to expedite data analysis.
* Develop a ShinyApp for further user interaction and analysis of patient data, with dropdown options to select the type of available heatmap and turn on/off the column annotations as preferred. 

<img width="1426" alt="image" src="https://github.com/estherteo/ShinyApp-Interactive-Complex-Heatmap/assets/104299126/6aaf26bf-ad2b-4f6e-9940-b0de61d86efd">
