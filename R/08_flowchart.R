library("DiagrammeRsvg")
library("rsvg") 
library("DiagrammeR")

flowchart <- grViz("
digraph flowchart {

  # a 'graph' statement
  graph [overlap = true, fontsize = 10]


  ####### ####### ####### ####### ####### ####### NODES  ####### ####### ####### #######
  ####### RAW 'NODE' STATEMENTS ####### 
  node [shape = box,
        fontname = Arial,
        fillcolor = '#1B325F',
        fontsize = 25,
        style = filled,
        fontcolor = white];
  nodesep= 1.25;
  
  'GSE13937_series_matrix.txt'; 'A-GEOD-8835.adf.txt'; 'A-GEOD-8835_comments.txt';
  
  
  ####### LOAD 'NODE' STATEMENTS  #######
  #CSV #FFF6CD
  node [shape = box,
        fontname = Arial,
        fillcolor = '#E9F2F9', 
        style = filled,
        fontcolor = black]
  
  'data_load.csv' 'probes_data_load.csv'
  
  #NO CSV
   node [shape = box,
        fontname = Arial,
        fillcolor = white,
        style = filled]
  
  raw_data ; 'raw_probes'; 'raw_comments'; 
  
  
  
  ####### CLEAN 'NODE' STATEMENTS  #######
  node [shape = box,
        fontname = Arial,
        fillcolor = '#9CC4E4',
        style = filled]
  
  'data_clean.csv' 
  'probes_data_clean.csv'
  
  ######## AUGMENT 'NODE' STATEMENTS  #######
  #CSV
  node [shape = box,
        fontname = Arial,
        fillcolor = '#3A89C9',
        stile = filled]
  
  'data_aug_wide.csv'
  'data_aug_long.csv'
  
  #NO CSV
  node [shape = box,
        fontname = Arial,
        fillcolor = white,
        style = filled]
  data_aug
  
  ######## RESULTS NODE  #######
  node [shape = box,
        fontname = Arial,
        style = filled, 
        fillcolor='#FDBBA8']
  'Results and visualizations' [fillcolor = '#F26C4F']
  
  'PCA'; 'Differential expression'; 'Kaplan-Meier'; 'Linear regression'
  
  
   
   ####### ####### ############### 'EDGE' STATEMENTS  ####### ####### ####### #######
  #Raw - load
  'GSE13937_series_matrix.txt' -> raw_data [label='', fontsize=15, fontname = Arial, fontcolor='#131A37']
  
  raw_data -> 'data_load.csv' 
  'A-GEOD-8835.adf.txt' -> 'raw_probes' 
  'A-GEOD-8835_comments.txt' -> 'raw_comments'
  {'raw_probes' 'raw_comments'} -> 'probes_data_load.csv'
  
  #Load - clean
  'data_load.csv' -> 'data_clean.csv'
  'probes_data_load.csv' -> 'probes_data_clean.csv'
  
  #Clean - augmented
  'data_clean.csv' -> data_aug
  {data_aug 'probes_data_clean.csv'} -> 'data_aug_wide.csv'
  'data_aug_long.csv' -> 'data_aug_wide.csv'
  'data_aug_wide.csv' -> 'data_aug_long.csv'


 #Augmented - Visualizations 
{'data_aug_wide.csv'} -> 'Results and visualizations' -> {'PCA' 'Differential expression' 'Kaplan-Meier' 'Linear regression'} 


  
{ rank = same; 'GSE13937_series_matrix.txt'; 'A-GEOD-8835.adf.txt'; 'A-GEOD-8835_comments.txt' }
{ rank = same; 'data_clean.csv'; 'probes_data_clean.csv' }
{ rank = same; 'data_aug_wide.csv'; 'data_aug_long.csv' }

}
")

flowchart = export_svg(flowchart)
flowchart = charToRaw(flowchart)
previous_wd = getwd()
setwd("./doc/images")
rsvg_png(flowchart, "flowchart.png", height = 900, width = 1500)
setwd(previous_wd)