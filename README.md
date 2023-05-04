# MetaAnalysis

## Code to calculate the models of the extracted effect sizes and create the figures. 

### Settings to be made by the user:
  - set the path to the data
  - choose the "correction_mode":
  
      Influences effect sizes calculated from dependen t-tests. To avoid overestimation the following formula was used: d=t √(2(1-r)/n) (Dunlap et al., 1996)
      - (1) chen: based on the data kindly provided by Chen et al. (2022), (r = .63, RS baseline vs. after neurofeedback: r = .92, block 1 vs block 10: r = .67))
      - (2) con: conservative -  based on the retest-reliability of cognitive EEG (McEvoy et al.,2000), (r = 0.97 )
      - (3) lib: liberal - based on the retest-reliability of cognitive EEG (McEvoy et al.,2000),(r = 0.74) 
  
  - Outliers and influential studies should be investigated manually before exclusion
  - marking of the influential outliers in the plot needs to be manually adjusted, if the size is not A4
