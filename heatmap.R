# yorig: labels
# ypred: labels to compare with
# filename: if you want to save pdf (optional)
# title:  title of the plot (optional)

createheatmap <- function(yorig, ypred, filename = 0, title = 0){
  x_table = table(yorig, ypred)
  x_table2 = x_table / rowSums(x_table) # Normalize the counts such that they are between 0 and 1 (to color the heatmap)
  x_table2 = data.frame(x_table2)
  x_table = data.frame(x_table)
  x_table2$Counts = x_table$Freq
  grey_col = grey.colors(5000, start = 1, end = 0) # Color of the text in the heatmap
  grey_col[1:2500] = grey_col[1] #low values are white
  grey_col[2501:5000] = grey_col[5000] #high values are black
  
  if (filename != 0){
    pdf(filename, width = 10, height = 9.5)
  }
  
  # Order of the columns 
  order_y <- rev(c('pbmc', 'Myeloid', 'DC', 'cDC2', 'ASDC', 'pDC', 'Monocyte', 'CD14 Mono', 'CD16 Mono',
                   'Lymphoid', 'B cell', 'B intermediate', 'B naive', 'B memory', 'Plasmablast', 'NK cell', 'NK', 'NK_CD56bright',
                   'NK Proliferating', 'T cell', 'CD8+ T cell', 'CD8 TCM', 'CD8 Naive', 'CD8 TEM', 'CD8 Proliferating',
                   'CD4+ T cell', 'Treg', 'CD4 Proliferating', 'CD4 Naive', 'CD4 TCM', 'CD4 CTL', 'CD4 TEM', 'gdT', 'MAIT', 'dnT',
                   'HSPC', 'Platelet', 'Eryth'))
  
  # Order of the rows
  order_x <- c('cDC2', 'ASDC', 'CD14 Mono', 'CD16 Mono', 'B intermediate', 'B naive', 'B memory', 
               'Plasmablast', 'NK', 'NK_CD56bright', 'NK Proliferating', 'CD8 TCM', 'CD8 Naive', 'CD8 TEM', 'Treg',
               'CD4 Proliferating', 'CD4 Naive', 'CD4 TCM', 'CD4 CTL', 'CD4 TEM', 'gdT', 
               'MAIT', 'dnT', 'HSPC', 'Platelet', 'Eryth')
  
  # Colormap
  hmcol<-viridis(500)
  grey_col[1] <- hmcol[1] #Plot 0 as invisible 
  p3 = ggplot(data = x_table2, aes(x = ypred, y = yorig, fill = Freq)) +
    geom_tile(colour="white",size=0.25) + 
    geom_text(aes(label = Counts, colour = Freq), size = 4) +
    scale_colour_gradientn(colours=grey_col, limits=c(0, 1.1), guide = FALSE) +
    scale_fill_gradientn(colours=hmcol, limits=c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1),
                         guide = FALSE) +
    scale_y_discrete(expand=c(0,0), limits = order_y) + 
    scale_x_discrete(expand=c(0,0), limits = order_x) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          strip.text.x = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          plot.background=element_blank(),
          panel.border=element_blank(),
          legend.position = 'bottom',
          legend.direction = 'horizontal') +
    #     coord_equal() +
    xlab('Azimuth') +
    ylab('HierscPred')
  
  if(title != 0){
    p3 = p3 + ggtitle(title)
  }
  
  print(p3)
  
  if (filename != 0){
    dev.off()
  }
  
}
