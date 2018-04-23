library(ggplot2)
library(RColorBrewer)

csv_files = grep(".csv", value = T, dir())

for( i in csv_files ){

        table = read.table(i, sep = ",", header = T,row.names = NULL)

        trimmed_by_pvalues = table[table$p.value <= 0.5,]

        samples = unique(trimmed_by_pvalues$Sample)

        png(gsub("csv", "png", i), width = 9, height = 6, units = 'in', res = 900)
        plots <- ggplot(trimmed_by_pvalues,
                        aes(x = Bases_from_end, y = R.squared, colour = Sample))+
                scale_color_manual(values = colorRampPalette(brewer.pal(9,"Blues"))(length(samples)))+
                geom_point(size = 1.3)+
                #stat_smooth(se = T, size = 2, span = 0.35)+
                coord_cartesian(xlim = table$Bases_from_end, ylim = c(0,1))+
                geom_line(size = 1.1)+
                theme_bw(base_size = 15)+
                labs(y = "R-squared",
                     x = "Length of adapter",
                     title = paste("Matches of",
                                   gsub(".csv", "", i) ,
                                   "adapter in all samples"))+
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.position = "none",
                      plot.title = element_text(size=19, face= "bold"))
        print(plots)
        dev.off()
}
