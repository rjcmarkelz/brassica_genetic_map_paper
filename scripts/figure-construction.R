library(qtl)
library(ggplot2)
plotMap(newmap)

plotMap
newmap
str(genetic_map)

genetic_map <- shiftmap(genetic_map)

out <- scanone(genetic_map, pheno.col = 2)
out2 <- out
head(out2)
out2$len  <- 1
out2$xout <- 1
out2$pos
str(out2)


map_plot <- ggplot(out2, aes(len, pos)) + geom_point() + facet_grid(. ~ chr)
map_plot <- map_plot + geom_errorbarh(aes(xmax = xout + 5, xmin = xout - 5, height = 0))
map_plot <- map_plot + theme(panel.background = NULL)
map_plot


map_plot <- ggplot(out2, aes(len, pos)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + theme(panel.background = NULL)
map_plot

map_plot <- ggplot(out2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + scale_y_reverse()
map_plot <- map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
map_plot <- map_plot + ylab("Genetic Position (cM)") + ggtitle("")
map_plot


map_plot <- ggplot(out2, aes(len, pos, ymin = 0, ymax = 40)) + facet_grid(. ~ chr)
map_plot <- map_plot + geom_segment(aes(x = len - 5, y = pos, xend = len + 5, yend = pos))
map_plot <- map_plot + geom_segment(aes(x = len, y = 0, xend = len, yend = pos))
map_plot <- map_plot + scale_y_reverse()
map_plot <- map_plot + theme(panel.background = NULL, axis.text.x = element_blank(),
                             axis.line = element_blank(), panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), axis.ticks.x = element_blank(), 
                             axis.title.x = element_blank())
map_plot <- map_plot + ylab("Genetic Position (cM)") + ggtitle("")
map_plot

ggsave("test_genetic.pdf", map_plot)

