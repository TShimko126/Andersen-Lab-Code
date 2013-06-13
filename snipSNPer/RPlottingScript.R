require(ggplot2)
setwd(cwd)
print(file.exists("cutsiteslist.csv"))
print(file.exists("positionlist.csv"))
cutsites = read.csv("cutsiteslist.csv")
positions = read.csv("positionlist.csv")
start = min(positions[,1])
sto = max(positions[,1])
for (i in seq(1,length(cutsites$Position))){
  cutsites$Position[i] = as.numeric(as.character(cutsites$Position[i]))
}
plot = ggplot() + geom_rect(data = cutsites, mapping = aes(xmin = min(Position), xmax = max(Position), ymin = 0, ymax = 1), fill = "#545454") + geom_rect(mapping = aes(xmin = start, xmax = sto, ymin = 0, ymax = 1), fill = "#238E23") + ylim(c(-1,2)) + geom_segment(data = cutsites, aes(x = Position, y = 0, xend = Position, yend = 1, col = Enzyme), size = 1.5) + geom_point(data = cutsites, aes(x = Position, y = 1, fill = Enzyme, col = Enzyme),  size = 7, shape = 25) + theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y =  element_blank())  
ggsave(filename = filename, plot = plot, width = 246.327083333, height = 137.31875, units = "mm")