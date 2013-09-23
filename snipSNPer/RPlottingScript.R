require(ggplot2)
setwd(cwd)
cutsites = read.csv("cutsiteslist.csv")
positions = read.csv("positionlist.csv")
start = min(positions[,1])
sto = max(positions[,1])
for (i in seq(1,length(cutsites$Position))){
  cutsites$Position[i] = as.numeric(as.character(cutsites$Position[i]))
}
if (length(positions[,1]) == 2){
  plot = ggplot() + geom_rect(data = cutsites, mapping = aes(xmin = min(Position), xmax = max(Position), ymin = 0, ymax = 1), fill = "#545454") + geom_rect(mapping = aes(xmin = start, xmax = sto, ymin = 0, ymax = 1), fill = "#238E23") + ylim(c(-1,2)) + geom_segment(data = cutsites, aes(x = Position, y = 0, xend = Position, yend = 1, col = Enzyme), size = 1.5, show_guide = FALSE) + geom_point(data = cutsites[seq(1,nrow(cutsites),2),], aes(x = Position, y = 1, fill = Enzyme, col = Enzyme),  size = 7, shape = 25) + geom_point(data = cutsites[seq(2,nrow(cutsites),2),], aes(x = Position, y = 0, fill = Enzyme, col = Enzyme),  size = 7, shape = 24) + geom_text(data = cutsites[seq(1,nrow(cutsites),2),], aes(x = Position, y = 1.2, label = seq(1,nrow(cutsites),2), col = Enzyme), size = 10, show_guide=FALSE) + geom_text(data = cutsites[seq(2,nrow(cutsites),2),], aes(x = Position, y = -.2, label = seq(2,nrow(cutsites),2), col = Enzyme), size = 10, show_guide=FALSE) + theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y =  element_blank()) + xlab("Position")
}

if (length(positions[,1]) == 1){
  plot = ggplot() + geom_rect(data = cutsites, mapping = aes(xmin = min(Position), xmax = max(Position), ymin = 0, ymax = 1), fill = "#545454") + geom_segment(mapping = aes(x = positions[1,1], y = 0, xend = positions[1,1], yend = 1), col = "#238E23", size = 2) + ylim(c(-1,2)) + geom_segment(data = cutsites, aes(x = Position, y = 0, xend = Position, yend = 1, col = Enzyme), size = 1.5, show_guide = FALSE) + geom_point(data = cutsites[seq(1,9,2),], aes(x = Position, y = 1, fill = Enzyme, col = Enzyme),  size = 7, shape = 25) + geom_point(data = cutsites[seq(2,10,2),], aes(x = Position, y = 0, fill = Enzyme, col = Enzyme),  size = 7, shape = 24) + geom_text(data = cutsites[seq(1,9,2),], aes(x = Position, y = 1.2, label = seq(1,9,2), col = Enzyme), size = 10, show_guide=FALSE) + geom_text(data = cutsites[seq(2,10,2),], aes(x = Position, y = -.2, label = seq(2,10,2), col = Enzyme), size = 10, show_guide=FALSE) + theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.text.y =  element_blank()) + xlab("Position")
}


ggsave(filename = filename, plot = plot, width = 246.327083333, height = 137.31875, units = "mm")

ggsave(filename = "ggplot.png", plot = plot, width = 246.327083333, height = 137.31875, units = "mm")