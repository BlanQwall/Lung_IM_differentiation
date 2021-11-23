barChart <- function(freqTable.list, 
                     sample.names = names(freqTable.list) 
                     # showTop = nrow(freqTable), 
                     # cols = NULL, 
                     # labels = TRUE, 
                     # num.labels = showTop, 
                     # title = NULL, 
                     # cex = 1, 
                     # show.pct = TRUE, 
                     # pct.inside = TRUE
                     ){
  # Load ggplot2
  require(ggplot2)
  require(RColorBrewer)
  
  if (!class(freqTable.list) == "list") {stop("freqTable.list must be a list of freqTable.")}
  if (is.null(sample.names)) {sample.names <- 1:length(freqTable.list)}
  if (! length(sample.names) == length(freqTable.list)) {stop("Sample.names is not equal to the number of freqTable in list.")}
  
# Create Data
data <- data.frame()

for (i in 1:length(freqTable.list)) {
  freqTable.list[[i]]$sample <- rep(sample.names[i], nrow(freqTable.list[[i]]))
  data <- rbind(data, freqTable.list[[i]])
}


# data[[1]] <- as.character(data[[1]])

# show the rest: 
# if (showTop < lgth) {
#   data[showTop+1, 2] <- sum(data[(showTop+1):lgth,2])
#   data[showTop+1, 1] <- "Others"
#   data <- data[1:(showTop+1),]
#  rownames(data)[showTop+1] <- "Others"
# }
data$sample <- factor(data$sample, levels = sample.names)
#levels(data[[1]]) <- as.character(data[[1]])

#  return(data)
ggplot(data = data, aes_string(x = colnames(data)[3], y = colnames(data)[2], fill = colnames(data)[1] )) + 
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw()




# create palette: 
# if (is.null(cols)){
#   palette.RColorBrewer <- "PuRd"
# if (showTop < 3) {pal <- brewer.pal(3, palette.RColorBrewer)[3:1] 
# pal <- pal[showTop:1]}
# else{
#   pal <- brewer.pal(ifelse(showTop>9, 9, showTop), palette.RColorBrewer)
#   if (showTop > 9) {pal <- colorRampPalette(pal)(showTop)}
#   } 
# pal <- pal[showTop:1]} else 
# {pal <- cols}
# 
# if (showTop < lgth) {pal <- c(pal,"#e6e6e6")} # add color for other
# 
# 
# 
# names(pal) <- rownames(data)
# if (labels) { 
#   labels <- rep("", nrow(data))
#   labels[1:num.labels] <- data[1:num.labels,1]
#   if (num.labels == showTop) {labels <- data[,1]} 
#   } else {labels <- NULL}
# 
#  print(labels)
# 
# if (show.pct) {
#   prop <- data[,2]/sum(data[,2])
#   pct <- signif(prop*100, 3)
#   pct <- paste(pct, "%", sep = "")
#   pct[(labels == "")] <- ""
#   if (!pct.inside) {
#     pct <- paste("(", pct, ")", sep = "")
#     labels <- paste(labels, pct)
#     labels[labels == " ()"] <- ""
#   }
#   
# #  print(pct) 
# }
# 
# # pie chart: 
# # ----
# pie(data[[2]],
#     labels = labels,
#     border="white",
#     col=pal,
#     init.angle = 90,
#     clockwise = TRUE, 
#     main = title, 
#     cex = cex)
# # -------
# 
# #print(prop)
# # add percentages text: 
# if(pct.inside)
# { radi = 0.6
#   prop.accumu <- numeric(length(labels))
#   
#   for (i in 1:length(labels)){
#     prop.accumu[i] <- sum(prop[1:i])
#   }
#   
#   #cat("prop.accu \n")
#   #print(prop.accumu)
#   
#   prop.accumu <- prop.accumu - 0.5*prop
#   x = sin(prop.accumu*2*pi)*radi
#   y = cos((prop.accumu)*2*pi)*radi
#   text(x, y, pct, cex = 0.7)}

# Basic piechart
# ggplot(data, aes(x="", y=Freq, fill=clonotype_id)) +
#   geom_bar(stat="identity", width=1, color="white") +
#   coord_polar("y", start=0) + 
#   theme_void() + # remove background, grid, numeric labels
#   scale_fill_manual(values = pal)
}
