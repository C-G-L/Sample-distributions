##############################################################
# load dependencies
##############################################################
library(ggplot2)
library(maps)
library(mapdata)
library(plyr)
library(alphahull)
library(G1DBN)

##############################################################
# distribution function
##############################################################
distribution <- function(data,
                         info,
                         title = NULL,
                         xlim = NULL,
                         ylim = NULL,
                         axes = FALSE,
                         border = "white",
                         bg = "white",
                         polygon = TRUE,
                         points = TRUE, 
                         alphaFill = 0.5,
                         alphaPoints = 0.5,
                         alphaShape = 1){
  
  map("worldHires", xlim = xlim, ylim = ylim, col="gray70", border = border, fill = TRUE, bg = bg, lforce="e")
  
  title(main = title, line = 1)
  mtext('longitude', side = 1, line = 2.5, cex = 0.8)
  mtext('latitude', side = 2, line = 2.5, cex = 0.8)
  
  if(axes){map.axes()}
  
  e <- 1
  alphaS <- alphaShape
  
  while(e <= nrow(info)){
    
    Loc2 <- merge(data, info, by = 'taxon_nombre')
    
    sub <- subset(Loc2, Loc2$order == e)

    LocUnique <- unique(sub[,c("lon_dec", "lat_dec")])
    
    colS <- as.vector(col2rgb(unique(sub$color)))/255
    
    if(polygon){
      
      n10 <- ashape(LocUnique, alpha = alphaS)
      n10g = graph.edgelist(cbind(as.character(n10$edges[, "ind1"]), as.character(n10$edges[, "ind2"])), 
                            directed = FALSE)
    
      if(!is.connected(n10g) || (clusters(n10g)$no > 1)) {  # any(degree(n10g) != 2) ||
        
        alphaS <- alphaS + 0.1
        print(paste('increasing alpha with 0.1 for', info$taxon_nombre[e], sep = ' '))
      }
      
      else{
        
        cutg = n10g - E(n10g)[1]
        # find chain end points
        ends = names(which(degree(cutg) == 1))
        path <- tryCatch(get.shortest.paths(cutg, ends[1], ends[2])[[1]],
                         warning = function(err) {
                           print(paste('increasing alpha', info$taxon_nombre[e], sep = ' '))
                           return(NULL)
                           })
        
        if(is.null(path)){
          alphaS <- alphaS + 0.1
          next
          }
        
        # this is an index into the points
        pathX = as.numeric(V(n10g)[path[[1]]]$name)
        # join the ends
        pathX = c(pathX, pathX[1])
        
        ashapem <- as.matrix(n10$x[pathX, ])
        
        ashapem[,2][ashapem[,2] < ylim[1]] <- ylim[1]
        ashapem[,2][ashapem[,2] > ylim[2]] <- ylim[2]
        ashapem[,1][ashapem[,1] < xlim[1]] <- xlim[1]
        ashapem[,1][ashapem[,1] > xlim[2]] <- xlim[2]
        
        polygon(ashapem, col = rgb(colS[1], colS[2], colS[3], alpha = alphaFill), border = NA)
        
        
        if(points){
          
          points(x = LocUnique$lon_dec, y = LocUnique$lat_dec, col = rgb(colS[1], colS[2], colS[3], alpha = alphaPoints), 
                 pch = 19, cex = 1)
        }
        
        e <- e + 1
        alphaS <- alphaShape
      }
    }

    if(points == TRUE && polygon == FALSE){

      points(x = LocUnique$lon_dec, y = LocUnique$lat_dec, col = rgb(colS[1], colS[2], colS[3], alpha = alphaPoints), 
             pch = 17, cex = 1)
    
      e <- e + 1
    }

  }
}

##############################################################
# Run example for H. erato data
##############################################################

# Load locality files
files <- list.files(path='Localities/H_erato', pattern='.csv', full.names = T)

# Combine files
Loc <- c()
for(f in 1:length(files)){
  table <- read.csv(files[f], h=T)
  Loc <-rbind(Loc,table)
}

# Remove NA records
Loc <- Loc[!is.na(Loc$lon_dec),]
Loc <- Loc[!is.na(Loc$lat_dec),]

# Extract species/race names
LocSpec <- unique(Loc[,c("taxon_nombre")])

# Define color vector and transform to RGB colors (should have the same length as number of populations)
colVec <- c("#c851b5", "#5cb248", "#7f66d2", "#b6b444", "#6984c8", 
            "#80da8b", "#4cbcd1", "#cd4e34", "#54a978", "#c54b6c", 
            "#737a32", "#ba71a9", "#c27e5a", "#c851b5", "#000000")

# Show color scheme
plot(0, xlim=c(0,15), ylim=c(0,1), pch = '')
for(e in 1:length(colVec)){
  rect(e-1, 0, e, 1, col = colVec[e])
}

# Define order for plotting (should have the same length as number of populations)
order <- c(5,3,7,13,2,1,4,6,12,11,10,9,8,14,15)

# Combine colors and species/races
Spec <- as.data.frame(cbind(as.character(LocSpec), colVec, order))
colnames(Spec) <- c('taxon_nombre', 'color', 'order')

# Plot distribution
layout(matrix(c(1:3), nrow=1, byrow=TRUE))
layout.show(n=3)

par(mar=c(1,1,1,1), oma=c(1,2,1,1))

distribution(data = Loc, info = Spec, axes = TRUE, 
             title = substitute(paste(italic("H. erato "), "points ", sep=" ")),
             xlim = c(-90,-35), ylim = c(-35, 15), 
             polygon = FALSE, points = TRUE, alphaPoints = 0.5)

distribution(data = Loc, info = Spec, axes = TRUE,
             title = substitute(paste(italic("H. erato "), "alpha Hull = 100", sep=" ")),
             xlim = c(-90,-35), ylim = c(-35, 15), 
             polygon = TRUE, points = TRUE, alphaFill = 0.5, alphaShape = 100)

distribution(data = Loc, info = Spec, axes = TRUE,
             title = substitute(paste(italic("H. erato "), "alpha Hull = 1", sep=" ")),
             xlim = c(-90,-35), ylim = c(-35, 15), 
             polygon = TRUE, points = TRUE, alphaFill = 0.5, alphaShape = 4)