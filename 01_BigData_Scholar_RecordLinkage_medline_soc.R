library(tibble)
library(reshape2)
library(ggplot2)



## Code and data plot scholar "big data"
bigdat <- data.frame(years = 2000:2018,
   count = rev(c(75800, 95700, 87700, 79600, 52600, 31500, 12400, 4700, 2770, 2380,
      2070, 1900, 1930, 1210, 1160, 1010, 859, 703, 578)))

# Plot
p1 <- ggplot(data = bigdat, aes(y = count, x = (years)))
#p1 <- p1 + geom_smooth(se=F, method="loess", color="#004c93") #  Optional smoother
p1 <- p1 + geom_point(size=4,color="#004c93") #  Dots
p1 <- p1 + geom_line(size=1,color="#004c93") #  Line
p1 <- p1 + xlab("Year") + ylab("No. of Google Scholar Hits\n for 'Big Data'")
p1 <- p1 + scale_y_continuous(breaks = seq(0,90000,10000), labels = format(seq(0,90000,10000), scientific = FALSE)) # Breaks values
p1 <- p1 + scale_x_continuous(breaks = seq(2000,2018,2)) # Breaks years
p1   + theme_bw() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
    panel.spacing = unit(0.85, "lines"), text = element_text(size=35),
    legend.justification=c(0.9,0.05), legend.key = element_blank(), 
    legend.key.size = unit(5,"char"), 
    strip.background = element_rect(colour='white',fill = 'white')) # Look
ggsave("./results/Scholar_Bigdata.pdf", 
    height=11, width=15, units="in", device=cairo_pdf()) # Fig 1.2


### Plot Medline vs sc Abstracts

# Sociological abstracts:
# ti("record linkage")
# Only journals
# 1980-2019
# Export als XLS, Extraktion der Jahre mit Notepad++
socabstracts <- c(1980,
1986,
1987,
1988,
1991,
1996,
2000,
2000,
2000,
2005,
2005,
2006,
2006,
2006,
2007,
2007,
2011,
2011,
2013,
2013,
2013,
2014,
2014,
2016,
2017,
2017,
2018)
# Create table with counts
soc1 <- data.frame(years= names(table(socabstracts)), 
    count = as.numeric(table(socabstracts)))
# Create table with zeroes for all years
soc2 <- data.frame(years= factor(1980:2018), 
    count = rep(0,length(factor(1980:2018))))
# Fill zeroes with the years that have counts in table soc1
soc2[soc2$years %in% soc1$years,2] <- soc1[soc1$years %in% soc2$years,2]
rm(soc1) # remove temporary table
# Add source for plotting and make year numeric
soc2$source = "Sociological Abstracts"
soc2$years =  as.numeric(as.character(soc2$years))


# Medline:
#  "record linkage".m_titl. 
# Only journal articles
# 1980-2018 einzeln abgefragt
# Zahl je Jahr herausgeschrieben
medline <- tribble(
~years, ~count,
2019,4  ,
2018,63 ,
2017,78 ,
2016,59 ,
2015,65 ,
2014,52 ,
2013,35 ,
2012,47 ,
2011,46 ,
2010,31 ,
2009,30 ,
2008,24 ,
2007,14 ,
2006,22 ,
2005,24 ,
2004,12 ,
2003,8  ,
2002,12 ,
2001,13 ,
2000,13 ,
1999,13 ,
1998,20 ,
1997,15 ,
1996,7  ,
1995,17 ,
1994,10 ,
1993,6  ,
1992,5  ,
1991,9  ,
1990,11 ,
1989,17 ,
1988,5  ,
1987,7  ,
1986,11 ,
1985,11 ,
1984,6  ,
1983,4  ,
1982,1  ,
1981,6  ,
1980,5  
) 
# Add source for plotting
medline$source = "Medline"

# Make plottable data frame, delete 2018
plotdat <- subset(rbind(medline,soc2), as.numeric(years) < 2019)
plotdat <- melt(plotdat, id.vars = c("years","source"))

# Plot
p1 <- ggplot(data = plotdat, aes(y = value, x = (years), group=source, color=source))
#p1 <- p1 + geom_smooth(se=F, method="loess", color="#004c93") #  optional smoother
p1 <- p1 + geom_point(size=4) #  points
p1 <- p1 + geom_line(size=1) #  lines
p1 <- p1 + xlab("Year") + ylab("No. of Results for 'Record Linkage'\n") # labels
p1 <- p1 + scale_color_manual(name = "Database", labels=
    c("Medline","Sociological Abstracts"), 
    values=c("#004c93","#56B1F7")) # Set colors
p1 <- p1 + scale_x_continuous(breaks = seq(1980,2018,2)) # set year breaks
p1   + theme_bw() +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
    panel.spacing = unit(0.85, "lines"), text = element_text(size=35),
    legend.justification=c(0.9,0.05), legend.key = element_blank(), 
    legend.key.size = unit(5,"char"), 
    strip.background = element_rect(colour='white',fill = 'white')) # look
ggsave("./results/RL_years.pdf", 
    height=11, width=18, units="in", device=cairo_pdf()) # Fig 1.1
