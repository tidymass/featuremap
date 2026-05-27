library(ggplot2)
library(dplyr)

setwd("~/Documents/AtNTUlife/ShenLab/1_MyProject/featureFMA/featuremap_case/4_manuscript")

# Simulated HMDB class composition (replace counts with real data if you have it)
donut <- tibble(
  class = c("Lipids", "Organic acids", "Organoheterocyclics",
            "Benzenoids", "Other"),
  count = c(38, 24, 18, 12, 28)
) %>%
  mutate(
    class = factor(class, levels = class),          # keep order
    frac  = count / sum(count),
    ymax  = cumsum(frac),
    ymin  = lag(ymax, default = 0)
  )

# Colours: four classes + grey for "Other" (matches your reference)
donut_cols <- c("Lipids"             = "#2c6fbb",   # blue
                "Organic acids"      = "#9c7a1e",   # yellow-brown
                "Organoheterocyclics"= "#e8821e",   # orange
                "Benzenoids"         = "#7a3fb0",   # purple
                "Other"              = "#9e9e9e")   # grey

p <- ggplot(donut, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = class)) +
  geom_rect(colour = "white", linewidth = 1.2) +
  coord_polar(theta = "y") +
  xlim(c(1.5, 4)) +                                  # inner radius -> donut hole
  scale_fill_manual(values = donut_cols) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text  = element_text(size = 11),
        plot.background = element_rect(fill = "white", colour = NA))

pdf("./Figures/Figure1/class_profile.pdf", width = 4, height = 4)
p
dev.off()
