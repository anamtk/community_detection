subsamples <- readRDS(here("05_visualizations",
                           "viz_data",
                           "ABUR1_bray_samples.RDS"))

summary <- readRDS(here("05_visualizations",
                        "viz_data",
                        "ABUR1_bray_summary.RDS"))

raw_bray <- readRDS(here("05_visualizations",
                       "viz_data",
                       "sbc_ABUR1_raw_bray.RDS"))

ggplot() +
  geom_line(data = subsamples, aes(x = year, y = bray, group = iter), color = "blue", alpha = 0.1) +
  geom_line(data = raw_bray, aes(x = year, y= raw_bray), color = "yellow") +
  geom_line(data = summary, aes(x = year, y = mean_bray), color = "darkblue") +
  theme_bw()
