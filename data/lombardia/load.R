dat <- readRDS("data_lombardia.rds")

ggplot(dat$lombardy_data_4may) +
  geom_col(aes(x = date, y = cases), fill = "seagreen") +
  geom_col(aes(x = date, y = deaths), fill = "firebrick")
