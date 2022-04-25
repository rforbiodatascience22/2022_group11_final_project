tidy_data = read_csv("./data/data_tidy.csv")

data = tidy_data %>% 
  select(survival_days, alcohol_consumption, smoking, death_due_to_cancer) %>%
  drop_na()

data %>%
  ggplot(aes(x = survival_days, y = smoking, colour = death_due_to_cancer)) +
  geom_point()
