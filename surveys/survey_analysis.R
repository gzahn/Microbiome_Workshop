library(tidyverse)
library(readxl)
library(patchwork)
library(wordcloud)

theme_set(theme_minimal() +
            theme(axis.text.x = element_text(angle=90,hjust=1)))

dat <- read_xlsx("./surveys/survey_data.xlsx")
dat_full <- dat
dat <- dat %>% filter(pre_post == 'pre')

# basic info about participants ####
p1 <- dat %>% 
  ggplot(aes(x=fields)) +
  geom_bar()

p2 <- dat %>% 
  ggplot(aes(x=career_stage)) +
  geom_bar()

p3 <- dat %>% 
  ggplot(aes(x=operating_system)) +
  geom_bar()

p4 <- dat %>% 
  select(gui_stats_use, coding_use, databse_use,version_control_use,command_line_use) %>% 
  pivot_longer(1:5) %>% 
  mutate(value=factor(value,levels = c("never","less than once per year","several times per year","monthly","weekly"))) %>% 
  ggplot(aes(x=value)) +
  geom_bar() +
  facet_wrap(~name)

(p1 + p2) / (p3 + p4)
ggsave("./surveys/survey_participants.png",height = 8,width = 6)

# pre- vs post-workshop values
dat_full %>% 
  select(participant_id, pre_post, starts_with("q")) %>% 
  pivot_longer(starts_with("q")) %>% 
  mutate(pre_post=factor(pre_post,levels=c("pre","post"))) %>% 
  ggplot(aes(x=pre_post,y=value)) +
  geom_point() +
  geom_line(aes(group=participant_id)) +
  facet_wrap(~name)
ggsave("./surveys/pre-post_attitudes.png")

# free response info ####
words <- dat_full$notes %>% 
  str_replace_all("/"," ") %>% 
  str_split(" ") %>% 
  unlist() %>% 
  str_remove_all(";|,|\\.") %>% 
  table() %>% 
  as.data.frame()
names(words) <- c("word","freq")
words <- words %>% 
  mutate(word = as.character(word))

words <- words %>% 
  filter(nchar(word) > 3)

wordcloud(words = words$word, freq = words$freq, min.freq = 1,           
          max.words=200, random.order=FALSE, rot.per=0.35,            
          colors=brewer.pal(8, "Dark2"),
          scale = c(3,.25))
ggsave("./surveys/wordcloud_free_response.png")
