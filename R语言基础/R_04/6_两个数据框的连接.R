test1 <- data.frame(name = c('jimmy','nicker','Damon','Sophie'), 
                    blood_type = c("A","B","O","AB"))
test1
test2 <- data.frame(name = c('Damon','jimmy','nicker','tony'),
                    group = c("group1","group1","group2","group2"),
                    vision = c(4.2,4.3,4.9,4.5))
test2
library(dplyr)
inner_join(test1,test2,by="name")
right_join(test1,test2,by="name")
full_join(test1,test2,by="name")
semi_join(test1,test2,by="name")
anti_join(test1,test2,by="name")
