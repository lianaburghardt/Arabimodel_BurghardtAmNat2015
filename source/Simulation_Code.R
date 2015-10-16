###### Here is the code to generate all the results found in the Dryad results folder so you can see how they were made #########

###### Source in the functions from Arabi_Base_final ##########
#    setwd("/Users/X/X/")
#    source("./Arabi_Base_final.r")

####################################VALENCIA########################################
env<-read.csv(file="./environmental_data_Valencia_60yrs.csv",header=TRUE)

##### low FRI/ Low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000                   
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######High FRI/ low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/high dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

###### High Fri/High dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                          
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

#####High FRI/ mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                                         
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

rm(env)

################################HALLE#######################################
env<-read.csv(file="./environmental_data_Halle_60yrs.csv.csv",header=TRUE)

#####low FRI/ low Dor
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000                 
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

####### high FRI/ Low Dor
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/high dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

###### High Fri/High dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                          
                               )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######High FRI/ mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                                         
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

rm(env)

################################# OULU################################
env<-read.csv(file="./environmental_data_Oulu_60yrs.csv",header=TRUE)

######Low Fri/low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                	moisture=env[,"moisture"],
                               	day_length=env[,"day_length"],
                               	day_hours=env[,"day_hours"],
                                 	max.cohort.store=10000                 
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######High Fri/low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                time.step=1,
                                temp=env[,"temperature"],
                                moisture=env[,"moisture"],
                                day_length=env[,"day_length"],
                                day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                              , Fi=.737                      
                               )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/high dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

###### High Fri/High dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                              , Fi=.737
                                , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                          
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######High FRI/ mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                                         
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

#####Low Fri/mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

rm(env)

######################NORWICH################################
env<-read.csv(file="./environmental_data_Valencia_60yrs.csv",header=TRUE)

#####Low Fri/Low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000                  
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######High Fri/low dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

######Low Fri/high dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

###### High Fri/High dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=2.5  # It takes 100 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                          
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

#####High FRI/ mid dormancy

rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , Fi=.737
                                , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                                         
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

#####Low Fri/mid dormancy
rc <- iteratePop(Tmax=(24*365*60),
                                 time.step=1,
                                 temp=env[,"temperature"],
                                 moisture=env[,"moisture"],
                                 day_length=env[,"day_length"],
                                 day_hours=env[,"day_hours"],
                                 max.cohort.store=10000
                               , psi.mean=1.25  # It takes 50 days of after-ripening at 22C for half of seeds to germinate at 22 like cold matured Col                      
                                )

life<-whatistheLifecycle.output(rc=rc,first.year=15,last.year=60)     ######## data summarization function set to summarize the last 45 years of output- save this and/or the rc file for working with...
rm(rc)
rm(life)

rm(env)
