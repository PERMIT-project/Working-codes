#otKD diagnostic codes for shorter running with larger datasets

#CKD IS STAGED HERE BASED ON MAINTENANCE OF EACH EGFR LEVEL FOR A MINIMUM OF 3 MONTHS
#CKD HAS BEEN DERIVED HERE FROM CKD-EPI EGFR- DECIDE IN ADVANCE WHICH EQUATION TO USE.

load("crea.rep.rda") #table featuring 1 row per patient per day with maximum creatinine retained

#Reduce dataset to individuals eligible for CKD and AKI staging before looping;
#Identify potential creatinine based CKD patients (individuals with any egfr<60)
#Isolate the first qualifying test for each patient

crea.rep$KDmark<-ifelse(!is.na(crea.rep$CKDEPIeGFR) & crea.rep$CKDEPIeGFR<60,1,0)
CKs<-crea.rep[crea.rep$KDmark==1,c("PatientID","event.date","KDmark")] #Kidney Injury flagged tests

#Create a dataset including all data from patients with at least 1 KD flag:
CKpot<-crea.rep[crea.rep$PatientID %in% CKs$PatientID,] 
#CKpot$EntryDate2<-as.Date(as.character(CKpot$event.date),format="%Y-%m-%d")

#CKpot<-CKpot[CKpot$entrydate>="20070101",]

#Mark a 90 day lookback period
CKpot$EntryDate1<-CKpot$event.date-90
CKpot<-CKpot[,c("PatientID","event.date","KDmark","EntryDate1","CKDEPIeGFR")] 

############################################################################################## CHECKED

#Diagnose CKD

#If eGFR under 60 is not sustained for 90 days, recode KDmark as 0
library("doParallel")
library("foreach")

cl <- makeCluster(10)
registerDoParallel(cl)

system.time(tmp <- foreach(x = 1:10, .combine=c("rbind"), .packages = "dplyr", .export = ls()) %dopar% {
  
  crea.rep[crea.rep$PatientID==CKpot$PatientID[x] & # select data within 90 days from the ith eventfor current patient
                crea.rep$event.date>CKpot$EntryDate1[x] & 
                  crea.rep$event.date<=CKpot$event.date[x],] %>%
        summarize(CKD=min(KDmark), MaxCKD=max(CKDEPIeGFR),MinCKD=min(CKDEPIeGFR)) %>%
          as.data.frame %>%
            bind_cols(CKpot[x, ])
 
})

stopCluster(cl)

#CKD is 1 if the test qualifies and there is no normal test within 3 months prior
#At this stage some are temporarily falsely positively identified that have no lookback test
#Near date match EntryDate1 (90 days prior date)- to entries from the full dataset
#neardate preferably matches to a prior entry if one is available
#if the closest match is before EntryDate1 (data from more than 90 days prior available), retain row.
indx1<-neardate(CKpot$PatientID,crea.rep$PatientID,CKpot$EntryDate1,crea.rep$event.date,best="prior",nomatch=NA_integer_)
CKpot$Lookback<-crea.rep[indx1,"event.date"]
CKpot$CKDGStage<-ifelse(CKpot$Lookback<=CKpot$EntryDate1 & CKpot$CKDE==1,1,0)
CKpot<-unique(CKpot[CKpot$CKDGStage==1,])
#CKpot is a subset table of CKD qualifying tests and their markers to be merged onto crea.rep

crea.rep<-merge(crea.rep,CKpot[,c(1,2,8:10)],all.x=TRUE)
crea.rep$CKDGStage<-ifelse(!is.na(crea.rep$CKDGStage),2,0)
#Stage 1 is skipped here as it cannot be identified from creatinine only, start by coding all as stage 2

crea.rep$CKDGStage<-ifelse(crea.rep$CKDGStage>0&crea.rep$MaxCKD>=30&crea.rep$MaxCKD<45,3.5,crea.rep$CKDGStage)
crea.rep$CKDGStage<-ifelse(crea.rep$CKDGStage>0&crea.rep$MaxCKD>=45&crea.rep$MaxCKD<60,3,crea.rep$CKDGStage)
crea.rep$CKDGStage<-ifelse(crea.rep$CKDGStage>0&crea.rep$MaxCKD>=15&crea.rep$MaxCKD<30,4,crea.rep$CKDGStage)
crea.rep$CKDGStage<-ifelse(crea.rep$CKDGStage>0&crea.rep$MaxCKD<15,5,crea.rep$CKDGStage)
crea.rep$CKDGStage<-ifelse(is.na(crea.rep$CKDGStage),0,crea.rep$CKDGStage)
table(crea.rep$CKDGStage)

crea.rep$CKDGStage<-ifelse(crea.rep$CKDGStage==3.5,paste("3b"),paste(crea.rep$CKDGStage))
crea.rep$MaxCKDGStage<-ifelse(crea.rep$MaxCKDGStage==3.5,paste("3b"),paste(crea.rep$MaxCKDGStage))

#########################################################################################################

#Mark start of first CKD diagnosis based on creatinine only

CKs<-CKpot %>%
group_by(PatientID)%>%
slice(which.min(EntryDate1)) %>%
as.data.frame
CKs<-CKs[,c(1,2)]
names(CKs)<-c("PatientID","CKDG_Date")

crea.rep<-merge(crea.rep,CKs,all.x=TRUE)
crea.rep$TimeSinceCKD<-difftime(as.Date(as.character(crea.rep$EntryDate),format="%Y%m%d"),crea.rep$CKDG_Date,unit="days")

########################################################################################################

#Incorporate Urine Albumin to creatinine ratio data if available:

crea.rep$CKDAStage<-ifelse(crea.rep$UACratio<3,1,0)
crea.rep$CKDAStage<-ifelse(crea.rep$UACratio>=3&crea.rep$UACratio<30&!crea.rep$CKDGStage=="0",2,crea.rep$CKDAStage)
crea.rep$CKDAStage<-ifelse(crea.rep$UACratio>30&!crea.rep$CKDGStage=="0",3,crea.rep$CKDAStage)


#Stage CKD based on both UAC and creatinine data:

crea.rep$CKDPrognosis<-crea.rep$CKDGStage

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage<=2 & crea.rep$CKDAStage==2,1,crea.rep$CKDPrognosis)

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage<=2 & crea.rep$CKDAStage==3,2,crea.rep$CKDPrognosis)

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage>=3,3,crea.rep$CKDPrognosis)

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage==3.5 & crea.rep$CKDAStage==1,2,crea.rep$CKDPrognosis)

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage==3 & crea.rep$CKDAStage==1,1,crea.rep$CKDPrognosis)

crea.rep$CKDPrognosis<-ifelse(crea.rep$CKDGStage==3 & crea.rep$CKDAStage==2,2,crea.rep$CKDPrognosis)



#Where CKD diagnosed based on creatinine, add metric for summarised eGFR range for further phenotyping

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MaxCKDE>=60 & crea.rep$MaxCKDE<90,1,NA)

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 & crea.rep$MaxCKDE>=60 & crea.rep$MaxCKDE<90,1,NA)

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MaxCKDE>=45 & crea.rep$MaxCKDE<60,2,crea.rep$CustomeGFR)

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MaxCKDE>=30 & crea.rep$MaxCKDE<45,3,crea.rep$CustomeGFR)

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MaxCKDE>=15 & crea.rep$MaxCKDE<30,4,crea.rep$CustomeGFR)

crea.rep$CustomeGFR<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MaxCKDE<15,5,crea.rep$CustomeGFR)



crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MinCKDE>=60 & crea.rep$MinCKDE<90,1,NA)

crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 & crea.rep$MinCKDE>=60 & crea.rep$MDRDeGFR<90,1,NA)

crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MinCKDE>=45 & crea.rep$MDRDeGFR<60,2,crea.rep$CustomeGFR)

crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MinCKDE>=30 & crea.rep$MDRDeGFR<45,3,crea.rep$CustomeGFR)

crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MinCKDE>=15 & crea.rep$MDRDeGFR<30,4,crea.rep$CustomeGFR)

crea.rep$CustomeGFR2<-ifelse(crea.rep$CKDPrognosis>0 &crea.rep$MinCKDE<15,5,crea.rep$CustomeGFR)
