#!/usr/bin/env Rscript

library(data.table)

dt <- fread("merged_CHG_symmetric_all.txt.gz")
chromatin_states <- fread("chromatin_states.txt",header = FALSE)
out_suffix <- "_chg"

colnames(chromatin_states) <- c("chr","start","end","chromatin_state")
chromatin_states[, chromatin_state := as.character(chromatin_state)]

dt <- dt[cov_CS > 10 & cov_CSxP > 10 & cov_P > 10]

dt[, subgenome := sub("^chr[0-9]+([ABD])$","\\1",chr)]
dt <- dt[subgenome %in% c("A","B","D")]

dt[, start_pos := pos]
dt[, end_pos := pos]

setkey(chromatin_states,chr,start,end)
setkey(dt,chr,start_pos,end_pos)

dt <- foverlaps(dt,chromatin_states,by.x = c("chr","start_pos","end_pos"),by.y = c("chr","start","end"),type = "within",nomatch = 0)

dt <- dt[chromatin_state %in% c("1-4","13")]

dt[, A := pct_CS / 100]
dt[, H := pct_CSxP / 100]
dt[, B := pct_P / 100]

dt[, x := H - A]
dt[, y := H - B]
dt[, radius := sqrt(x^2 + y^2)]
dt[, angle_deg := (atan2(y,x) * 180 / pi) %% 360]

circ_dist <- function(a,b) {
  d <- abs(a - b)
  pmin(d,360 - d)
}

nearest_sector_center <- function(angle_deg) {
  centers <- c(0,45,90,135,180,225,270,315)
  centers[which.min(circ_dist(angle_deg,centers))]
}

sector_to_class <- function(sector_center) {
  if (sector_center %in% c(0,180)) {
    "P_dominant"
  } else if (sector_center %in% c(90,270)) {
    "CS_dominant"
  } else if (sector_center %in% c(135,315)) {
    "additive"
  } else if (sector_center == 45) {
    "overdominant"
  } else if (sector_center == 225) {
    "underdominant"
  } else {
    NA_character_
  }
}

dt[, sector_center := vapply(angle_deg,nearest_sector_center,numeric(1))]
dt[, max_mC := pmax(A,H,B,na.rm = TRUE)]
dt <- dt[max_mC >= 0.01]

dt[, category := ifelse(radius < 0.1,"conserved_mC",vapply(sector_center,sector_to_class,character(1)))]
dt[, category := factor(category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))]

summary_dt <- dt[, .N, by = .(chromatin_state,category)]
summary_dt[, proportion := N / sum(N), by = chromatin_state]
summary_dt[, chromatin_state := factor(chromatin_state,levels = c("1-4","13"))]
summary_dt[, category := factor(category,levels = c("conserved_mC","additive","CS_dominant","P_dominant","overdominant","underdominant"))]
setorder(summary_dt,chromatin_state,category)

write.table(summary_dt,file = paste0("percent_met_chromatin_state",out_suffix,".tsv"),sep = "\t",quote = FALSE,row.names = FALSE)
