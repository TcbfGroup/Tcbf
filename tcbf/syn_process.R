library(GenomicRanges)
library(tidyverse)
library(feather)
library(arrow)
###本脚本用来处理将TAD边界比对到另一个物种的基因组上的共线性结果，
# 检查在另一个物种的基因组上比对区域是否也是TAD边界。
###



args <- commandArgs(trailingOnly  = TRUE)

syn <- args[1]

tad <- args[2]

bound <- args[3]

#output file name
network <- args[4]
conserved <- args[5]

minioverlap <- 3000


get_overlap <- function (query,subject){
  hits <- findOverlaps(query,subject,
                       minoverlap = minioverlaps)
  overlaps <- pintersect(query[queryHits(hits)],subject[subjectHits(hits)])
  tibble(as.data.frame(hits)) %>%
    mutate(width = width(overlaps))
}

# 读取基因组map的结果
synteny <-  arrow::read_feather(syn) %>%
  select(seq_id,seq_id2,start2,end2) %>%
  dplyr::rename(chromosome = seq_id2,
         start = start2,
         end = end2) %>%
  unique() %>%
  mutate(subjectHits = 1:dim(.)[1])

#读取reference genome的TAD 信息
target.tad <- read_csv(tad) %>%
  select(tad_name,chromosome,left_start,left_end,
         right_start,right_end) %>%
  mutate(left_bound =str_c(left_start,left_end,sep = "-"),
         right_bound =str_c(right_start,right_end,sep = "-") ) %>%
  select(-c(left_start,left_end,
            right_start,right_end)) %>%
  pivot_longer(cols = c(left_bound,right_bound),
               names_to = "bound",
               values_to = "region") %>%
  separate(region,c("start","end"),sep = "-") %>%
  mutate(start = as.integer(start),
         end = as.integer(end)) %>%
  mutate(queryHits  = 1:dim(.)[1])

target.bound <- read_tsv(bound)%>%
  mutate(queryHits  = 1:dim(.)[1])

D_bound.bed <- makeGRangesFromDataFrame(target.bound)
target.tad.bed <- makeGRangesFromDataFrame(target.tad)
synteny.bed <- makeGRangesFromDataFrame(synteny)

result1 <- tibble(as.data.frame(findOverlaps(target.tad.bed,
                                            synteny.bed,
                                            minoverlap = 3000)))

s1 <- select(target.tad,c(tad_name,bound,queryHits))
s2 <- select(synteny,c(seq_id,subjectHits))

final1 <- result1 %>%
  left_join(s1) %>%
  left_join(s2) %>%
  select(3:5) %>%
  mutate(bound = if_else(bound == "left_bound",0,1))


write_feather(final1,conserved)




minioverlaps <-  3000
syn <- arrow::read_feather(syn) %>%
  select(seq_id,seq_id2,start2,end2) %>%
  rename(chromosome = seq_id2,
         start = start2,
         end = end2) %>%
  unique()
bound <- read_tsv(bound) %>%
  mutate(subjectHits = 1:nrow(.))
bound.bed <- makeGRangesFromDataFrame(bound)


m <- split(syn,syn$seq_id)

fun <- function(x){
  name <- names(m)[x]
  x <- m[[x]]
  syn.bed <- GenomicRanges::reduce(makeGRangesFromDataFrame(x))
  res <- get_overlap(syn.bed,bound.bed) %>%

    group_by(subjectHits) %>%
    summarise(total = sum(width)) %>%
    dplyr::left_join(bound)%>%
    select(tad_name,total)
  res["query_tad"] <- name
  res[,c("tad_name","query_tad","total")]
}
library(data.table)

result <- rbindlist(lapply(seq_along(m),fun))
result <- result %>%
  as_tibble() %>%
  dplyr::select(query_tad,tad_name,total)

write_tsv(result,network)