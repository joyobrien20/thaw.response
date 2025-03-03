# Loop examples 

# Loop example
word <- c("c","o","f","f","e","e")
i < -3
word[i]
nextletter <- word[i + 1]
print(nextletter)
for (i in 1:3) { 
  nextletter <- (word[i + 1])
  print(nextletter)
}

# for each column in dorm_subtract: 
# for each element in column: 
# if element !=0: 
# add row ID to list. 
# remove all rows in row ID from sample in vector) 

badOTUs <- c()
for (i in 1:23) {
  for (j in 1:5725) {
    if (otu_table(dorm_subtract)[j,i] != 0) {
      badOTUs <- append(badOTUs,row.names(otu_table(dorm_subtract))[j])
    }
  }
}
dorm_taxaremoved <- subset_taxa(dorm1, !(rowname %in% badOTUs))