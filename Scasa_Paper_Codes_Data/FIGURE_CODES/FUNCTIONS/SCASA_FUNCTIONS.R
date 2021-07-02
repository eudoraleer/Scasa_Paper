## getting true counts for Scasa 
## beware of paralogs: either gene or isoforms
## input scasa = scasaAle_iso or scasaAle_gene
## output:  singleton genes/isoforms in true + paralogs in scasa,
##          NOT equalized dimension as input
##
scasa_truemat = function(scasa, tcount){
  para = rownames(scasa)
  paralist = strsplit(para,' ')   ##list of tx or genes
  npara = sapply(paralist, length)
  tgenes = rownames(tcount) ## genes in true mat
  # singletons observed in scasa
  para1 = para[npara==1]
    mat1 = tcount[(tgenes %in% para1),]
  # making paralogs
  mat2 = sapply(paralist[npara>1], parasum, count_isoform=tcount)
    colnames(mat2) = para[npara>1]
  # other true singles not observed in any form in scasa
  pick = !(rownames(tcount) %in% unlist(paralist))
  mat3 = NULL
  if (sum(pick)>0) mat3 = tcount[pick,]
  mat = rbind(mat1, t(mat2), mat3)
  return(mat)
}

##
parasum = function(tx, count_isoform) {
     ncell = ncol(count_isoform)
     pick = tx %in% rownames(count_isoform)
     if (sum(pick)>1) {return(colSums(count_isoform[tx[pick],]))}
     if (sum(pick)==1){return(count_isoform[tx[pick],])}
     if (sum(pick)==0){return(rep(0, ncell))}
 }

## convert isoform/paralog-level--> gene-level
## a <- scasa_genemat(scasaAle_iso, txmap=genes.tx.map.all.final)
##
scasa_genemat = function(scasa_iso, txmap)
{
  tx = rownames(scasa_iso)
  genes = txmap[tx]
  # rename the member genes outside the paralog to paralog-names
  gpara = unique(genes)
  glist = strsplit(gpara,' ')
  ngene = sapply(glist,length)
# member genes outside the paralogs
  gene_out = gpara[which(gpara %in% unlist(glist[ngene>1]))]
#check: 
#for (gene in gene_out){loc = grep(gene, gpara[ngene>1]);print(gpara[ngene>1][loc])}
#
  for (gene in gene_out){
#   loc = grep(gene, gpara[ngene>1]);
    loc = which(sapply(glist[ngene>1],function(x) gene%in%x ));
   genes[which(genes==gene)] = gpara[ngene>1][loc]
  }
## assigning transcripts
  txlist = tapply(tx, genes, c)  # list of tx/paralogs
  ntx = sapply(txlist, length)
  tx1 = unlist(txlist[ntx==1])
  mat1 = scasa_iso[tx1,]
    rownames(mat1) = names(tx1)
  mat2 = sapply(txlist[ntx>1], parasum, count_isoform=scasa_iso)
    colnames(mat2) = names(txlist[ntx>1])
  gcount = rbind(mat1, t(mat2))
    colnames(gcount) = colnames(scasa_iso)
  return(gcount)
}

#
# ..... expand APEs of gene paralogs to indiv genes
#
APE_extgene = function(APE){
  para = rownames(APE); 
  npara = sapply(strsplit(para,' '),length)
  para.rep = rep(para, npara)
  genes = unlist(strsplit(para,' '))   ##list of gene-level 
  out = APE[para.rep,]   ## repeat
    rownames(out) = genes
  return(out)
}

#
# ..... expand APEs of gene paralogs to indiv genes
# Old NOTE: is a gene appears both as single and in paralogs
# ... APE is averaged
APE_extgene_old = function(APE){
  para = rownames(APE); 
  npara = sapply(strsplit(para,' '),length)
  para.rep = rep(para, npara)
  genes = unlist(strsplit(para,' '))   ##list of gene-level 
  tmp1 = APE[para.rep,]   ## repeat
    rownames(tmp1) = genes
  ngenes= table(genes)
  sgenes = names(which(ngenes==1))  ## singles
  dgenes = names(which(ngenes>1))   ## multiples should exist
  tmp2 = NULL
  for (dgene in dgenes){
    a =tmp1[genes %in% dgene,]
    tmp2 = rbind(tmp2,colMeans(a))
  }
  rownames(tmp2) = dgenes
  colnames(tmp2) = colnames(APE)
  out = rbind(tmp1[sgenes,], tmp2)
  return(out)
}


##
## equalizing rownames of two matrices (not dataframe)
## output: each matrix with union of rownames
## values = 0 if original rowname is absent
## assume cols A is the same set of cols B
equrow = function(A,B){
  rowA = rownames(A)
  rowB = rownames(B)
  colA = colnames(A)
  uni = union(rowA, rowB)
  int = intersect(rowA,rowB)
  AminB = setdiff(rowA, rowB)
  BminA = setdiff(rowB,rowA)
  extA = matrix(0, nrow=length(BminA), ncol=ncol(A))
    rownames(extA) = BminA
  extB = matrix(0, nrow=length(AminB), ncol=ncol(B))
    rownames(extB) = AminB
  A = rbind(A, extA)[uni,colA]
  B = rbind(B, extB)[uni,colA]
  return(list(A=A, B=B))
}

##
## equalizing rows and cols of matlist = a list of matrices
## based on union of row and col names
## 
equrow_col = function(matlist){
  # find union of row and colnames
  ucols= Reduce(union, sapply(matlist,colnames))
  urows= Reduce(union, sapply(matlist,rownames))
  # padding each mat by zeros
  pad = function(mat){
    newcol = setdiff(ucols, colnames(mat))
    newrow = setdiff(urows, rownames(mat))
    xrow = matrix(0, nrow=length(newrow), ncol=ncol(mat))
      rownames(xrow) = newrow
    Mat = rbind(mat,xrow)
    xcol = matrix(0, ncol=length(newcol), nrow=nrow(Mat))
      colnames(xcol) = newcol
    Mat = cbind(Mat,xcol)
    Mat = Mat[urows,ucols]   ## equalize here
    return(Mat)
  }
  Matlist = lapply(matlist, pad)
  return(Matlist)
}


## ........    getting gene-level expression matrix: from Nghia
genemat = function(isoform_count,txmap)
{
SCASA=isoform_count
dim(SCASA)
pick=which(rowSums(SCASA)>0)
SCASA=SCASA[pick,]
dim(SCASA)

SCASA.libsize=colSums(SCASA)
summary(SCASA.libsize)

#######
allGenes=txmap[rownames(SCASA)]
isoform_count=SCASA
dim(isoform_count)

dup=duplicated(allGenes)
pick=which(allGenes %in% allGenes[dup])
gene_count1=isoform_count[-pick,]
rownames(gene_count1)=allGenes[-pick]

gene_count2=isoform_count[pick,]
allGenes=allGenes[pick]
uGenes=unique(allGenes)
gene_count3=matrix(0,nrow=length(uGenes),ncol=ncol(gene_count2))
cat("\n #of genes for processing: ",nrow(gene_count3))
for (i in 1:nrow(gene_count3)){
#	cat(" ",i)
	pick=which(allGenes %in% uGenes[i])
	gene_count3[i,]=colSums(gene_count2[pick,])
}
rownames(gene_count3)=uGenes
dim(gene_count3)

gene_count1=as.matrix(gene_count1)
gene_count=rbind(gene_count1,gene_count3)
gene_count=gene_count[,order(colnames(gene_count))]
dim(gene_count)

return(gene_count)
}



